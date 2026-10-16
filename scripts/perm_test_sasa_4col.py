#!/usr/bin/env python3
"""
Permutation test: are modified rRNA sites more solvent-exposed than other sites
of the same base?

Input CSV: four columns of SASA values, in this order (columns may be different
lengths; leave the unused cells blank):

    1. 28S modified sites of interest
    2. 28S all other sites of that base
    3. 18S modified sites of interest
    4. 18S all other sites of that base

The header names can be anything. The first row must be a header row.
If your columns are in a different order, pass --cols with the four column
names in the order: 28S_modified 28S_other 18S_modified 18S_other.

Because each background already contains only the modifiable base, the random
draws are automatically base-matched.

Two separate tests are run: one for 28S and one for 18S, each against its own
background.

Usage:
    python perm_test_sasa_4col.py sasa.csv --n-perm 100000 --stat median --out results
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.stats import mannwhitneyu
except ImportError:  # scipy is optional; only used for the cross-check columns
    mannwhitneyu = None


def load_columns(path, cols):
    df = pd.read_csv(path)
    if cols is None:
        if df.shape[1] < 4:
            raise SystemExit(f"Expected 4 columns, found {df.shape[1]}.")
        cols = list(df.columns[:4])
        print("Using columns in this order:")
        for label, c in zip(["28S modified", "28S other", "18S modified", "18S other"], cols):
            print(f"  {label:13s} <- '{c}'")
    else:
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise SystemExit(f"Column(s) not found in CSV: {missing}. Found: {list(df.columns)}")

    out = []
    for c in cols:
        raw = df[c]
        num = pd.to_numeric(raw, errors="coerce")
        bad = int(raw.notna().sum() - num.notna().sum())
        if bad:
            print(f"WARNING: {bad} non-numeric value(s) in column '{c}' were ignored.")
        out.append(num.dropna().to_numpy(dtype=float))
    return dict(zip(["m28", "b28", "m18", "b18"], out))


def perm_test(strata, n_perm, stat_fn, alternative, rng):
    """strata: list of (modified_values, background_values). Returns obs, null, p."""
    for mod, bg in strata:
        if len(mod) == 0:
            raise SystemExit("A modified-sites column is empty.")
        if len(bg) < len(mod):
            raise SystemExit("A background column has fewer values than its modified-sites column.")

    obs = stat_fn(np.concatenate([m for m, _ in strata]))
    null = np.empty(n_perm)
    for i in range(n_perm):
        draw = np.concatenate([rng.choice(bg, size=len(mod), replace=False)
                               for mod, bg in strata])
        null[i] = stat_fn(draw)

    p_greater = (1 + np.sum(null >= obs)) / (1 + n_perm)
    p_less = (1 + np.sum(null <= obs)) / (1 + n_perm)
    if alternative == "greater":
        p = p_greater
    elif alternative == "less":
        p = p_less
    else:
        p = min(1.0, 2 * min(p_greater, p_less))
    return obs, null, p


def mwu_crosscheck(strata, alternative):
    """Mann-Whitney U (modified vs. all other sites) and AUC, as a sanity check."""
    if mannwhitneyu is None:
        return np.nan, np.nan
    mod = np.concatenate([m for m, _ in strata])
    bg = np.concatenate([b for _, b in strata])
    res = mannwhitneyu(mod, bg, alternative=alternative)
    auc = res.statistic / (len(mod) * len(bg))
    return res.pvalue, auc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("--cols", nargs=4, metavar=("MOD28", "OTHER28", "MOD18", "OTHER18"),
                    help="column names, in the order: 28S_modified 28S_other 18S_modified 18S_other")
    ap.add_argument("--n-perm", type=int, default=100_000)
    ap.add_argument("--stat", choices=["mean", "median"], default="median")
    ap.add_argument("--alternative", choices=["greater", "less", "two-sided"], default="greater",
                    help="'greater' tests whether modified sites are MORE exposed (default)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="perm_results")
    args = ap.parse_args()

    d = load_columns(args.csv, args.cols)
    stat_fn = np.median if args.stat == "median" else np.mean
    rng = np.random.default_rng(args.seed)

    sets = {
        "28S": [(d["m28"], d["b28"])],
        "18S": [(d["m18"], d["b18"])],
    }

    results = {name: perm_test(s, args.n_perm, stat_fn, args.alternative, rng)
               for name, s in sets.items()}

    # ---- summary table ----
    rows = []
    for name, (obs, null, p) in results.items():
        mw_p, auc = mwu_crosscheck(sets[name], args.alternative)
        rows.append({
            "set": name,
            "n_modified": int(sum(len(m) for m, _ in sets[name])),
            "n_background": int(sum(len(b) for _, b in sets[name])),
            f"observed_{args.stat}_SASA": obs,
            f"null_{args.stat}": np.median(null),
            "null_95pct_upper": np.percentile(null, 95),
            "perm_p": p,
            "MannWhitney_p": mw_p,
            "AUC": auc,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(f"{args.out}.csv", index=False)
    print(f"\nPermutation test ({args.n_perm:,} permutations, alternative = {args.alternative}):\n")
    print(summary.to_string(index=False))

    # ---- plot ----
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, (name, (obs, null, p)) in zip(axes, results.items()):
        ax.hist(null, bins=60, color="0.7", edgecolor="white")
        ax.axvline(obs, color="crimson", lw=2, label=f"Observed = {obs:.2f}")
        n_mod = sum(len(m) for m, _ in sets[name])
        ax.set_title(f"{name}  (n = {n_mod}, p = {p:.4g})")
        ax.set_xlabel(f"{args.stat.capitalize()} SASA of random site sets")
        ax.set_ylabel("Permutations")
        ax.legend(frameon=False, loc="upper left")
    fig.suptitle(f"Null distribution from {args.n_perm:,} random draws of same-base sites")
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=300)
    fig.savefig(f"{args.out}.pdf")
    print(f"\nSaved {args.out}.csv, {args.out}.png, {args.out}.pdf")


if __name__ == "__main__":
    main()
