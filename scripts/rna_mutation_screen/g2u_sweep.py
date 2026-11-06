#!/usr/bin/env python3
"""
g2u_sweep.py

Sweep every G->U point mutation in an RNA sequence and score its effect on a
FIXED reference secondary structure (rationale: see accompanying README).

Method
------
For a fixed reference structure S and sequence X:
    dG(X, S)  = ViennaRNA structure-evaluation energy (equivalent to RNAeval)
This script computes dG(WT, S) once, then for every position i where WT has a
'G', builds a mutant sequence with a single G->U substitution and computes
dG(mutant, S). It reports:
    ddG_i = dG(mutant_i, S) - dG(WT, S)

Higher (more positive) ddG = more destabilizing relative to the fixed
reference fold. Positions that are unpaired in S will typically have ddG
near zero and serve as an internal sanity check.

Requires the ViennaRNA python bindings (`pip install ViennaRNA` /
`conda install -c bioconda viennarna`). See README.md for full setup.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path

try:
    import RNA
except ImportError as e:
    sys.exit(
        "ERROR: could not import the ViennaRNA python bindings ('RNA' module).\n"
        "Install with:  pip install ViennaRNA   (or)   conda install -c bioconda viennarna\n"
        f"Original error: {e}"
    )

VALID_BASES = set("ACGU")


# --------------------------------------------------------------------------
# Input parsing + validation
# --------------------------------------------------------------------------

@dataclass
class LoadedInputs:
    name: str
    sequence: str
    structure: str
    pair_map: dict[int, int]  # 0-based position -> 0-based partner position


def _fail(msg: str) -> None:
    sys.exit(f"INPUT VALIDATION ERROR: {msg}")


def read_fasta(path: Path) -> tuple[str, str]:
    """Read a single-record FASTA file. Returns (header, sequence)."""
    text = path.read_text().strip().splitlines()
    if not text:
        _fail(f"{path} is empty.")
    if not text[0].startswith(">"):
        _fail(f"{path} does not start with a FASTA header line ('>...').")
    header = text[0][1:].strip()
    seq_lines = [ln.strip() for ln in text[1:] if ln.strip()]
    if not seq_lines:
        _fail(f"{path} has a header but no sequence.")
    # reject a second header line (multi-record FASTA not supported here)
    for ln in seq_lines:
        if ln.startswith(">"):
            _fail(
                f"{path} appears to contain more than one FASTA record. "
                "This script expects a single sequence per file."
            )
    sequence = "".join(seq_lines).upper()
    return header, sequence


def normalize_sequence(seq: str, source: str) -> str:
    """Upper-case, convert DNA 'T'->'U', and validate alphabet."""
    seq = seq.upper()
    if "T" in seq:
        n_t = seq.count("T")
        print(
            f"NOTE: {source} contains {n_t} 'T' character(s); "
            "converting to 'U' (treating as RNA).",
            file=sys.stderr,
        )
        seq = seq.replace("T", "U")
    bad_chars = sorted(set(seq) - VALID_BASES)
    if bad_chars:
        _fail(
            f"{source} contains non-ACGU character(s) after normalization: {bad_chars}. "
            "Ambiguity codes (N, R, Y, ...) and gaps are not supported by this script."
        )
    return seq


def read_dotbracket(path: Path) -> str:
    """Read a dot-bracket structure file.

    Accepts either:
      - a bare single line of './(' / ')' characters (what RNAcentral's
        "copy dot-bracket" button gives you), or
      - a two-line FASTA-like file (header + structure), which is also
        tolerated in case you paste it with a header by habit.
    """
    lines = [ln.strip() for ln in path.read_text().strip().splitlines() if ln.strip()]
    if not lines:
        _fail(f"{path} is empty.")
    # drop a leading '>' header line if present
    if lines[0].startswith(">"):
        lines = lines[1:]
    if not lines:
        _fail(f"{path} only contained a header line, no structure.")
    if len(lines) > 1:
        _fail(
            f"{path} has more than one non-header line. Expected a single "
            "dot-bracket string on one line."
        )
    structure = lines[0]
    allowed = set(".()")
    bad_chars = sorted(set(structure) - allowed)
    if bad_chars:
        _fail(
            f"{path} contains character(s) not in '.()': {bad_chars}. "
            "Pseudoknot bracket types ([], {{}}, etc.) and probability "
            "annotations are not supported -- this pipeline requires a "
            "pseudoknot-free nested structure."
        )
    return structure


def build_pair_map(structure: str, source: str) -> dict[int, int]:
    """Validate proper nesting of a dot-bracket string and return a 0-based
    position -> partner position map (unpaired positions are absent)."""
    stack: list[int] = []
    pair_map: dict[int, int] = {}
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                _fail(
                    f"{source}: unbalanced structure -- ')' at position "
                    f"{i + 1} has no matching '('."
                )
            j = stack.pop()
            pair_map[i] = j
            pair_map[j] = i
    if stack:
        unmatched = [p + 1 for p in stack]
        _fail(
            f"{source}: unbalanced structure -- unmatched '(' at position(s) "
            f"{unmatched[:5]}{'...' if len(unmatched) > 5 else ''}."
        )
    return pair_map


def load_and_validate(name: str, fasta_path: Path, dbn_path: Path) -> LoadedInputs:
    header, raw_seq = read_fasta(fasta_path)
    sequence = normalize_sequence(raw_seq, str(fasta_path))
    structure = read_dotbracket(dbn_path)

    if len(sequence) != len(structure):
        _fail(
            f"Length mismatch for '{name}': sequence is {len(sequence)} nt "
            f"({fasta_path}) but structure is {len(structure)} characters "
            f"({dbn_path}). These must be exactly the same length and "
            "positionally aligned."
        )

    pair_map = build_pair_map(structure, str(dbn_path))

    print(f"[{name}] loaded OK: {len(sequence)} nt, {len(pair_map) // 2} base pairs.")
    print(f"[{name}] FASTA header: '{header}'")
    return LoadedInputs(name=name, sequence=sequence, structure=structure, pair_map=pair_map)


# --------------------------------------------------------------------------
# Energetics
# --------------------------------------------------------------------------

def eval_energy(sequence: str, structure: str) -> float:
    """dG (kcal/mol) of `sequence` folded into the fixed `structure`.

    Uses RNA.OPTION_EVAL_ONLY: this skips allocating the O(n^2) MFE/partition
    function DP matrices that fold_compound() builds by default (which we
    never use here, since we only ever call eval_structure). Confirmed to
    give bit-identical energies to the default constructor; the only
    difference is allocation. This matters a lot for long sequences: for
    a ~5,000 nt rRNA, this drops the per-call cost from ~230ms to ~0.25ms
    (~1000x), which is the difference between a sweep finishing in seconds
    vs. taking the better part of ten minutes.
    """
    fc = RNA.fold_compound(sequence, None, RNA.OPTION_EVAL_ONLY)
    return fc.eval_structure(structure)


def run_sweep(
    data: LoadedInputs,
    from_base: str = "G",
    to_base: str = "U",
) -> tuple[float, list[dict]]:
    seq = data.sequence
    struct = data.structure
    pair_map = data.pair_map

    dg_wt = eval_energy(seq, struct)
    print(f"[{data.name}] reference dG (WT, full-length): {dg_wt:.2f} kcal/mol")

    positions = [i for i, b in enumerate(seq) if b == from_base]
    print(f"[{data.name}] found {len(positions)} '{from_base}' positions to mutate.")

    rows = []
    for i in positions:
        mutant = seq[:i] + to_base + seq[i + 1 :]
        dg_mut = eval_energy(mutant, struct)
        ddg = dg_mut - dg_wt
        partner = pair_map.get(i)
        rows.append(
            {
                "position": i + 1,  # 1-based for human readability
                "paired": partner is not None,
                "partner_position": (partner + 1) if partner is not None else "",
                "partner_base": seq[partner] if partner is not None else "",
                "ddG_kcal_per_mol": round(ddg, 4),
            }
        )

    return dg_wt, rows


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def sha256sum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(
    out_path: Path,
    data: LoadedInputs,
    dg_wt: float,
    rows: list[dict],
    from_base: str,
    to_base: str,
    fasta_path: Path,
    dbn_path: Path,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        f.write(f"# rna: {data.name}\n")
        f.write(f"# mutation: {from_base} -> {to_base} (single substitutions only)\n")
        f.write(f"# reference_dG_wt_full_length_kcal_per_mol: {dg_wt:.4f}\n")
        f.write(
            "#   (ddG below = dG(mutant, reference_structure) - reference_dG_wt; "
            "add reference_dG_wt back to recover the mutant's absolute dG)\n"
        )
        f.write(f"# sequence_length_nt: {len(data.sequence)}\n")
        f.write(f"# reference_base_pairs: {len(data.pair_map) // 2}\n")
        f.write(f"# fasta_source: {fasta_path.name} (sha256 {sha256sum(fasta_path)})\n")
        f.write(f"# structure_source: {dbn_path.name} (sha256 {sha256sum(dbn_path)})\n")
        f.write(f"# ViennaRNA_version: {RNA.__version__}\n")
        f.write(f"# generated: {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        writer = csv.DictWriter(
            f,
            fieldnames=["position", "paired", "partner_position", "partner_base", "ddG_kcal_per_mol"],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"[{data.name}] wrote {len(rows)} rows to {out_path}")


def print_summary(data: LoadedInputs, rows: list[dict]) -> None:
    paired = [r for r in rows if r["paired"]]
    unpaired = [r for r in rows if not r["paired"]]

    def stats(label: str, subset: list[dict]) -> None:
        if not subset:
            print(f"  {label}: (none)")
            return
        vals = [r["ddG_kcal_per_mol"] for r in subset]
        print(
            f"  {label}: n={len(vals)}  mean={sum(vals)/len(vals):.2f}  "
            f"min={min(vals):.2f}  max={max(vals):.2f} kcal/mol"
        )

    print(f"\n[{data.name}] summary:")
    stats("paired positions  ", paired)
    stats("unpaired positions", unpaired)

    top5 = sorted(rows, key=lambda r: r["ddG_kcal_per_mol"], reverse=True)[:5]
    print(f"  top 5 most destabilizing mutations:")
    for r in top5:
        partner = f"(pairs with {r['partner_base']}{r['partner_position']})" if r["paired"] else "(unpaired)"
        print(f"    pos {r['position']:>5}  ddG={r['ddG_kcal_per_mol']:>7.2f}  {partner}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, help="Short label for this RNA, e.g. '18S'")
    ap.add_argument("--fasta", required=True, type=Path, help="Path to single-record FASTA file")
    ap.add_argument("--structure", required=True, type=Path, help="Path to dot-bracket structure file")
    ap.add_argument("--out", required=True, type=Path, help="Output CSV path")
    ap.add_argument("--from-base", default="G", help="Wild-type base to mutate from (default: G)")
    ap.add_argument("--to-base", default="U", help="Mutant base to mutate to (default: U)")
    args = ap.parse_args()

    for p in (args.fasta, args.structure):
        if not p.exists():
            _fail(f"file not found: {p}")

    data = load_and_validate(args.name, args.fasta, args.structure)
    dg_wt, rows = run_sweep(data, from_base=args.from_base.upper(), to_base=args.to_base.upper())
    write_csv(args.out, data, dg_wt, rows, args.from_base.upper(), args.to_base.upper(), args.fasta, args.structure)
    print_summary(data, rows)


if __name__ == "__main__":
    main()
