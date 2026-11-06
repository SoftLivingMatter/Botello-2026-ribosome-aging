# /// script
# dependencies = [
#     "marimo",
#     "pandas==3.0.5",
#     "openpyxl",
#     "scipy==1.18.1",
#     "statsmodels",
#     "plotly==7.0.0",
#     "numpy==2.5.2",
# ]
# requires-python = ">=3.14"
# ///

"""
Age-dependent protein interactions from RPL10A and RPS17 pulldown proteomics.

Marimo notebook. For each protein, a one-way ANOVA across the T0, T24 and T72
time points is run separately for each pulldown, p-values are corrected by
Benjamini-Hochberg within each pulldown, and proteins are classified as
increasing (+), decreasing (-) or non-monotonic (0) from the group means.

Expects L10Aproteomics_08182026.xlsx and S17proteomics_08182026.xlsx in the same
directory as this file.  Writes results.csv to the same directory.

    marimo edit anova_trend.py      # interactive
    python anova_trend.py           # run as a script
"""

import marimo

__generated_with = "0.23.6"
app = marimo.App(width="medium")


@app.cell
def _():
    from pathlib import Path

    import marimo as mo
    import numpy as np
    import pandas as pd
    from scipy import stats

    DATA_DIR = Path(mo.notebook_dir() or ".")
    return DATA_DIR, np, pd, stats


@app.cell
def _(np, pd, stats):
    def run_anova(df, group_cols, labels=None, experiment_col="Experiment", gene_col="Gene Name"):
        from statsmodels.stats.multitest import multipletests
        if labels is None:
            labels = [f"G{i}" for i in range(len(group_cols))]
        records = []
        for (experiment, gene), grp in df.groupby([experiment_col, gene_col]):
            groups = []
            for cols in group_cols:
                vals = grp[cols].values.flatten()
                vals = vals[~np.isnan(vals)]
                groups.append(vals)

            if any(len(g) < 2 for g in groups):
                continue

            _, p = stats.f_oneway(*groups)
            means = [g.mean() for g in groups]

            if   all(a < b for a, b in zip(means, means[1:])): trend = "+"
            elif all(a > b for a, b in zip(means, means[1:])): trend = "-"
            else:                                               trend = "0"

            record = {"Experiment": experiment, "Gene": gene}
            for label, mean in zip(labels, means):
                record[f"Mean {label}"] = round(mean, 6)
            record["p-value"] = p
            record["Trend"]   = trend
            records.append(record)

        _df = pd.DataFrame(records).set_index(["Experiment", "Gene"]).dropna(subset=["p-value"])
        _df["q-value"] = np.nan
        for _exp in _df.index.get_level_values("Experiment").unique():
            _mask = _df.index.get_level_values("Experiment") == _exp
            _pvals = _df.loc[_mask, "p-value"]
            _, _q, _, _ = multipletests(_pvals, method="fdr_bh")
            _df.loc[_mask, "q-value"] = pd.Series(_q, index=_pvals.index)
        return _df.sort_values("p-value")

    return (run_anova,)


@app.cell
def _(DATA_DIR, pd):
    df_l10a = pd.read_excel(
        DATA_DIR / "L10Aproteomics_08182026.xlsx",
        index_col=0,
    )
    df_l10a
    return (df_l10a,)


@app.cell
def _(DATA_DIR, pd):
    df_s17 = pd.read_excel(
        DATA_DIR / "S17proteomics_08182026.xlsx",
        index_col=0,
    )
    df_s17
    return (df_s17,)


@app.cell
def _(df_l10a, df_s17, pd):
    df_merged = pd.concat([df_l10a.assign(Experiment="L10A"), df_s17.assign(Experiment="S17")])
    df_merged
    return (df_merged,)


@app.cell
def _(df_merged, run_anova):
    LOW_COLS  = ["T0-1",  "T0-2",  "T0-3"]
    MID_COLS  = ["T24-1", "T24-2", "T24-3"]
    HIGH_COLS = ["T72-1", "T72-2", "T72-3"]

    df_results = run_anova(df_merged, [LOW_COLS, MID_COLS, HIGH_COLS], labels=["T0", "T24", "T72"])
    df_results
    return (df_results,)


@app.cell
def _(df_results, df_top10, np):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    _COLORS = {
        "+": {"light": "#2a78d6", "dark": "#3987e5"},
        "-": {"light": "#eb6834", "dark": "#d95926"},
    }
    _LABELS = {"+": "Increasing (+)", "-": "Decreasing (−)"}

    _df_plot = df_results[
        (df_results["Trend"] != "0") &
        (df_results["Mean T0"] > 0) &
        (df_results["Mean T72"] > 0)
    ].copy()
    _df_plot["log2FC"]    = np.log2(_df_plot["Mean T72"] / _df_plot["Mean T0"])
    _df_plot["neg_log_p"] = -np.log10(_df_plot["p-value"])

    def _q_threshold_y(df_exp):
        sig = df_exp[df_exp["q-value"] <= 0.05]
        if sig.empty:
            return None
        return -np.log10(sig["p-value"].max())

    def _volcano_traces(df_exp, showlegend):
        traces = []
        for trend in ["+", "-"]:
            sub = df_exp[df_exp["Trend"] == trend]
            if sub.empty:
                continue
            traces.append(go.Scatter(
                x=sub["log2FC"],
                y=sub["neg_log_p"],
                mode="markers",
                name=_LABELS[trend],
                legendgroup=trend,
                showlegend=showlegend,
                marker=dict(color=_COLORS[trend]["light"], size=5,
                            opacity=0.65, line=dict(width=0)),
                customdata=np.stack([
                    sub.index.get_level_values("Gene"),
                    sub["Mean T0"].values, sub["Mean T24"].values,
                    sub["Mean T72"].values, sub["p-value"].values,
                    sub["q-value"].values,
                ], axis=1),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "log₂FC: %{x:.2f}<br>−log₁₀(p): %{y:.2f}<br>"
                    "Mean T0: %{customdata[1]:.4f}<br>"
                    "Mean T24: %{customdata[2]:.4f}<br>"
                    "Mean T72: %{customdata[3]:.4f}<br>"
                    "p-value: %{customdata[4]:.2e}<br>"
                    "q-value: %{customdata[5]:.2e}<extra></extra>"
                ),
            ))
        return traces

    _experiments = _df_plot.index.get_level_values("Experiment").unique().sort_values()

    fig = make_subplots(
        rows=1, cols=len(_experiments),
        subplot_titles=list(_experiments),
        shared_yaxes=True,
        horizontal_spacing=0.06,
    )

    for _col, _exp in enumerate(_experiments, start=1):
        _df_exp = _df_plot.xs(_exp, level="Experiment")
        for trace in _volcano_traces(_df_exp, showlegend=(_col == 1)):
            fig.add_trace(trace, row=1, col=_col)

        # q=0.05 threshold line
        _ty = _q_threshold_y(_df_exp)
        if _ty is not None:
            fig.add_hline(
                y=_ty, row=1, col=_col,
                line=dict(color="#898781", width=1.5, dash="dash"),
                annotation=dict(text="q = 0.05", font=dict(color="#898781", size=10),
                                xanchor="right", x=1),
            )

        # Label the top-10 proteins in this experiment
        _top_exp = df_top10.xs(_exp, level="Experiment") if _exp in df_top10.index.get_level_values("Experiment") else None
        if _top_exp is not None:
            _top_plot = _df_exp.loc[_df_exp.index.isin(_top_exp.index)]
            fig.add_trace(go.Scatter(
                x=_top_plot["log2FC"],
                y=_top_plot["neg_log_p"],
                mode="markers+text",
                text=_top_plot.index.tolist(),
                textposition="top center",
                textfont=dict(size=9, color="#0b0b0b"),
                marker=dict(color="rgba(0,0,0,0)", size=5,
                            line=dict(color="#0b0b0b", width=1.5)),
                showlegend=False,
                hoverinfo="skip",
            ), row=1, col=_col)

    _axis_style = dict(
        showgrid=True, gridcolor="#e1e0d9", gridwidth=1,
        zeroline=True, zerolinecolor="#c3c2b7", zerolinewidth=1,
        tickfont=dict(color="#52514e", size=11),
        title_font=dict(color="#52514e", size=12),
    )
    fig.update_layout(
        paper_bgcolor="#fcfcfb",
        plot_bgcolor="#fcfcfb",
        font=dict(family="system-ui, -apple-system, Segoe UI, sans-serif",
                  color="#0b0b0b"),
        legend=dict(
            title=dict(text="Trend", font=dict(color="#52514e", size=12)),
            font=dict(color="#0b0b0b", size=12),
            bgcolor="rgba(0,0,0,0)", borderwidth=0,
        ),
        height=500,
        margin=dict(t=60, b=60, l=70, r=20),
    )
    fig.update_xaxes(**_axis_style, title_text="log₂FC (T72 / T0)")
    fig.update_yaxes(**_axis_style, title_text="−log₁₀(p-value)", col=1)

    fig
    return go, make_subplots


@app.cell
def _(df_results, np):
    _df_scored = df_results[
        (df_results["Trend"] != "0") &
        (df_results["Mean T0"] > 0) &
        (df_results["Mean T72"] > 0)
    ].copy()
    _df_scored["log2FC"] = np.log2(_df_scored["Mean T72"] / _df_scored["Mean T0"])
    _df_scored["score"]  = -np.log10(_df_scored["p-value"]) * _df_scored["log2FC"].abs()

    df_top10 = (
        _df_scored
        .groupby(level="Experiment", group_keys=False)
        .apply(lambda _grp: _grp.nlargest(10, "score"))
    )
    df_top10[["log2FC", "p-value", "q-value", "Trend", "score"]]
    return (df_top10,)


@app.cell
def _(df_merged, df_top10, go, make_subplots, np, pd):
    _CAT8 = [
        "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
        "#e87ba4", "#008300", "#4a3aa7", "#e34948",
    ]
    _T_COLS = [
        ["T0-1",  "T0-2",  "T0-3"],
        ["T24-1", "T24-2", "T24-3"],
        ["T72-1", "T72-2", "T72-3"],
    ]
    _T_LABELS = ["T0", "T24", "T72"]
    _EXPS    = sorted(df_top10.index.get_level_values("Experiment").unique())
    _TRENDS  = ["+", "-"]
    _TR_LBL  = {"+": "Increasing", "-": "Decreasing"}

    fig_lines = make_subplots(
        rows=len(_TRENDS), cols=len(_EXPS),
        subplot_titles=[f"{_e} — {_TR_LBL[_t]}" for _t in _TRENDS for _e in _EXPS],
        shared_yaxes="rows",
        shared_xaxes=True,
        horizontal_spacing=0.10,
        vertical_spacing=0.14,
    )

    for _ri, _trend in enumerate(_TRENDS, start=1):
        for _ci, _exp in enumerate(_EXPS, start=1):
            _subset = (
                df_top10
                .xs(_exp, level="Experiment")
                .pipe(lambda _d: _d[_d["Trend"] == _trend])
            )
            for _i, (_gene, _row) in enumerate(_subset.iterrows()):
                _color = _CAT8[_i % 8]
                _dash  = "solid" if _i < 8 else "dash"

                _raw = df_merged.loc[_gene]
                if isinstance(_raw, pd.Series):
                    _raw = _raw.to_frame().T
                _raw = _raw[_raw["Experiment"] == _exp]

                _means, _sems = [], []
                for _cols in _T_COLS:
                    _vals = _raw[_cols].astype(float).values.flatten()
                    _vals = _vals[~np.isnan(_vals)]
                    _means.append(_vals.mean())
                    _sems.append(_vals.std(ddof=1) / np.sqrt(len(_vals)) if len(_vals) > 1 else 0)

                fig_lines.add_trace(go.Scatter(
                    x=_T_LABELS,
                    y=_means,
                    error_y=dict(type="data", array=_sems, visible=True,
                                 color=_color, thickness=1.5, width=4),
                    mode="lines+markers+text",
                    name=_gene,
                    text=["", "", _gene],
                    textposition="middle right",
                    textfont=dict(color=_color, size=9,
                                  family="system-ui, -apple-system, Segoe UI, sans-serif"),
                    line=dict(color=_color, width=2, dash=_dash),
                    marker=dict(color=_color, size=7,
                                line=dict(color="#fcfcfb", width=1.5)),
                    showlegend=False,
                    hovertemplate=(
                        f"<b>{_gene} ({_exp})</b><br>"
                        "%{x}: %{y:.4f} ± %{error_y.array:.4f}<extra></extra>"
                    ),
                ), row=_ri, col=_ci)

    _ax = dict(showgrid=True, gridcolor="#e1e0d9", gridwidth=1,
               zeroline=False,
               tickfont=dict(color="#52514e", size=10),
               title_font=dict(color="#52514e", size=11))

    fig_lines.update_layout(
        paper_bgcolor="#fcfcfb", plot_bgcolor="#fcfcfb",
        font=dict(family="system-ui, -apple-system, Segoe UI, sans-serif",
                  color="#0b0b0b"),
        showlegend=False,
        height=750,
        margin=dict(t=60, b=60, l=70, r=160),
    )
    fig_lines.update_xaxes(**_ax)
    fig_lines.update_yaxes(**_ax)
    for _ri in range(1, len(_TRENDS) + 1):
        fig_lines.update_yaxes(title_text="Mean abundance", row=_ri, col=1)
    for _ci in range(1, len(_EXPS) + 1):
        fig_lines.update_xaxes(title_text="Time", row=len(_TRENDS), col=_ci)

    fig_lines
    return


@app.cell
def _(DATA_DIR, df_results):
    _out = DATA_DIR / "results.csv"
    df_results.to_csv(_out)
    print(f"Saved {len(df_results)} rows to {_out}")
    return


if __name__ == "__main__":
    app.run()
