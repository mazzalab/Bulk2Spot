#!/usr/bin/env python3
"""Bulk2Spot HTML report: one self-contained, interactive HTML page built from
the files the R stages wrote (QC tables, phenoData, PCA/UMAP exports, DEG/GSEA
tables, deconvolution results). Nothing is re-analysed here, with one
exception: a quick "overview" PCA of the Q3-normalized data (before batch
correction), shown next to the final PCA/UMAP exported by stage 3.
"""

from __future__ import annotations

import argparse
import base64
import html
import re
from datetime import datetime
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs

PLOT_CONFIG = {"displaylogo": False, "responsive": True, "modeBarButtonsToRemove": ["lasso2d", "select2d"]}
PALETTE = {"UP": "#4575B4", "DOWN": "#D73027", "NS": "#b8b2a4"}
CATEGORICAL = ["#1f6f6b", "#ae3b2e", "#4575B4", "#e69f00", "#7570b3", "#66a61e", "#e7298a", "#a6761d"]
FONT = dict(family='-apple-system, "Segoe UI", Helvetica, Arial, sans-serif', size=12)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-name", required=True)
    p.add_argument("--preprocessing-dir", required=True, help="results/preprocessing")
    p.add_argument("--spe-dir", required=True, help="results/spe")
    p.add_argument("--deg-dir", required=True, help="results/deg_gsea")
    p.add_argument("--batch-dir", required=True, help="results/batch_correction")
    p.add_argument("--group-column", required=True, help="config comparisons.column, e.g. 'class'")
    p.add_argument("--p-cutoff", type=float, default=0.1, help="config thresholds.p_cutoff (volcano)")
    p.add_argument("--fc-soft", type=float, default=0.57, help="config thresholds.fc_soft (volcano)")
    p.add_argument("--batch-performed", action="store_true", help="Whether batch_correction.perform was true in config")
    p.add_argument("--output", required=True)
    p.add_argument("--top-genes", type=int, default=50)
    p.add_argument("--top-pathways-per-direction", type=int, default=10)
    p.add_argument("--deconv-dir", default=None, help="results/deconvolution (omit to skip that section)")
    p.add_argument("--deconv-id-cols", default="ROI",
                    help="Comma-separated metadata columns in proportions.txt that are NOT cell types "
                         "(config deconvolution.export_columns' keys, plus ROI) -- everything else in "
                         "that table is treated as a cell-type column.")
    p.add_argument("--logo", default=None, help="Bulk2Spot logo mark (SVG), shown in the header and as favicon")
    return p.parse_args(argv)


# Colours of docs/images/logo-mark.svg, swapped for CSS classes when the mark is
# inlined so it follows the report's light/dark theme.
LOGO_MAIN_FILL = 'fill="#0E7C7B"'
LOGO_ACCENT_FILL = 'fill="#F2994A"'


def load_logo(path: Optional[str]) -> tuple:
    """Return (inline header SVG, favicon <link>), or empty strings if there is no logo."""
    if not path or not Path(path).is_file():
        return "", ""
    svg = Path(path).read_text(encoding="utf-8")
    favicon = ('<link rel="icon" type="image/svg+xml" href="data:image/svg+xml;base64,'
               + base64.b64encode(svg.encode("utf-8")).decode("ascii") + '">')
    inline = re.sub(r"<title>.*?</title>\s*", "", svg, flags=re.S).strip()
    inline = re.sub(r'\s(?:width|height|role|aria-label)="[^"]*"', "", inline)
    inline = inline.replace("<svg ", '<svg class="b2s-logo" aria-hidden="true" focusable="false" ', 1)
    inline = inline.replace(LOGO_MAIN_FILL, 'class="b2s-logo-main"').replace(LOGO_ACCENT_FILL, 'class="b2s-logo-accent"')
    return inline, favicon


def read_table(path, **kwargs):
    p = Path(path)
    if not p.exists():
        return None
    return pd.read_csv(p, sep="\t", **kwargs)


def read_xlsx(path, **kwargs):
    p = Path(path)
    if not p.exists():
        return None
    return pd.read_excel(p, **kwargs)


def df_to_html(df, columns=None, n=None, round_cols=None, na_rep="-"):
    if df is None or df.empty:
        return '<p class="missing">No data available.</p>'
    d = df.copy()
    if columns:
        d = d[[c for c in columns if c in d.columns]]
    if n:
        d = d.head(n)
    if round_cols:
        for c in round_cols:
            if c in d.columns:
                d[c] = pd.to_numeric(d[c], errors="coerce").map(lambda v: f"{v:.3g}" if pd.notnull(v) else na_rep)
    return d.to_html(index=False, na_rep=na_rep, classes="data-table", border=0, escape=True)


def fmt_number(x):
    try:
        return f"{int(x):,}"
    except (TypeError, ValueError):
        return str(x)


# ============================================================================
# Plotly figure builders. Each returns an HTML snippet, or a "not available"
# note when its input file is missing.
# ============================================================================

def _layout(fig, height=440, **kwargs):
    kwargs.setdefault("margin", dict(l=60, r=30, t=40, b=50))
    legend = dict(bgcolor="rgba(0,0,0,0)")
    legend.update(kwargs.pop("legend", {}) or {})
    fig.update_layout(
        template="plotly_white",
        font=FONT,
        height=height,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=legend,
        **kwargs,
    )
    return fig


def to_div(fig, div_id):
    return fig.to_html(full_html=False, include_plotlyjs=False, config=PLOT_CONFIG, div_id=div_id)


def missing_plot(label):
    return f'<p class="missing">Figure not available: {label}</p>'


def fig_wrap(div_html, caption=None):
    fig = f'<div class="plot">{div_html}</div>'
    if caption:
        fig = f"<figure>{fig}<figcaption>{caption}</figcaption></figure>"
    return fig


def plot_qc_funnel(qc_log: pd.DataFrame, div_id: str):
    if qc_log is None or qc_log.empty:
        return missing_plot("QC funnel")
    fig = go.Figure(go.Bar(
        x=qc_log["Step"], y=qc_log["Samples"], marker_color="#1f6f6b",
        customdata=qc_log["Features"],
        hovertemplate="<b>%{x}</b><br>Segments: %{y}<br>Features: %{customdata}<extra></extra>",
    ))
    fig.update_xaxes(tickangle=-30, title="")
    fig.update_yaxes(title="Segments remaining")
    _layout(fig, height=380)
    return fig_wrap(to_div(fig, div_id), "Segments remaining after each preprocessing step")


def plot_feature_funnel(qc_log: pd.DataFrame, div_id: str):
    if qc_log is None or qc_log.empty:
        return missing_plot("feature funnel")
    fig = go.Figure(go.Bar(
        x=qc_log["Step"], y=qc_log["Features"], marker_color="#ae3b2e",
        hovertemplate="<b>%{x}</b><br>Features: %{y}<extra></extra>",
    ))
    fig.update_xaxes(tickangle=-30, title="")
    fig.update_yaxes(title="Features (probes/genes) remaining")
    _layout(fig, height=380)
    return fig_wrap(to_div(fig, div_id), "Features remaining after each preprocessing step")


def plot_qc_flags(qc_summary: pd.DataFrame, div_id: str):
    if qc_summary is None or qc_summary.empty:
        return missing_plot("segment QC flags")
    d = qc_summary[qc_summary.index != "TOTAL FLAGS"] if "TOTAL FLAGS" in qc_summary.index else qc_summary
    fig = go.Figure()
    fig.add_trace(go.Bar(x=d.index, y=d["Pass"], name="Pass", marker_color="#1f6f6b"))
    fig.add_trace(go.Bar(x=d.index, y=d["Warning"], name="Warning", marker_color="#ae3b2e"))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(tickangle=-30, title="")
    fig.update_yaxes(title="Segments")
    _layout(fig, height=380)
    return fig_wrap(to_div(fig, div_id), "Segment QC flag outcomes by check (a segment can fail more than one check)")


def plot_composition(pheno: pd.DataFrame, column: str, div_id: str):
    if pheno is None or column not in pheno.columns:
        return missing_plot(f"{column} composition")
    counts = pheno[column].value_counts()
    if counts.empty:
        return missing_plot(f"{column} composition")
    fig = go.Figure(go.Bar(
        x=counts.index.astype(str), y=counts.values, marker_color="#1f6f6b",
        hovertemplate="<b>%{x}</b><br>%{y} segment(s)<extra></extra>",
    ))
    fig.update_xaxes(tickangle=-30, title="")
    fig.update_yaxes(title="Segments")
    _layout(fig, height=360)
    return fig_wrap(to_div(fig, div_id), f"Segments by {column}")


def compute_pca(expr: pd.DataFrame, pheno: pd.DataFrame):
    """log2(Q3-normalized + 1), mean-center genes, SVD -- a lightweight QC
    PCA, not a re-derivation of any pipeline statistical result."""
    common = expr.columns.intersection(pheno.index)
    if len(common) < 3:
        return None
    mat = np.log2(expr[common].to_numpy(dtype=float) + 1)
    mat = mat - mat.mean(axis=1, keepdims=True)
    # genes x segments -> PCA over segments: transpose to segments x genes
    u, s, vt = np.linalg.svd(mat.T, full_matrices=False)
    var_explained = (s ** 2) / np.sum(s ** 2) * 100
    n_pc = min(4, u.shape[1])
    pcs = u[:, :n_pc] * s[:n_pc]
    pca_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(n_pc)], index=common)
    pca_df = pca_df.join(pheno.loc[common])
    return pca_df, var_explained[:n_pc]


def plot_pca(pca_result, color_col: str, div_id: str):
    if pca_result is None:
        return missing_plot("PCA")
    pca_df, var_explained = pca_result
    if color_col not in pca_df.columns:
        return missing_plot("PCA")
    fig = go.Figure()
    groups = sorted(pca_df[color_col].dropna().unique(), key=str)
    for i, g in enumerate(groups):
        sub = pca_df[pca_df[color_col] == g]
        fig.add_trace(go.Scatter(
            x=sub["PC1"], y=sub["PC2"], mode="markers", name=str(g),
            marker=dict(size=9, color=CATEGORICAL[i % len(CATEGORICAL)], opacity=0.85, line=dict(width=1, color="white")),
            text=sub.index,
            hovertemplate=f"<b>%{{text}}</b><br>{color_col}={g}<br>PC1=%{{x:.2f}}, PC2=%{{y:.2f}}<extra></extra>",
        ))
    fig.update_xaxes(title=f"PC1 ({var_explained[0]:.1f}%)", zeroline=True, zerolinecolor="#e2ddd0")
    fig.update_yaxes(title=f"PC2 ({var_explained[1]:.1f}%)", zeroline=True, zerolinecolor="#e2ddd0")
    _layout(fig, height=440)
    return fig_wrap(to_div(fig, div_id), f"PCA of log2(Q3-normalized) expression, colored by {color_col}")


def read_pca_final(batch_dir: Path):
    """Read 03_batch_correction.R's own final-PCA export -- the actual
    SpatialExperiment DEG uses (RUV4(+limma)-corrected if batch correction
    ran, TMM-normalized only if it didn't), not a re-derivation. Same
    (pca_df, var_explained) shape as compute_pca() so plot_pca() below
    renders either one identically."""
    pca_df = read_table(batch_dir / "PCA_final.txt", index_col=0)
    if pca_df is None or pca_df.empty:
        return None
    pc_cols = [c for c in pca_df.columns if re.match(r"^PC\d+$", c)]
    if not pc_cols:
        return None
    var_df = read_table(batch_dir / "PCA_final_variance.txt")
    var_lookup = dict(zip(var_df["PC"], var_df["percentVar"])) if var_df is not None and not var_df.empty else {}
    var_explained = [var_lookup.get(f"PC{i+1}", float("nan")) for i in range(len(pc_cols))]
    return pca_df, var_explained


def read_umap_final(batch_dir: Path):
    umap_df = read_table(batch_dir / "UMAP_final.txt", index_col=0)
    if umap_df is None or umap_df.empty or "UMAP1" not in umap_df.columns or "UMAP2" not in umap_df.columns:
        return None
    return umap_df


def plot_umap(umap_df, color_col: str, div_id: str):
    if umap_df is None or color_col not in umap_df.columns:
        return missing_plot("UMAP")
    fig = go.Figure()
    groups = sorted(umap_df[color_col].dropna().unique(), key=str)
    for i, g in enumerate(groups):
        sub = umap_df[umap_df[color_col] == g]
        fig.add_trace(go.Scatter(
            x=sub["UMAP1"], y=sub["UMAP2"], mode="markers", name=str(g),
            marker=dict(size=9, color=CATEGORICAL[i % len(CATEGORICAL)], opacity=0.85, line=dict(width=1, color="white")),
            text=sub.index,
            hovertemplate=f"<b>%{{text}}</b><br>{color_col}={g}<br>UMAP1=%{{x:.2f}}, UMAP2=%{{y:.2f}}<extra></extra>",
        ))
    fig.update_xaxes(title="UMAP1", zeroline=True, zerolinecolor="#e2ddd0")
    fig.update_yaxes(title="UMAP2", zeroline=True, zerolinecolor="#e2ddd0")
    _layout(fig, height=440)
    return fig_wrap(to_div(fig, div_id), f"UMAP of the post-normalization SpatialExperiment, colored by {color_col}")


def plot_volcano(deg: pd.DataFrame, fc_soft: float, p_cutoff: float, div_id: str, title: str):
    if deg is None or deg.empty:
        return missing_plot("volcano")
    d = deg.copy()
    d["neglog10padj"] = -np.log10(d["adj.P.Val"].clip(lower=np.finfo(float).tiny))
    d["status"] = "NS"
    d.loc[(d["adj.P.Val"] < p_cutoff) & (d["logFC"] > fc_soft), "status"] = "UP"
    d.loc[(d["adj.P.Val"] < p_cutoff) & (d["logFC"] < -fc_soft), "status"] = "DOWN"
    label_col = d.columns[0]  # gene symbol, from the R script's rowNames=TRUE export
    fig = go.Figure()
    for status in ("NS", "DOWN", "UP"):
        sub = d[d["status"] == status]
        if sub.empty:
            continue
        fig.add_trace(go.Scattergl(
            x=sub["logFC"], y=sub["neglog10padj"], mode="markers", name=status,
            marker=dict(size=5, color=PALETTE.get(status, "#999"), opacity=0.65),
            text=sub[label_col],
            customdata=np.stack([sub["AveExpr"], sub["adj.P.Val"]], axis=-1),
            hovertemplate=("<b>%{text}</b><br>logFC=%{x:.2f}<br>AveExpr=%{customdata[0]:.2f}<br>"
                            "adj.P.Val=%{customdata[1]:.2e}<extra></extra>"),
        ))
    fig.add_vline(x=fc_soft, line_dash="dash", line_color="#9b9384")
    fig.add_vline(x=-fc_soft, line_dash="dash", line_color="#9b9384")
    fig.add_hline(y=-np.log10(p_cutoff), line_dash="dash", line_color="#9b9384")
    fig.update_xaxes(title="log2 fold change")
    fig.update_yaxes(title="-log10 adjusted p-value")
    _layout(fig, height=480, title=dict(text=title, font=dict(size=13)))
    return fig_wrap(to_div(fig, div_id), f"Volcano plot &mdash; {html.escape(title)} (hover for gene detail; drag to zoom)")


def plot_gsea(gsea: pd.DataFrame, top_n: int, div_id: str, label: str):
    if gsea is None or gsea.empty or "NES" not in gsea.columns:
        return missing_plot(f"{label} GSEA")
    valid = gsea.dropna(subset=["NES"])
    up = valid[valid["NES"] > 0].sort_values("NES", ascending=False).head(top_n)
    down = valid[valid["NES"] < 0].sort_values("NES").head(top_n)
    d = pd.concat([up, down]).sort_values("NES")
    if d.empty:
        return missing_plot(f"{label} GSEA")
    colors = [PALETTE["UP"] if v > 0 else PALETTE["DOWN"] for v in d["NES"]]
    fig = go.Figure(go.Bar(
        x=d["NES"], y=d["Description"], orientation="h", marker_color=colors,
        customdata=np.stack([d["p.adjust"], d["setSize"]], axis=-1),
        hovertemplate="<b>%{y}</b><br>NES=%{x:.2f}<br>padj=%{customdata[0]:.3g}<br>gene set size=%{customdata[1]:.0f}<extra></extra>",
    ))
    fig.add_vline(x=0, line_color="#9b9384")
    fig.update_xaxes(title="Normalized enrichment score (NES)")
    fig.update_yaxes(tickfont=dict(size=9), automargin=True)
    _layout(fig, height=max(260, 24 * len(d)), margin=dict(l=320, r=30, t=20, b=50))
    return fig_wrap(to_div(fig, div_id), f"Top {label} pathways by NES (positive = up in test group, negative = up in reference group)")


def plot_celltype_stacked_bar(prop: pd.DataFrame, celltype_cols: List[str], sort_col: str, div_id: str):
    """Per-ROI stacked bar of estimated cell-type proportions from
    05_deconvolution.R's proportions.txt (sum-to-one within each ROI, since
    SpatialDecon's prop_of_nontumor is already a composition)."""
    if prop is None or prop.empty or not celltype_cols:
        return missing_plot("cell-type composition")
    d = prop.sort_values([sort_col, "ROI"]) if sort_col in prop.columns else prop.sort_values("ROI")
    fig = go.Figure()
    for i, ct in enumerate(celltype_cols):
        fig.add_trace(go.Bar(x=d["ROI"].astype(str), y=d[ct], name=ct, marker_color=CATEGORICAL[i % len(CATEGORICAL)]))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(title="ROI", tickfont=dict(size=8), tickangle=-60)
    fig.update_yaxes(title="Estimated proportion")
    _layout(fig, height=480, legend=dict(orientation="h", y=-0.55, font=dict(size=9)))
    return fig_wrap(to_div(fig, div_id), "Estimated cell-type composition per ROI (SpatialDecon, safeTME signature)")


def plot_celltype_mean_by_group(prop: pd.DataFrame, celltype_cols: List[str], group_col: str, div_id: str):
    if prop is None or prop.empty or group_col not in prop.columns or not celltype_cols:
        return missing_plot(f"mean cell-type composition by {group_col}")
    means = prop.groupby(group_col)[celltype_cols].mean()
    if means.empty:
        return missing_plot(f"mean cell-type composition by {group_col}")
    fig = go.Figure()
    for i, ct in enumerate(celltype_cols):
        fig.add_trace(go.Bar(x=means.index.astype(str), y=means[ct], name=ct, marker_color=CATEGORICAL[i % len(CATEGORICAL)]))
    fig.update_layout(barmode="stack")
    fig.update_xaxes(title=group_col)
    fig.update_yaxes(title="Mean estimated proportion")
    _layout(fig, height=420, legend=dict(orientation="h", y=-0.4, font=dict(size=9)))
    return fig_wrap(to_div(fig, div_id), f"Mean cell-type composition by {group_col}")


def plot_gsva_heatmap(gsva_long: pd.DataFrame, div_id: str):
    """Signature x ROI heatmap, z-scored per signature (row) -- same axis
    convention as 07_deconvolution_gsva.R's own pheatmap(..., scale='row')."""
    if gsva_long is None or gsva_long.empty:
        return missing_plot("GSVA heatmap")
    wide = gsva_long.pivot_table(index="Signature", columns="ROI", values="Score", aggfunc="first")
    if wide.empty:
        return missing_plot("GSVA heatmap")
    mat = wide.to_numpy(dtype=float)
    mu = np.nanmean(mat, axis=1, keepdims=True)
    sd = np.nanstd(mat, axis=1, keepdims=True)
    sd[sd == 0] = 1
    mat_z = (mat - mu) / sd
    fig = go.Figure(go.Heatmap(
        z=mat_z, x=wide.columns.astype(str), y=wide.index.astype(str),
        colorscale="RdBu_r", zmid=0,
        colorbar=dict(title="z-score"),
        hovertemplate="<b>%{y}</b><br>%{x}<br>z-score=%{z:.2f}<extra></extra>",
    ))
    fig.update_xaxes(showticklabels=False, title="ROI")
    fig.update_yaxes(tickfont=dict(size=9), automargin=True)
    _layout(fig, height=max(320, 22 * len(wide)), margin=dict(l=220, r=30, t=20, b=50))
    return fig_wrap(to_div(fig, div_id), "GSVA enrichment heatmap across ROIs (z-scored per signature)")


def read_stats_dir(stats_dir: Path):
    """Concatenate every Stats_<gcol>.txt found in a deconvolution stats
    directory (06_deconvolution_stats.R's or 07_deconvolution_gsva.R's own
    output), tagging each row with which grouping column it came from."""
    if stats_dir is None or not stats_dir.exists():
        return None
    frames = []
    for f in sorted(stats_dir.glob("Stats_*.txt")):
        gcol = f.stem[len("Stats_"):]
        df = read_table(f)
        if df is None or df.empty:
            continue
        df = df.copy()
        df.insert(0, "grouping_column", gcol)
        frames.append(df)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True).sort_values("pvalue")


CSS = """
:root {
  --bg: #fbf9f5; --fg: #211f1b; --muted: #726a5c; --card: #f2efe8; --border: #e2ddd0;
  --accent: #1f6f6b; --accent-ink: #ffffff; --up: #ae3b2e; --down: #2f5fa8;
  --serif: "Iowan Old Style", "Palatino Linotype", Palatino, Georgia, "Times New Roman", serif;
  --sans: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root { --bg: #13181a; --fg: #e9e5da; --muted: #9b9384; --card: #1b2124; --border: #2a3134; --accent: #4fb3ab; --accent-ink: #0b1414; }
}
:root[data-theme="dark"] { --bg: #13181a; --fg: #e9e5da; --muted: #9b9384; --card: #1b2124; --border: #2a3134; --accent: #4fb3ab; --accent-ink: #0b1414; }
:root[data-theme="light"] { --bg: #fbf9f5; --fg: #211f1b; --muted: #726a5c; --card: #f2efe8; --border: #e2ddd0; --accent: #1f6f6b; --accent-ink: #ffffff; }
* { box-sizing: border-box; }
body {
  counter-reset: figure;
  background: var(--bg); color: var(--fg); font-family: var(--sans);
  max-width: 1000px; margin: 0 auto; padding: 3rem 1.5rem 6rem; line-height: 1.55;
  font-variant-numeric: tabular-nums;
}
h1, h2, h3 { font-family: var(--serif); text-wrap: balance; font-weight: 600; }
h1 { font-size: 2.1rem; margin: .15rem 0 .3rem; letter-spacing: -.01em; }
h2 { font-size: 1.4rem; border-bottom: 1px solid var(--border); padding-bottom: .5rem; margin-top: 3.5rem; }
h3 { font-size: 1.15rem; color: var(--accent); margin-top: 2.25rem; font-style: italic; }
h4 { font-size: .92rem; font-family: var(--sans); font-weight: 700; margin: 1.25rem 0 .5rem; }
.eyebrow { font-family: var(--sans); font-size: .75rem; font-weight: 700; letter-spacing: .12em;
           text-transform: uppercase; color: var(--accent); margin: 0; }
.subtitle { color: var(--muted); margin-top: 0; font-family: var(--serif); font-style: italic; font-size: 1.05rem; }
.meta { color: var(--muted); font-size: .82rem; }
.card { background: var(--card); border: 1px solid var(--border); border-radius: 6px; padding: 1.1rem 1.3rem; margin: 1.25rem 0; }
.stat-grid { display: flex; flex-wrap: wrap; gap: .75rem; margin: 1.25rem 0; }
.stat { background: var(--card); border: 1px solid var(--border); border-radius: 6px; padding: .8rem 1.15rem; min-width: 150px; flex: 1; }
.stat .n { font-size: 1.55rem; font-weight: 700; display: block; font-variant-numeric: tabular-nums; }
.stat .label { font-size: .74rem; color: var(--muted); text-transform: uppercase; letter-spacing: .06em; }
.stat.up .n { color: var(--up); } .stat.down .n { color: var(--down); }
table.data-table { border-collapse: collapse; width: 100%; font-size: .82rem; font-variant-numeric: tabular-nums; }
table.data-table th, table.data-table td { border-bottom: 1px solid var(--border); padding: .4rem .65rem; text-align: left; white-space: nowrap; }
table.data-table th { color: var(--muted); font-weight: 700; font-size: .72rem; text-transform: uppercase; letter-spacing: .04em;
                       position: sticky; top: 0; background: var(--bg); }
table.data-table tr:hover { background: var(--card); }
.table-scroll { overflow-x: auto; max-height: 480px; overflow-y: auto; border: 1px solid var(--border); border-radius: 6px; }
figure { margin: 1.5rem 0; text-align: center; counter-increment: figure; }
figcaption { color: var(--muted); font-size: .8rem; margin-top: .5rem; font-family: var(--serif); font-style: italic; }
figcaption::before { content: "Figure " counter(figure) ". "; font-weight: 700; font-style: normal; color: var(--fg); }
.missing { color: var(--muted); font-style: italic; }
.plot { width: 100%; }
.grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 1.25rem; }
@media (max-width: 720px) { .grid-2 { grid-template-columns: 1fr; } }
.toc { display: flex; flex-wrap: wrap; gap: .4rem .9rem; font-family: var(--sans); font-size: .85rem; }
.toc a { color: var(--accent); text-decoration: none; border-bottom: 1px solid transparent; }
.toc a:hover, .toc a:focus-visible { border-bottom-color: var(--accent); }
a:focus-visible, button:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
code, .mono { font-family: var(--mono); font-size: .85em; }
.badge { display: inline-block; background: var(--card); border: 1px solid var(--border); border-radius: 4px;
         font-family: var(--mono); color: var(--muted);
         padding: .15rem .55rem; font-size: .76rem; margin: 0 .3rem .3rem 0; }
footer { color: var(--muted); font-size: .8rem; margin-top: 4rem; border-top: 1px solid var(--border); padding-top: 1.25rem; }
@media (prefers-reduced-motion: no-preference) { a, .stat { transition: border-color .15s ease; } }
dl.glossary dt { font-weight: 700; font-family: var(--sans); margin-top: .9rem; }
dl.glossary dd { margin: .2rem 0 0; color: var(--muted); }
.header-row { display: flex; justify-content: space-between; align-items: flex-start; gap: 1rem; flex-wrap: wrap; }
.header-text { flex: 1; min-width: 0; }
.brand { display: flex; align-items: center; gap: .55rem; }
.b2s-logo { width: 34px; height: 34px; flex: none; }
.b2s-logo-main { fill: var(--accent); }
.b2s-logo-accent { fill: #e8914a; }
.theme-toggle {
  font-family: var(--sans); font-size: .78rem; font-weight: 600; white-space: nowrap;
  background: var(--card); color: var(--fg); border: 1px solid var(--border);
  border-radius: 999px; padding: .45rem 1rem .45rem .8rem; cursor: pointer;
  display: inline-flex; align-items: center; gap: .45rem; margin-top: .2rem;
}
.theme-toggle:hover, .theme-toggle:focus-visible { border-color: var(--accent); color: var(--accent); }
.theme-toggle svg { width: 15px; height: 15px; flex: none; }
"""

THEME_SCRIPT = r"""
function b2sReportColors(dark) {
  return dark
    ? { font: '#c9c3b4', grid: '#2a3134', line: '#3a4144', zero: '#454c4f' }
    : { font: '#4a453c', grid: '#e6e1d6', line: '#d8d2c4', zero: '#c9c2b2' };
}
function b2sReportRethemePlots(dark) {
  var c = b2sReportColors(dark);
  document.querySelectorAll('.plotly-graph-div').forEach(function (gd) {
    if (!gd.layout || typeof Plotly === 'undefined') return;
    var update = { 'font.color': c.font };
    if (gd.layout.legend) update['legend.font.color'] = c.font;
    Object.keys(gd.layout).forEach(function (key) {
      if (/^(xaxis|yaxis)\d*$/.test(key)) {
        update[key + '.gridcolor'] = c.grid;
        update[key + '.zerolinecolor'] = c.zero;
        update[key + '.linecolor'] = c.line;
        update[key + '.tickfont.color'] = c.font;
      }
    });
    try { Plotly.relayout(gd, update); } catch (e) {}
    if (gd.data) {
      gd.data.forEach(function (trace, i) {
        if (trace.colorbar) {
          try {
            Plotly.restyle(gd, {
              'colorbar.tickfont.color': [c.font],
              'colorbar.title.font.color': [c.font],
            }, [i]);
          } catch (e) {}
        }
      });
    }
  });
}
function b2sReportUpdateButton(theme) {
  var btn = document.getElementById('theme-toggle');
  if (!btn) return;
  btn.setAttribute('aria-pressed', theme === 'dark' ? 'true' : 'false');
  var label = btn.querySelector('.theme-toggle-label');
  if (label) label.textContent = theme === 'dark' ? 'Light mode' : 'Dark mode';
}
function b2sReportSetTheme(theme) {
  document.documentElement.setAttribute('data-theme', theme);
  try { localStorage.setItem('bulk2spot-report-theme', theme); } catch (e) {}
  b2sReportUpdateButton(theme);
  b2sReportRethemePlots(theme === 'dark');
}
(function () {
  var stored = null;
  try { stored = localStorage.getItem('bulk2spot-report-theme'); } catch (e) {}
  var theme = stored || (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
  document.documentElement.setAttribute('data-theme', theme);
  document.addEventListener('DOMContentLoaded', function () {
    b2sReportUpdateButton(theme);
    b2sReportRethemePlots(theme === 'dark');
    var btn = document.getElementById('theme-toggle');
    if (btn) btn.addEventListener('click', function () {
      var current = document.documentElement.getAttribute('data-theme');
      b2sReportSetTheme(current === 'dark' ? 'light' : 'dark');
    });
  });
})();
"""

THEME_TOGGLE_BUTTON = """
<button type="button" id="theme-toggle" class="theme-toggle" aria-pressed="false">
<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
<circle cx="12" cy="12" r="4"></circle>
<path d="M12 2v2M12 20v2M4.93 4.93l1.41 1.41M17.66 17.66l1.41 1.41M2 12h2M20 12h2M4.93 19.07l1.41-1.41M17.66 6.34l1.41-1.41"></path>
</svg>
<span class="theme-toggle-label">Dark mode</span>
</button>
"""

GLOSSARY = [
    ("AOI / ROI / Segment", "In GeoMx DSP, a Region of Interest (ROI) is a user-drawn area on the tissue "
     "slide; within it, an Area of Illumination (AOI) is the specific sub-region UV-exposed and collected "
     "for one barcoded readout &mdash; what this report calls a \"segment\". Segmentation strategy (e.g. "
     "whole geometric shapes vs. marker-positive/negative compartments like \"PanCK+\"/\"PanCK-\") is "
     "chosen per experiment, not fixed by the platform."),
    ("WTA panel", "Whole Transcriptome Atlas &mdash; a GeoMx probe panel targeting the whole protein-coding "
     "transcriptome (~18,000+ genes), as opposed to a smaller targeted/CTA panel."),
    ("Negative control probes / NegGeoMean", "Probes with no complementary target in the panel, used to "
     "estimate background signal per segment. Their geometric mean (NegGeoMean) sets the Limit of "
     "Quantification (LOQ) &mdash; the level a real signal must clear to be called \"detected\" above noise."),
    ("LOQ (Limit of Quantification)", "Per-segment, per-module threshold (NegGeoMean &times; a configured "
     "SD multiplier, floored at a minimum) above which a gene's counts are considered detected rather than "
     "background noise. Genes/segments that don't clear it enough of the time get filtered out."),
    ("Segment / probe QC flags", "NanoString's standard per-segment technical QC (read depth, trimming/"
     "stitching/alignment rate, sequencing saturation, negative control counts, no-template-control "
     "counts, nuclei count, area) and per-probe QC (outlier detection via Grubb's test). Thresholds are "
     "dataset/panel-specific &mdash; always check real per-segment distributions before trusting a "
     "template's defaults, not just copy them from another project."),
    ("Q3 normalization", "GeomxTools' recommended normalization for WTA data: each segment's counts are "
     "scaled so its 75th-percentile (Q3) value matches a common target, correcting for differences in "
     "overall signal intensity between segments before comparing expression levels."),
    ("SpatialExperiment", "A Bioconductor data structure (via the <code>standR</code> package here) that "
     "extends the standard expression-matrix-plus-metadata model with spatial coordinates, used as the "
     "common object downstream analyses (batch correction, DEG) are built on."),
    ("RUV4 / RUV4 batch correction", "Remove Unwanted Variation, 4th formulation &mdash; a method that "
     "estimates and removes technical (batch) variation using negative control genes, while explicitly "
     "preserving a specified biological factor of interest. Only meaningful when a genuinely independent "
     "batch variable exists in the data (not confounded with the biology being studied) &mdash; skipped "
     "entirely, honestly, when that's not the case rather than fabricated."),
    ("TMM normalization", "Trimmed Mean of M-values (edgeR/standR) &mdash; a normalization method "
     "accounting for RNA composition differences between samples, applied here independently of the "
     "earlier Q3 normalization as part of the batch-correction/diagnostics stage."),
    ("UMAP", "Uniform Manifold Approximation and Projection &mdash; a nonlinear dimensionality-reduction "
     "technique, complementary to PCA. Where PCA axes are directly interpretable (variance explained per "
     "component), UMAP prioritizes preserving local neighborhood structure over global distances or axis "
     "meaning &mdash; useful for spotting cluster structure, but distances between well-separated clusters "
     "and the absolute scale of either axis are not meaningful the way PCA's are."),
    ("limma-voom", "A differential expression method that models the mean-variance relationship of "
     "log-counts explicitly (voom) before standard linear-model (limma) testing &mdash; used here per "
     "segment type, per comparison, rather than DESeq2/edgeR's own count-based models."),
    ("log2 fold change (logFC)", "How much a gene's expression differs between the two groups being "
     "compared, on a log2 scale. +1 means doubled expression in the test group; -1 means halved."),
    ("p-value / adjusted p-value (adj.P.Val)", "The p-value is the probability of seeing a difference this "
     "large by chance alone if there were truly no effect. Because thousands of genes are tested at once, "
     "raw p-values are corrected for multiple testing (Benjamini&ndash;Hochberg) to control the false "
     "discovery rate &mdash; that corrected value, <code>adj.P.Val</code>, is what should be used to call "
     "a gene significant, not the raw p-value."),
    ("Volcano plot", "Each dot is one gene: horizontal position is log2 fold change (effect size), "
     "vertical position is statistical significance (-log10 adjusted p-value, so higher = more "
     "significant). Genes in the upper-left/upper-right corners changed a lot AND are statistically "
     "confident."),
    ("GSEA (Gene Set Enrichment Analysis) / NES", "Instead of testing one gene at a time, GSEA asks "
     "whether a whole curated group of biologically related genes (a GO term or KEGG pathway) tends to "
     "sit consistently toward one end of the full ranked gene list. NES (normalized enrichment score) is "
     "the effect size: positive means the pathway is shifted toward the test group, negative toward the "
     "reference group."),
    ("GO Biological Process / KEGG", "Gene Ontology Biological Process terms describe what larger "
     "biological process a gene participates in; KEGG pathways describe curated molecular pathways/"
     "networks. Both are hierarchical and overlapping by design &mdash; a single real biological signal "
     "routinely lights up several closely related terms at once, which is expected, not a sign of error."),
    ("SpatialDecon / safeTME", "SpatialDecon estimates the relative abundance of different cell types "
     "within each ROI's mixed signal, by fitting the ROI's observed expression as a weighted combination "
     "of reference cell-type expression profiles. safeTME is SpatialDecon's own built-in reference "
     "signature (18 immune/stromal cell types), purpose-built for GeoMx WTA data. Results are "
     "<i>estimated relative proportions</i>, not a direct cell count &mdash; and are only as accurate as "
     "how well the real tissue's cell types are represented in the reference signature."),
    ("GSVA (Gene Set Variation Analysis)", "A different technique from the GSEA used in the DEG section "
     "above: rather than testing whether a gene set is enriched among genes ranked by a group comparison, "
     "GSVA computes one enrichment score <i>per sample</i> (here, per ROI) for each gene signature, "
     "independent of any group comparison. That per-ROI score is what downstream significance testing "
     "and the heatmap below are built on."),
    ("Custom gene signature", "A user-supplied list of genes believed to jointly mark a cell state or "
     "process (e.g. a published tumor-microenvironment signature), scored per ROI via GSVA above &mdash; "
     "distinct from safeTME, which is a pre-built reference for estimating cell-type <i>proportions</i> "
     "rather than scoring arbitrary gene sets."),
]


def build_glossary_html():
    items = "".join(f"<dt>{term}</dt><dd>{definition}</dd>" for term, definition in GLOSSARY)
    return f"""
<h2 id="glossary">Glossary</h2>
<p class="meta">Brief, non-technical explanations of the terms and plot types used throughout this report.</p>
<dl class="glossary">{items}</dl>
"""


def build_html(args: argparse.Namespace) -> str:
    prep_dir = Path(args.preprocessing_dir)
    spe_dir = Path(args.spe_dir)
    deg_dir = Path(args.deg_dir)
    batch_dir = Path(args.batch_dir)
    group_col = args.group_column

    qc_log = read_table(prep_dir / "QC_log.txt")
    qc_summary = read_table(prep_dir / "qc_summary_table.txt", index_col=0)
    ntc_summary = read_table(prep_dir / "ntc_count_summary.txt")
    data_q3 = read_table(prep_dir / "data_Q3.txt", index_col=0)

    pheno = read_table(spe_dir / "phenoData.txt", index_col=0)

    pca_result_final = read_pca_final(batch_dir)
    umap_df = read_umap_final(batch_dir)

    deg_summary = read_table(deg_dir / "DEG_GSEA_summary.txt")

    deconv_dir = Path(args.deconv_dir) if args.deconv_dir else None
    deconv_prop = read_table(deconv_dir / "proportions.txt") if deconv_dir else None
    deconv_id_cols = [c.strip() for c in args.deconv_id_cols.split(",") if c.strip()]
    deconv_celltype_cols = [c for c in deconv_prop.columns if c not in deconv_id_cols] if deconv_prop is not None else []
    deconv_gsva_long = read_table(deconv_dir / "gsva_long.txt") if deconv_dir else None
    deconv_stats = read_stats_dir(deconv_dir / "stats") if deconv_dir else None
    deconv_gsva_stats = read_stats_dir(deconv_dir / "GSVA_stats") if deconv_dir else None
    deconv_toc = '<a href="#deconv">Deconvolution</a>' if deconv_dir else ""

    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    logo_svg, favicon = load_logo(args.logo)
    sections: List[str] = []

    # ---------------------------------------------------------------- header
    sections.append(f"""
<div class="header-row">
<div class="header-text">
<div class="brand">{logo_svg}<p class="eyebrow">Bulk2Spot &middot; GeoMx DSP spatial transcriptomics report</p></div>
<h1>{html.escape(args.project_name)}</h1>
<p class="subtitle">Segment/probe QC &middot; Q3 normalization &middot; SpatialExperiment &middot; limma-voom DEG &middot; GO/KEGG GSEA</p>
<p class="meta">Generated {generated}</p>
</div>
{THEME_TOGGLE_BUTTON}
</div>
<nav class="toc card" aria-label="Table of contents">
<a href="#overview">Overview</a>
<a href="#qc">Preprocessing QC</a>
<a href="#composition">Sample composition</a>
<a href="#batch">Batch correction</a>
<a href="#deg">DEG &amp; GSEA</a>
{deconv_toc}
<a href="#methods">Methods &amp; parameters</a>
<a href="#glossary">Glossary</a>
</nav>
""")

    # ---------------------------------------------------------------- overview
    n_segments_final = int(qc_log["Samples"].iloc[-1]) if qc_log is not None and not qc_log.empty else None
    n_segments_raw = int(qc_log["Samples"].iloc[0]) if qc_log is not None and not qc_log.empty else None
    n_genes_final = int(qc_log["Features"].iloc[-1]) if qc_log is not None and not qc_log.empty else None
    n_segment_types = pheno["Segment"].nunique() if pheno is not None and "Segment" in pheno.columns else None
    n_comparisons = deg_summary["Comparison"].nunique() if deg_summary is not None and not deg_summary.empty else 0
    stat_html = "".join([
        f'<div class="stat"><span class="n">{fmt_number(n_segments_raw)}</span><span class="label">Segments (raw)</span></div>' if n_segments_raw is not None else "",
        f'<div class="stat"><span class="n">{fmt_number(n_segments_final)}</span><span class="label">Segments (QC-passed)</span></div>' if n_segments_final is not None else "",
        f'<div class="stat"><span class="n">{fmt_number(n_genes_final)}</span><span class="label">Genes (final)</span></div>' if n_genes_final is not None else "",
        f'<div class="stat"><span class="n">{fmt_number(n_segment_types)}</span><span class="label">Segment type(s)</span></div>' if n_segment_types is not None else "",
        f'<div class="stat"><span class="n">{n_comparisons}</span><span class="label">Comparison(s)</span></div>',
        f'<div class="stat"><span class="n">{"Applied" if args.batch_performed else "Skipped"}</span><span class="label">Batch correction</span></div>',
    ])
    sections.append(f"""
<h2 id="overview">Overview</h2>
<div class="stat-grid">{stat_html}</div>
""")

    # ---------------------------------------------------------------- preprocessing QC
    sections.append(f"""
<h2 id="qc">Preprocessing QC</h2>
<div class="grid-2">
  {plot_qc_funnel(qc_log, "plot-qc-funnel")}
  {plot_feature_funnel(qc_log, "plot-feature-funnel")}
</div>
{plot_qc_flags(qc_summary, "plot-qc-flags")}
""")
    if qc_log is not None:
        sections.append(f"""
<h3>Preprocessing step log</h3>
<div class="table-scroll">{df_to_html(qc_log)}</div>
""")
    if ntc_summary is not None:
        sections.append(f"""
<h3>No-template-control (NTC) read counts</h3>
<div class="table-scroll">{df_to_html(ntc_summary)}</div>
""")

    # ---------------------------------------------------------------- sample composition
    comp_cols = [c for c in [group_col, "Segment", "slide.name", "region", "pathology"] if pheno is not None and c in pheno.columns]
    comp_html = "".join(plot_composition(pheno, c, f"plot-comp-{c}") for c in comp_cols[:4])
    sections.append(f"""
<h2 id="composition">Sample composition (QC-passed segments)</h2>
<div class="grid-2">{comp_html}</div>
""")

    # ---------------------------------------------------------------- batch correction / PCA / UMAP
    pca_result = compute_pca(data_q3, pheno) if data_q3 is not None and pheno is not None else None
    pca_color_cols = [c for c in [group_col, "Segment", "slide.name"] if pheno is not None and c in pheno.columns]
    pca_html = "".join(plot_pca(pca_result, c, f"plot-pca-overview-{c}") for c in pca_color_cols[:3])

    # Final PCA/UMAP: colour by every metadata column stage 3 exported
    # (group, Segment, Scan, Patient and the batch column, if present)
    pc_pattern = re.compile(r"^PC\d+$")
    final_pca_color_cols = (
        [c for c in pca_result_final[0].columns if not pc_pattern.match(c)] if pca_result_final else []
    )
    final_pca_html = "".join(
        plot_pca(pca_result_final, c, f"plot-pca-final-{c}") for c in final_pca_color_cols
    )
    umap_pattern = re.compile(r"^UMAP\d+$")
    umap_color_cols = [c for c in umap_df.columns if not umap_pattern.match(c)] if umap_df is not None else []
    umap_html = "".join(plot_umap(umap_df, c, f"plot-umap-{c}") for c in umap_color_cols)

    batch_note = (
        "RUV4 batch correction (with limma residual pass, if configured) was applied to the SpatialExperiment "
        "used for DEG below."
        if args.batch_performed else
        "Batch correction was <b>not</b> applied (either turned off in config, or no genuinely independent "
        "batch variable exists in this dataset that isn't confounded with the biology being studied) &mdash; "
        "DEG below uses the TMM-normalized, uncorrected SpatialExperiment, and so do the \"Final\" "
        "plots underneath."
    )
    sections.append(f"""
<h2 id="batch">Batch correction diagnostics</h2>
<p class="meta">{batch_note}</p>
<h3>Overview PCA (before batch correction)</h3>
<div class="grid-2">{pca_html}</div>
<h3>Final PCA and UMAP (data used for DEG)</h3>
<div class="grid-2">{final_pca_html}</div>
<div class="grid-2">{umap_html}</div>
""")

    # ---------------------------------------------------------------- DEG / GSEA
    deg_blocks = []
    if deg_summary is not None and not deg_summary.empty:
        for _, row in deg_summary.iterrows():
            segment, comparison = row["Segment"], row["Comparison"]
            cmp_dir = deg_dir / "DEG" / segment / comparison
            deg_tbl = read_xlsx(cmp_dir / f"TableDEG_{comparison}_{segment}.xlsx", index_col=0)
            gsea_bp = read_xlsx(cmp_dir / f"GSEA_BP_results_{comparison}_{segment}.xlsx")
            gsea_kegg = read_xlsx(cmp_dir / f"GSEA_KEGG_results_{comparison}_{segment}.xlsx")

            slug = f"{segment}-{comparison}".replace(" ", "_")
            title = f"{segment} &mdash; {comparison}"

            stat_row = f"""
<div class="stat-grid">
  <div class="stat up"><span class="n">{fmt_number(row.get('N_up'))}</span><span class="label">Up in {html.escape(str(row.get('Group_Test', '')))}</span></div>
  <div class="stat down"><span class="n">{fmt_number(row.get('N_down'))}</span><span class="label">Down (up in {html.escape(str(row.get('Group_Ref', '')))})</span></div>
  <div class="stat"><span class="n">{fmt_number(row.get('N_tested'))}</span><span class="label">Genes tested</span></div>
  <div class="stat"><span class="n">{fmt_number(row.get('N_GO_BP_terms'))}</span><span class="label">GO BP terms</span></div>
  <div class="stat"><span class="n">{fmt_number(row.get('N_KEGG_terms'))}</span><span class="label">KEGG terms</span></div>
</div>"""

            deg_table_html = ""
            if deg_tbl is not None and not deg_tbl.empty:
                d = deg_tbl.reset_index().rename(columns={deg_tbl.index.name or "index": "Symbol"})
                d = d.sort_values("adj.P.Val")
                deg_table_html = df_to_html(
                    d, columns=["Symbol", "logFC", "AveExpr", "P.Value", "adj.P.Val", "DetectionRate"],
                    n=args.top_genes, round_cols=["logFC", "AveExpr", "P.Value", "adj.P.Val", "DetectionRate"],
                )

            if gsea_bp is not None and not gsea_bp.empty:
                cols = [c for c in ["ID", "Description", "setSize", "NES", "pvalue", "p.adjust", "qvalue"] if c in gsea_bp.columns]
                gsea_bp_html = f"""
<h4>GO Biological Process</h4>
{plot_gsea(gsea_bp, args.top_pathways_per_direction, f"plot-gsea-bp-{slug}", "GO BP")}
<div class="table-scroll">{df_to_html(gsea_bp, columns=cols, round_cols=["NES", "pvalue", "p.adjust", "qvalue"])}</div>
"""
            else:
                gsea_bp_html = '<h4>GO Biological Process</h4><p class="missing">No enriched GO BP terms found for this comparison.</p>'

            if gsea_kegg is not None and not gsea_kegg.empty:
                cols = [c for c in ["ID", "Description", "setSize", "NES", "pvalue", "p.adjust", "qvalue"] if c in gsea_kegg.columns]
                gsea_kegg_html = f"""
<h4>KEGG</h4>
{plot_gsea(gsea_kegg, args.top_pathways_per_direction, f"plot-gsea-kegg-{slug}", "KEGG")}
<div class="table-scroll">{df_to_html(gsea_kegg, columns=cols, round_cols=["NES", "pvalue", "p.adjust", "qvalue"])}</div>
"""
            else:
                gsea_kegg_html = '<h4>KEGG</h4><p class="missing">No enriched KEGG pathways found for this comparison.</p>'

            deg_blocks.append(f"""
<h3 id="cmp-{slug}">{title}</h3>
{stat_row}
{plot_volcano(deg_tbl.reset_index() if deg_tbl is not None else None, args.fc_soft, args.p_cutoff, f"plot-volcano-{slug}", title)}
<h4>Top genes by adjusted p-value</h4>
<div class="table-scroll">{deg_table_html}</div>
{gsea_bp_html}
{gsea_kegg_html}
""")

    sections.append(f"""
<h2 id="deg">DEG &amp; GSEA</h2>
{''.join(deg_blocks) if deg_blocks else '<p class="missing">No segment x comparison results found.</p>'}
""")

    # ---------------------------------------------------------------- deconvolution
    # Only shown when --deconv-dir is given
    if deconv_dir:
        n_celltypes = len(deconv_celltype_cols)
        n_rois = len(deconv_prop) if deconv_prop is not None else None
        n_signatures = deconv_gsva_long["Signature"].nunique() if deconv_gsva_long is not None and not deconv_gsva_long.empty else None
        n_sig_celltype = int((deconv_stats["pvalue"] < 0.05).sum()) if deconv_stats is not None and not deconv_stats.empty else 0
        n_sig_signature = int((deconv_gsva_stats["pvalue"] < 0.05).sum()) if deconv_gsva_stats is not None and not deconv_gsva_stats.empty else 0

        deconv_stat_html = "".join([
            f'<div class="stat"><span class="n">{n_celltypes}</span><span class="label">safeTME cell type(s)</span></div>',
            f'<div class="stat"><span class="n">{fmt_number(n_rois)}</span><span class="label">ROIs deconvolved</span></div>' if n_rois is not None else "",
            f'<div class="stat"><span class="n">{fmt_number(n_signatures)}</span><span class="label">GSVA signature(s)</span></div>' if n_signatures is not None else "",
            f'<div class="stat"><span class="n">{n_sig_celltype}</span><span class="label">Significant cell-type comparisons (p&lt;0.05)</span></div>',
            f'<div class="stat"><span class="n">{n_sig_signature}</span><span class="label">Significant signature comparisons (p&lt;0.05)</span></div>',
        ])

        deconv_stats_html = (
            f'<div class="table-scroll">{df_to_html(deconv_stats, columns=["grouping_column", "segment", "celltype", "pvalue"], n=50, round_cols=["pvalue"])}</div>'
            if deconv_stats is not None and not deconv_stats.empty else
            '<p class="missing">No cell-type-proportion significance tests found (or none configured).</p>'
        )
        deconv_gsva_stats_html = (
            f'<div class="table-scroll">{df_to_html(deconv_gsva_stats, columns=["grouping_column", "Segment", "Signature", "pvalue"], n=50, round_cols=["pvalue"])}</div>'
            if deconv_gsva_stats is not None and not deconv_gsva_stats.empty else
            '<p class="missing">No GSVA signature significance tests found (or none configured).</p>'
        )

        sections.append(f"""
<h2 id="deconv">Deconvolution</h2>
<p class="meta">Cell-type composition estimated per ROI (SpatialDecon, safeTME signature) and custom-signature enrichment (GSVA), independent analyses computed from the same normalized expression.</p>
<div class="stat-grid">{deconv_stat_html}</div>
<h3>Cell-type composition</h3>
{plot_celltype_stacked_bar(deconv_prop, deconv_celltype_cols, group_col, "plot-deconv-stacked")}
{plot_celltype_mean_by_group(deconv_prop, deconv_celltype_cols, group_col, "plot-deconv-mean-by-group")}
<h4>Significant cell-type comparisons (all grouping columns, sorted by p-value)</h4>
{deconv_stats_html}
<h3>Custom-signature GSVA enrichment</h3>
{plot_gsva_heatmap(deconv_gsva_long, "plot-gsva-heatmap")}
<h4>Significant signature comparisons (all grouping columns, sorted by p-value)</h4>
{deconv_gsva_stats_html}
""")

    # ---------------------------------------------------------------- methods
    badges = [
        f'<span class="badge">comparison column: {html.escape(group_col)}</span>',
        f'<span class="badge">normalization: Q3 (GeomxTools)</span>',
        f'<span class="badge">batch correction: {"RUV4 (+ limma, if configured)" if args.batch_performed else "not applied"}</span>',
        '<span class="badge">DEG: limma-voom</span>',
        '<span class="badge">GSEA: gseGO (BP) + gseKEGG</span>',
    ]
    if deconv_dir:
        badges.append('<span class="badge">Deconvolution: SpatialDecon (safeTME)</span>')
        badges.append('<span class="badge">Signature enrichment: GSVA</span>')
    sections.append(f"""
<h2 id="methods">Methods &amp; parameters</h2>
<div class="card">
{"".join(badges)}
</div>
""")

    sections.append(build_glossary_html())

    sections.append(f"""
<footer>
Report assembled by Bulk2Spot (<code>generate_report.py</code>) from pipeline outputs in
<code>{html.escape(str(prep_dir))}</code>, <code>{html.escape(str(spe_dir))}</code>,
<code>{html.escape(str(deg_dir))}</code>{f', and <code>{html.escape(str(deconv_dir))}</code>' if deconv_dir else ''}.
Plots are interactive: hover for detail, drag to zoom, double-click to reset, click legend entries to toggle series.
</footer>
""")

    plotlyjs_script = f"<script>{get_plotlyjs()}</script>"

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<script>{THEME_SCRIPT}</script>
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(args.project_name)} &mdash; Bulk2Spot report</title>
{favicon}
<style>{CSS}</style>
{plotlyjs_script}
</head><body>
{''.join(sections)}
</body></html>"""


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    document = build_html(args)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(document, encoding="utf-8")
    print(f"Report written to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
