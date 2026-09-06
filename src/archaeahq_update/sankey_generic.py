#!/usr/bin/env python3
"""
Generic Sankey diagram for genomic (or any) database construction pipelines.

Reads a tab-separated input file (.tsv) and produces:
  • An interactive HTML file (always)
  • A 300 DPI PNG image      (always, requires kaleido + Pillow)

Usage:
    python3 sankey_generic.py <input.tsv> [output_stem]

    <input.tsv>    path to the tab-separated input file
    [output_stem]  optional base name for outputs (without extension).
                   Defaults to the input file name without its extension.

Outputs:
    <output_stem>.html   interactive Plotly diagram
    <output_stem>.png    300 DPI raster image for publication / print

Dependencies:
    pip install plotly "kaleido==0.2.1" Pillow
"""

import sys
import re
import io
import warnings
from pathlib import Path

import plotly.graph_objects as go

# Optional: Pillow is used to embed 300 DPI metadata in the PNG
try:
    from PIL import Image
    HAS_PILLOW = True
except ImportError:
    HAS_PILLOW = False

# ──────────────────────────────────────────────────────────────────────────────
# PNG EXPORT SETTINGS
# ──────────────────────────────────────────────────────────────────────────────

PNG_DPI    = 300          # target DPI for the output PNG
SCREEN_DPI = 96           # Plotly's effective screen resolution
PNG_SCALE  = PNG_DPI / SCREEN_DPI   # ≈ 3.125  →  multiplies pixel dimensions

# ──────────────────────────────────────────────────────────────────────────────
# COLOR PALETTE  (cycles automatically for datasets with > 10 groups)
# ──────────────────────────────────────────────────────────────────────────────

PALETTE_HEX = [
    '#2E7D32',   # dark green
    '#B71C1C',   # dark red
    '#E65100',   # dark orange
    '#6A1B9A',   # dark purple
    '#0277BD',   # dark blue
    '#F57F17',   # dark amber
    '#00695C',   # dark teal
    '#4527A0',   # deep purple
    '#AD1457',   # dark pink
    '#558B2F',   # olive green
]
PALETTE_RGBA = [
    'rgba(46,125,50,{a})',
    'rgba(183,28,28,{a})',
    'rgba(230,81,0,{a})',
    'rgba(106,27,154,{a})',
    'rgba(2,119,189,{a})',
    'rgba(245,127,23,{a})',
    'rgba(0,105,92,{a})',
    'rgba(69,39,160,{a})',
    'rgba(173,20,87,{a})',
    'rgba(85,139,47,{a})',
]

REJECT_HEX  = '#546E7A'                 # blue-grey  – rejected / removed nodes
REJECT_RGBA = 'rgba(84,110,122,{a})'

ALL_HEX     = '#1565C0'                 # deep blue  – start / end nodes

LINK_ALPHA  = 0.4                       # transparency of flow ribbons


# ──────────────────────────────────────────────────────────────────────────────
# INPUT PARSING
# ──────────────────────────────────────────────────────────────────────────────

def parse_input(path: Path):
    """
    Parse a tab-separated input file.

    File format
    -----------
    Optional metadata lines (start with #):
        # TITLE:       <diagram title>
        # START_LABEL: <label for the first merged node>
        # END_LABEL:   <label for the final merged node>

    Header row (first non-comment line):
        group <TAB> <stage_0_name> <TAB> <stage_1_name> <TAB> ...

    Data rows (one per group):
        <group_name> <TAB> <count_0> <TAB> <count_1> <TAB> ...

    Rules:
        • Counts must be non-negative integers.
        • Counts must be non-increasing across stages (items can only be
          removed, never added, between consecutive stages).
        • At least 2 stage columns are required (one initial + one filter).

    Returns
    -------
    meta   : dict  – title, start_label, end_label
    groups : list  – group names in file order
    stages : list  – stage column names (all columns after 'group')
    data   : dict  – {group_name: [count_stage0, count_stage1, ...]}
    """
    meta = {
        'title':       'Database Construction Pipeline',
        'start_label': 'All Items',
        'end_label':   'Final Dataset',
    }
    groups, stages, data = [], [], {}

    with open(path, encoding='utf-8') as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip('\n')
            if not line.strip():
                continue

            # ── Metadata lines ────────────────────────────────────────────────
            if line.startswith('#'):
                if ':' in line:
                    key, _, value = line[1:].partition(':')
                    key   = key.strip().lower().replace(' ', '_')
                    value = value.strip()
                    if key in meta:
                        meta[key] = value
                continue

            cols = line.split('\t')

            # ── Header row ────────────────────────────────────────────────────
            if not stages:
                if len(cols) < 2:
                    raise ValueError(
                        f"Line {lineno}: header must have at least 2 "
                        "tab-separated columns ('group' + one stage)."
                    )
                stages = [c.strip() for c in cols[1:]]
                continue

            # ── Data rows ─────────────────────────────────────────────────────
            if len(cols) != len(stages) + 1:
                raise ValueError(
                    f"Line {lineno}: expected {len(stages) + 1} columns "
                    f"(group + {len(stages)} stage counts), got {len(cols)}."
                )
            group  = cols[0].strip()
            counts = []
            for s_idx, raw_val in enumerate(cols[1:], 1):
                try:
                    counts.append(int(raw_val.replace(',', '').strip()))
                except ValueError:
                    raise ValueError(
                        f"Line {lineno}, column {s_idx + 1}: "
                        f"'{raw_val.strip()}' is not a valid integer."
                    )

            for i in range(1, len(counts)):
                if counts[i] > counts[i - 1]:
                    raise ValueError(
                        f"Group '{group}': count at stage {i} ({counts[i]:,}) "
                        f"exceeds count at stage {i - 1} ({counts[i - 1]:,}). "
                        "Counts must be non-increasing (items can only be removed)."
                    )

            groups.append(group)
            data[group] = counts

    if not groups:
        raise ValueError("No data rows found in the input file.")
    if len(stages) < 2:
        raise ValueError(
            "At least 2 stage columns are required "
            "(one 'Initial' column + at least one filter column)."
        )

    return meta, groups, stages, data


# ──────────────────────────────────────────────────────────────────────────────
# LAYOUT HELPERS
# ──────────────────────────────────────────────────────────────────────────────

def y_centers(heights, gap_frac=0.02):
    """
    Compute y-center positions (0 = top, 1 = bottom) for a column of nodes
    stacked proportionally with uniform gaps.

    heights  : node heights in top-to-bottom display order
    gap_frac : gap size as a fraction of the total heights sum
    """
    total_h = sum(heights)
    gap_h   = gap_frac * total_h
    denom   = total_h + (len(heights) - 1) * gap_h
    ys, cum = [], 0.0
    for h in heights:
        ys.append((cum + h / 2.0) / denom)
        cum += h + gap_h
    return ys


def weighted_center(heights, ys):
    """Weighted mean of y positions; returns 0.5 if total weight is zero."""
    total = sum(heights)
    if total == 0:
        return 0.5
    return sum(h * y for h, y in zip(heights, ys)) / total


def short_name(full_name):
    """
    Derive a short display label by keeping only text before the first
    dash, en-dash, parenthesis, or slash.

    'Methanobacteriati-Euryarchaeota' → 'Methanobacteriati'
    'Nanobdellati-DPANN'              → 'Nanobdellati'
    'Bacteria (all)'                  → 'Bacteria'
    """
    parts = re.split(r'[-–/(]', full_name, maxsplit=1)
    return parts[0].strip() if parts else full_name


# ──────────────────────────────────────────────────────────────────────────────
# SANKEY BUILDER
# ──────────────────────────────────────────────────────────────────────────────

def build_sankey(meta, groups, stages, data):
    """
    Build and return a Plotly Figure containing the Sankey diagram.

    Node index scheme
    -----------------
    0              → start node  (All Items)
    1 … n_g        → initial group nodes  (stage index 0)

    For each filter step f  (f = 1 … n_f,  where n_f = len(stages) - 1):
        base(f) = 1 + n_g + (f - 1) × (n_g + 1)
        base(f) … base(f) + n_g - 1  → passed-group nodes at step f
        base(f) + n_g                 → rejected node at step f

    last = 1 + n_g + n_f × (n_g + 1)  → end node  (Final Dataset)
    """
    n_g = len(groups)       # number of groups
    n_f = len(stages) - 1   # number of filter steps

    # One color per group; cycle palette if more groups than colors
    g_hex  = [PALETTE_HEX[i % len(PALETTE_HEX)]  for i in range(n_g)]
    g_rgba = [PALETTE_RGBA[i % len(PALETTE_RGBA)] for i in range(n_g)]

    # Display order within every column: largest group at top
    disp = sorted(range(n_g), key=lambda g: data[groups[g]][0], reverse=True)

    # ── Node index helpers ────────────────────────────────────────────────────
    def base(f):
        return 1 + n_g + (f - 1) * (n_g + 1)

    N       = 1 + n_g + n_f * (n_g + 1) + 1
    end_idx = N - 1

    labels  = [''] * N
    x_pos   = [0.0] * N
    y_pos   = [0.5] * N
    n_color = ['#888'] * N

    # ── Evenly spaced x positions ─────────────────────────────────────────────
    # Columns: start | initial groups | filter_1 | … | filter_n_f | end
    n_slots = n_f + 3
    step    = 0.98 / (n_slots - 1)
    x_slots = [round(0.01 + i * step, 4) for i in range(n_slots)]
    x_slots[0]  = 0.01
    x_slots[-1] = 0.99

    # ── Stage 0: initial groups ───────────────────────────────────────────────
    heights0 = [data[groups[g]][0] for g in disp]
    ys0      = y_centers(heights0)
    g_ys     = {0: {g: ys0[rank] for rank, g in enumerate(disp)}}

    for rank, g in enumerate(disp):
        i          = 1 + g
        labels[i]  = f'{groups[g]}<br><b>{data[groups[g]][0]:,}</b>'
        x_pos[i]   = x_slots[1]
        y_pos[i]   = ys0[rank]
        n_color[i] = g_hex[g]

    # Start node – vertically centered on the initial group distribution
    total_start = sum(data[groups[g]][0] for g in range(n_g))
    labels[0]   = f'{meta["start_label"]}<br><b>{total_start:,}</b>'
    x_pos[0]    = x_slots[0]
    y_pos[0]    = weighted_center(heights0, ys0)
    n_color[0]  = ALL_HEX

    # ── Filter stages ─────────────────────────────────────────────────────────
    for f in range(1, n_f + 1):
        b         = base(f)
        curr      = [data[groups[g]][f]     for g in range(n_g)]
        prev      = [data[groups[g]][f - 1] for g in range(n_g)]
        rejected  = [p - c for p, c in zip(prev, curr)]
        total_rej = sum(rejected)

        heights_f = [data[groups[g]][f] for g in disp] + [total_rej]
        ys_f      = y_centers(heights_f)
        g_ys[f]   = {g: ys_f[rank] for rank, g in enumerate(disp)}

        for rank, g in enumerate(disp):
            i          = b + g
            labels[i]  = f'{short_name(groups[g])}<br><b>{data[groups[g]][f]:,}</b>'
            x_pos[i]   = x_slots[f + 1]
            y_pos[i]   = ys_f[rank]
            n_color[i] = g_hex[g]

        rej_i          = b + n_g
        labels[rej_i]  = f'Rejected<br><b>{total_rej:,}</b>'
        x_pos[rej_i]   = x_slots[f + 1]
        y_pos[rej_i]   = ys_f[n_g]
        n_color[rej_i] = REJECT_HEX

    # End node – vertically centered on the last filter's group distribution
    heights_end = [data[groups[g]][n_f] for g in disp]
    ys_end      = [g_ys[n_f][g] for g in disp]
    total_end   = sum(heights_end)

    labels[end_idx]  = f'{meta["end_label"]}<br><b>{total_end:,}</b>'
    x_pos[end_idx]   = x_slots[-1]
    y_pos[end_idx]   = weighted_center(heights_end, ys_end)
    n_color[end_idx] = ALL_HEX

    # ── Links ─────────────────────────────────────────────────────────────────
    src, tgt, val, l_col, l_lbl = [], [], [], [], []

    def link(s, t, v, rgba_tmpl, lbl=''):
        src.append(s)
        tgt.append(t)
        val.append(v)
        l_col.append(rgba_tmpl.format(a=LINK_ALPHA))
        l_lbl.append(lbl)

    # Start → initial groups
    for g in range(n_g):
        link(0, 1 + g, data[groups[g]][0], g_rgba[g], groups[g])

    # Between consecutive stages
    for f in range(n_f):
        for g in range(n_g):
            s        = (1 + g) if f == 0 else (base(f) + g)
            passed   = data[groups[g]][f + 1]
            rejected = data[groups[g]][f] - passed
            link(s, base(f + 1) + g,   passed,   g_rgba[g],   f'{groups[g]} – passed')
            link(s, base(f + 1) + n_g, rejected, REJECT_RGBA, f'{groups[g]} – rejected')

    # Last filter stage → end node
    for g in range(n_g):
        link(base(n_f) + g, end_idx, data[groups[g]][n_f], g_rgba[g], groups[g])

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = go.Figure(go.Sankey(
        arrangement='fixed',
        node=dict(
            pad=8,
            thickness=22,
            line=dict(color='white', width=0.5),
            label=labels,
            color=n_color,
            x=x_pos,
            y=y_pos,
            hovertemplate='%{label}<extra></extra>',
        ),
        link=dict(
            source=src,
            target=tgt,
            value=val,
            color=l_col,
            label=l_lbl,
            hovertemplate='%{label}: <b>%{value:,}</b><extra></extra>',
        ),
    ))

    # ── Stage header annotations (above the diagram) ──────────────────────────
    header_labels = (
        [(x_slots[0],       meta['start_label'])]
        + [(x_slots[i + 1], stages[i]) for i in range(len(stages))]
        + [(x_slots[-1],    meta['end_label'])]
    )
    for x_ann, text in header_labels:
        fig.add_annotation(
            x=x_ann, y=1.07,
            xref='paper', yref='paper',
            text=f'<b>{text}</b>',
            showarrow=False,
            font=dict(size=10, color='#333333'),
            align='center',
            bgcolor='rgba(255,255,255,0.85)',
        )

    # ── Retention rate annotations (below the diagram) ────────────────────────
    for f in range(1, n_f + 1):
        total_f = sum(data[groups[g]][f] for g in range(n_g))
        pct     = 100 * total_f / total_start if total_start else 0
        fig.add_annotation(
            x=x_slots[f + 1], y=-0.07,
            xref='paper', yref='paper',
            text=f'Passed: {total_f:,} ({pct:.1f}%)',
            showarrow=False,
            font=dict(size=9, color='#555555'),
            align='center',
        )

    # ── Layout ────────────────────────────────────────────────────────────────
    fig.update_layout(
        title=dict(
            text=f'<b>{meta["title"]}</b>',
            font=dict(size=18, color='#1a1a1a'),
            x=0.5,
            xanchor='center',
            y=0.97,
        ),
        font=dict(size=11, family='Arial, sans-serif'),
        paper_bgcolor='white',
        plot_bgcolor='white',
        width=max(1200, 280 * (n_f + 3)),   # scales with number of filter stages
        height=max(500,  130 * n_g),         # scales with number of groups
        margin=dict(t=110, b=60, l=10, r=10),
    )

    return fig, total_start, total_end


# ──────────────────────────────────────────────────────────────────────────────
# IMAGE EXPORT
# ──────────────────────────────────────────────────────────────────────────────

def save_png(fig, path: Path, dpi: int = PNG_DPI):
    """
    Export the figure as a PNG image at the requested DPI.

    The pixel dimensions are scaled proportionally from the Plotly figure size
    (which is defined at SCREEN_DPI = 96 dpi).  Pillow is then used to embed
    the correct DPI metadata so the image opens at the right physical size in
    Word, Illustrator, ImageJ, etc.

    Parameters
    ----------
    fig  : plotly Figure
    path : output file path  (.png)
    dpi  : target resolution in dots per inch  (default 300)
    """
    scale = dpi / SCREEN_DPI               # ≈ 3.125 for 300 DPI
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=DeprecationWarning)
        img_bytes = fig.to_image(format='png', scale=scale)

    if HAS_PILLOW:
        img = Image.open(io.BytesIO(img_bytes))
        img.save(str(path), dpi=(dpi, dpi))
    else:
        # Pillow not available: write the raw bytes (no DPI metadata embedded)
        path.write_bytes(img_bytes)
        print(
            "  Note: Pillow not installed – DPI metadata was not embedded.\n"
            "        Install with:  pip install Pillow"
        )


# ──────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ──────────────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    input_path = Path(sys.argv[1])
    stem       = Path(sys.argv[2]) if len(sys.argv) > 2 else input_path.with_suffix('')
    html_path  = stem.with_suffix('.html')
    png_path   = stem.with_suffix('.png')

    if not input_path.exists():
        print(f"Error: input file not found: {input_path}")
        sys.exit(1)

    print(f"Reading : {input_path}")
    meta, groups, stages, data = parse_input(input_path)

    n_f = len(stages) - 1
    print(f"Groups  : {len(groups)}")
    print(f"Filters : {n_f}  ({', '.join(stages[1:])})")

    fig, total_start, total_end = build_sankey(meta, groups, stages, data)

    # ── HTML ──────────────────────────────────────────────────────────────────
    fig.write_html(str(html_path), include_plotlyjs='cdn')
    print(f"Saved   : {html_path}")

    # ── PNG ───────────────────────────────────────────────────────────────────
    try:
        save_png(fig, png_path, dpi=PNG_DPI)
        layout  = fig.layout
        px_w    = round(layout.width  * PNG_SCALE)
        px_h    = round(layout.height * PNG_SCALE)
        print(f"Saved   : {png_path}  ({px_w} × {px_h} px @ {PNG_DPI} DPI)")
    except Exception as exc:
        print(f"Warning : PNG export failed – {exc}")
        print("          Install kaleido with:  pip install 'kaleido==0.2.1'")

    # ── Pipeline summary ──────────────────────────────────────────────────────
    col_w = max(len(s) for s in stages) + 2
    print(f"\n{'Stage':<{col_w}}  {'Remaining':>10}  {'Removed':>10}  {'% of start':>10}")
    print('-' * (col_w + 36))

    prev_total = total_start
    for i, stage in enumerate(stages):
        curr_total   = sum(data[groups[g]][i] for g in range(len(groups)))
        removed      = prev_total - curr_total if i > 0 else 0
        removed_str  = f'-{removed:,}' if removed else '—'
        pct          = 100 * curr_total / total_start if total_start else 0
        print(f"{stage:<{col_w}}  {curr_total:>10,}  {removed_str:>10}  {pct:>9.1f}%")
        prev_total   = curr_total

    print('-' * (col_w + 36))
    print(f"Overall retention: {100 * total_end / total_start:.1f}%")


if __name__ == '__main__':
    main()
