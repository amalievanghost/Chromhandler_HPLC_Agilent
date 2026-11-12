# =============================================================================
# Broken-axis plots for an EnzymeML document (pyenzyme v2)
#
# This file ONLY defines settings and functions - it does not run anything.
# Use it from the notebook:
#
#     import sys
#     sys.path.append("Input")
#     from enzml_broken_axis import to_broken_axes, plot_replicates
#
#     fig_all, _  = plot_replicates(timecourse_enzymeML)            # title = doc name
#     fig_mean, _ = plot_replicates(timecourse_enzymeML, mean_sd=True)
#
# After editing this file, restart the kernel OR reload it in the notebook,
# otherwise the notebook keeps using the old version:
#
#     import importlib, enzml_broken_axis
#     importlib.reload(enzml_broken_axis)
#     from enzml_broken_axis import to_broken_axes, plot_replicates
# =============================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from brokenaxes import brokenaxes
from matplotlib.patches import Rectangle

# -----------------------------------------------------------------------------
# PLOT SETTINGS  (x-limits are in the document's own time unit!)
# -----------------------------------------------------------------------------
LEFT_XLIM      = (0, 125)      # tmmcGAS: (0, 125),  thscGAS: (0, 200)
RIGHT_XLIM     = (1430, 1450)
LEFT_TICK_STEP = 20            # thscGAS: 50, tmmcGAS: 20
RIGHT_TICKS    = [1440]        # tick(s) shown in the narrow right panel

TITLE_FONTSIZE = 12
FRAME_LW = 1.2
TICK_LW  = 1.2
TICK_LEN = 6
Y_START_AT_ZERO = True   # y-axis starts exactly at 0 (axes meet at 0,0)
EDGE_PAD = 0.04          # let markers at 0 poke out past the axes (axes fraction)

# --- frame / break appearance ------------------------------------------------
FRAME = "box"            # "box"    : one closed frame around the whole plot
                         #            (left spine on the first panel, right spine
                         #            on the last, top+bottom across both)
                         # "open"   : L-shaped, only left + bottom axis lines
                         # "panels" : a separate closed box per panel
BREAK_GAP  = 0.1         # width of the gap between the two panels, as a fraction
                         # of the average panel width. The break is horizontal,
                         # so this is GridSpec's *wspace* - hspace does nothing
                         # here, which is why the gap used to be the default 0.2.
BREAK_MARK_SIZE = 0.012  # length of the diagonal break slashes (figure fraction)
BREAK_MARK_TILT = 45     # angle of the slashes, degrees

# --- dashed connector across the break ---------------------------------------
BRIDGE_BREAK = False     # True : draw a dashed line from the last left point to
                         #        the first right point, running into the gap
                         # False: no dashed lines
BRIDGE_STYLE = dict(linestyle="--", linewidth=1.3, marker="None")

# Colours looked up by species NAME (or id); others get matplotlib defaults
SPECIES_COLORS = {
    "cGAMP": "#FD1FBE",
    "ADP":   "#C185C7",
    "GDP":   "#00F0D4",
    "ATP":   "#8E50D0",
    "GTP":   "#92D050",
}

# One marker per replicate (measurement)
REPLICATE_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
def style_frame(bax, lw=FRAME_LW, color="black"):
    """Spines + ticks, following the FRAME setting. In "box" mode the vertical
    spines are drawn only on the outer edges, so the two panels read as one
    frame interrupted by the break rather than as two separate boxes."""
    axs = bax.axs
    last = len(axs) - 1
    for i, ax in enumerate(axs):
        visible = {
            "bottom": True,
            "top":    FRAME in ("box", "panels"),
            "left":   i == 0 if FRAME in ("box", "open") else True,
            "right":  i == last if FRAME == "box" else FRAME == "panels",
        }
        for side, show in visible.items():
            sp = ax.spines[side]
            sp.set_visible(show)
            sp.set_linewidth(lw)
            sp.set_color(color)
        ax.tick_params(which="both", direction="out", top=False, right=False,
                       width=TICK_LW, length=TICK_LEN, color=color)


# backwards compatibility with the old name
close_brokenaxes_box = style_frame


def _draw_break_marks(bax, d=BREAK_MARK_SIZE, tilt=BREAK_MARK_TILT,
                      lw=FRAME_LW, color="black"):
    """Diagonal slashes where a panel edge meets the break. Drawn by hand (in
    figure coordinates) so they only appear on the axis lines that are actually
    visible - brokenaxes' own marks are switched off via diag_color='none'.
    With a closed frame the break interrupts the top line as well, so the
    slashes are drawn at both ends of each panel edge."""
    fig = bax.axs[0].figure
    w, h = fig.get_size_inches()
    dx = d * np.cos(np.deg2rad(tilt))
    dy = d * np.sin(np.deg2rad(tilt)) * w / h

    n = len(bax.axs)
    for i, ax in enumerate(bax.axs):
        x0, y0, bw, bh = ax.get_position().bounds
        xs = ([x0] if i > 0 else []) + ([x0 + bw] if i < n - 1 else [])
        ys = [y0] + ([y0 + bh] if FRAME in ("box", "panels") else [])
        for x in xs:
            for y in ys:
                fig.add_artist(Line2D((x - dx, x + dx), (y - dy, y + dy),
                                      transform=fig.transFigure, color=color,
                                      lw=lw, clip_on=False))


def _style_bax(bax, xlabel, ylabel):
    bax.set_xlabel(xlabel, labelpad=30)
    bax.set_ylabel(ylabel, labelpad=40)
    bax.axs[0].set_xticks(list(range(LEFT_XLIM[0], LEFT_XLIM[1] + 1, LEFT_TICK_STEP)))
    bax.axs[-1].set_xticks(RIGHT_TICKS)

    if Y_START_AT_ZERO:
        for ax in bax.axs:
            ax.set_ylim(bottom=0)
        _unclip_bottom_left(bax)

    style_frame(bax)
    _draw_break_marks(bax)


def _unclip_bottom_left(bax, pad=EDGE_PAD):
    """Points sitting exactly on y=0 or x=0 would be cut in half by the axes
    edge. Enlarge each panel's clip area slightly downwards (and to the left
    for the first panel) so those markers/error bars are drawn in full,
    without letting data leak into the break gap on the right."""
    for i, ax in enumerate(bax.axs):
        left = -pad if i == 0 else 0.0
        clip = Rectangle((left, -pad), 1 - left, 1 + pad,
                         transform=ax.transAxes, facecolor="none", edgecolor="none")
        for artist in list(ax.lines) + list(ax.collections):
            artist.set_clip_path(clip)


def split_by_break(x, y, yerr=None):
    """Split points into left/right of the break so no line spans the gap."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    e = np.asarray(yerr, dtype=float) if yerr is not None else None
    segs = []
    for mask in (x <= LEFT_XLIM[1], x >= RIGHT_XLIM[0]):
        if mask.any():
            segs.append((x[mask], y[mask], e[mask] if e is not None else None))
    return segs


def bridge_break(bax, x, y, color, **kwargs):
    """Dashed connector between the last point before the break and the first
    point after it. brokenaxes draws it on both panels, so it runs into the gap
    from the left and comes back out of it on the right."""
    if not BRIDGE_BREAK:
        return
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    left = np.flatnonzero(x <= LEFT_XLIM[1])
    right = np.flatnonzero(x >= RIGHT_XLIM[0])
    if not len(left) or not len(right):
        return
    i, j = left[-1], right[0]
    style = dict(BRIDGE_STYLE)
    style.update(kwargs)
    bax.plot(x[[i, j]], y[[i, j]], color=color, label=None, **style)


def _species_names(doc):
    """species_id -> species name."""
    names = {}
    for group in ("small_molecules", "proteins", "complexes"):
        for sp in getattr(doc, group, None) or []:
            if getattr(sp, "name", None):
                names[sp.id] = sp.name
    return names


def _color_map(names, colors):
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    return {n: colors.get(n, cycle[i % len(cycle)]) for i, n in enumerate(names)}


# -----------------------------------------------------------------------------
# 1) ONE PANEL PER MEASUREMENT (redraws pe.plot output on broken axes)
# -----------------------------------------------------------------------------
def to_broken_axes(axs, doc=None, colors=None, title=None,
                   panel_size=(7, 5), close_source=True):
    colors = SPECIES_COLORS if colors is None else colors
    names = _species_names(doc) if doc is not None else {}
    title = title or getattr(doc, "name", None)

    src_axes = [a for a in np.atleast_1d(axs).ravel() if a.get_visible()]
    if not src_axes:
        raise ValueError("No visible axes to convert.")

    _, ncols = src_axes[0].get_subplotspec().get_gridspec().get_geometry()
    nrows = int(np.ceil(len(src_axes) / ncols))

    fig = plt.figure(figsize=(panel_size[0] * ncols, panel_size[1] * nrows))
    gs = GridSpec(nrows, ncols, figure=fig, wspace=0.6, hspace=0.45)

    baxes = []
    for i, src in enumerate(src_axes):
        bax = brokenaxes(xlims=(LEFT_XLIM, RIGHT_XLIM), subplot_spec=gs[i],
                         wspace=BREAK_GAP, fig=fig, diag_color="none")

        for line in src.get_lines():
            x, y = line.get_xdata(), line.get_ydata()
            if len(x) == 0:
                continue
            raw_label = line.get_label()
            is_fit = raw_label.endswith(" Fit")
            sid = raw_label[:-4] if is_fit else raw_label
            name = names.get(sid, sid)
            label = None if raw_label.startswith("_") else (f"{name} Fit" if is_fit else name)
            color = colors.get(name, colors.get(sid, line.get_color()))
            style = dict(color=color,
                         linestyle=line.get_linestyle(),
                         linewidth=line.get_linewidth(),
                         marker=line.get_marker(),
                         markersize=line.get_markersize())

            if line.get_marker() not in (None, "None", "", " "):
                for k, (xs, ys, _) in enumerate(split_by_break(x, y)):
                    bax.plot(xs, ys, label=label if k == 0 else None, **style)
                bridge_break(bax, x, y, color)
            else:   # model fit: also cut at the break so nothing spans the gap
                for k, (xs, ys, _) in enumerate(split_by_break(x, y)):
                    bax.plot(xs, ys, label=label if k == 0 else None, **style)

        bax.set_title(src.get_title(), fontweight="bold")
        _style_bax(bax, src.get_xlabel() or "Time", src.get_ylabel())
        bax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize="small")
        baxes.append(bax)

    if title:
        fig.suptitle(title, fontweight="bold", fontsize=TITLE_FONTSIZE + 2)

    if close_source:
        plt.close(src_axes[0].figure)
    return fig, baxes


# -----------------------------------------------------------------------------
# 2) ALL REPLICATES IN ONE PLOT (or mean +/- SD)
# -----------------------------------------------------------------------------
def _collect_data(doc, species=None):
    """Long table: one row per data point -> Replicate, Species, Time, Value."""
    names = _species_names(doc)
    records = []
    for m in doc.measurements:
        rep = m.name or m.id
        for sd in m.species_data:
            t = list(sd.time or [])
            y = list(sd.data or [])
            if not y or len(t) != len(y):
                continue
            name = names.get(sd.species_id, sd.species_id)
            if species is not None and name not in species and sd.species_id not in species:
                continue
            records += [{"Replicate": rep, "Species": name, "Time": float(ti), "Value": float(yi)}
                        for ti, yi in zip(t, y)]
    return pd.DataFrame(records, columns=["Replicate", "Species", "Time", "Value"])


def plot_replicates(doc, mean_sd=False, species=None, colors=None, connect=True,
                    title=None, xlabel="Time (min)", ylabel="Concentration (mM)",
                    figsize=(8, 5)):
    """mean_sd=False: all replicates overlaid (colour = species, marker = replicate)
       mean_sd=True : mean +/- SD across replicates.  Title defaults to doc.name."""
    colors = SPECIES_COLORS if colors is None else colors
    title = title or getattr(doc, "name", None)

    df = _collect_data(doc, species)
    if df.empty:
        raise ValueError("No measurement data found in the EnzymeML document.")

    sp_names = list(dict.fromkeys(df["Species"]))
    reps = list(dict.fromkeys(df["Replicate"]))
    col = _color_map(sp_names, colors)
    ls = "-" if connect else "none"

    fig = plt.figure(figsize=figsize)
    bax = brokenaxes(xlims=(LEFT_XLIM, RIGHT_XLIM), wspace=BREAK_GAP, fig=fig,
                     diag_color="none")

    if mean_sd:
        for name in sp_names:
            sub = df[df["Species"] == name].copy()
            sub["Time"] = sub["Time"].round(6)
            g = sub.groupby("Time")["Value"]
            mean = g.mean()
            sd = g.std(ddof=1).fillna(0.0)
            for k, (xs, ys, es) in enumerate(split_by_break(mean.index, mean.values, sd.values)):
                bax.errorbar(xs, ys, yerr=es, color=col[name], marker="o", linestyle=ls,
                             linewidth=2, markersize=6, capsize=4,
                             label=name if k == 0 else None)
            if connect:
                bridge_break(bax, mean.index.values, mean.values, col[name], linewidth=2)
        handles, labels = bax.axs[0].get_legend_handles_labels()
    else:
        for (rep, name), sub in df.groupby(["Replicate", "Species"], sort=False):
            sub = sub.sort_values("Time")
            mk = REPLICATE_MARKERS[reps.index(rep) % len(REPLICATE_MARKERS)]
            for xs, ys, _ in split_by_break(sub["Time"], sub["Value"]):
                bax.plot(xs, ys, color=col[name], marker=mk, linestyle=ls,
                         linewidth=1.5, markersize=6, alpha=0.85)
            if connect:
                bridge_break(bax, sub["Time"].values, sub["Value"].values,
                             col[name], alpha=0.85)
        handles = [Line2D([], [], color=col[n], marker="o", linestyle=ls, linewidth=1.5)
                   for n in sp_names]
        handles += [Line2D([], [], color="grey", linestyle="none",
                           marker=REPLICATE_MARKERS[i % len(REPLICATE_MARKERS)])
                    for i in range(len(reps))]
        labels = sp_names + reps

    _style_bax(bax, xlabel, ylabel)
    bax.legend(handles=handles, labels=labels, loc="upper left",
               bbox_to_anchor=(1.01, 1), fontsize="small")
    if title:
        bax.set_title(title, fontweight="bold", fontsize=TITLE_FONTSIZE)
    return fig, bax