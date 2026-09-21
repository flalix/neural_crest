"""Pseudotime and fate figures for a fitted scFates object.

Why not scf.pl.trajectory / scf.pl.trends directly: both draw the principal
graph over whatever `basis` you pass, so a graph fitted in harmony/PCA space
appears as a zigzag over a UMAP, and their gene labels collide past ~15
features. This module keeps the embedding panels label-light and draws the
graph only when asked for the space it was actually fitted in.

"Fate" in scFates is not a probability: a cell's fate is the terminal branch
of the principal graph it sits on (`obs['seg']`, `obs['milestones']`). For
probabilistic fate assignment use CellRank, which scFates can import via
`scf.tl.cellrank_to_tree`.

    from libs.plot_fates import plot_pseudotime_fates
    fig = plot_pseudotime_fates(IN, "results/figures/pseudotime_fates.png",
                                label_key="cell.types", basis="X_umap_subset")
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def terminal_fates(adata):
    """Map each cell to the terminal branch it belongs to.

    Returns (labels, info) where labels is a Series over obs_names naming the
    fate ("fate: <dominant label>" for terminal segments, "progenitor//internal"
    otherwise), and info lists per-segment size and whether it is terminal.
    """
    g = adata.uns["graph"]
    tips = {int(t) for t in np.ravel(g["tips"])}
    seg = adata.obs["seg"].astype(str)
    rows, naming = [], {}
    for _, r in g["pp_seg"].iterrows():
        s = str(int(r["n"]))
        is_term = int(r["to"]) in tips
        members = adata.obs.loc[seg == s]
        dom = (members["cell.types"].astype(str).value_counts().idxmax()
               if len(members) and "cell.types" in members else "n/a")
        naming[s] = f"fate: {dom}" if is_term else "transit"
        rows.append(dict(seg=s, n=int((seg == s).sum()), terminal=is_term,
                         dominant=dom, length=float(r["d"])))
    return seg.map(naming).rename("fate"), pd.DataFrame(rows)


def plot_pseudotime_fates(adata, out_png, label_key="cell.types",
                          basis="X_umap_subset", t_key="t", draw_graph=False,
                          title=None, sizes=(8, 7, 6), figsize=(7.2, 4.6)):
    """Four panels: cell types, pseudotime, fate/branch, pseudotime per type.

    Two independent knobs, easy to confuse:

    figsize : (width, height) in INCHES -- the canvas, passed to plt.figure
    sizes   : (base, secondary, tick) in POINTS -- the font ladder
              (titles/axis-labels, legend/annotation, tick labels)

    Enlarging figsize alone makes text relatively smaller, so scale the
    ladder with it: at figsize=(12, 8) use sizes=(13, 12, 10). Marker areas
    are scaled automatically from the width.

    `draw_graph=True` overlays the principal points -- only meaningful when
    `basis` is the representation the graph was fitted in (uns['ppt']/['epg']
    coordinates live in that space, not in a UMAP).
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    k = figsize[0] / 7.2          # scale markers off the reference width
    pt = 3 * k ** 1.6
    fs_base, fs_small, fs_tick = sizes

    try:  # figure-style skill, when loaded
        apply_figure_style(sizes=sizes)  # noqa: F821
        frame = set_frame  # noqa: F821
        letter = panel_letter  # noqa: F821
    except NameError:
        mpl.rcParams.update({"font.size": fs_small, "axes.titlesize": fs_base,
                             "axes.labelsize": fs_base, "xtick.labelsize": fs_tick,
                             "ytick.labelsize": fs_tick, "figure.dpi": 110})
        frame = lambda ax: [ax.spines[s].set_visible(False) for s in ("top", "right")]  # noqa: E731
        letter = lambda ax, l: ax.set_title(l, loc="left", fontweight="bold")  # noqa: E731, E741

    emb = np.asarray(adata.obsm[basis])[:, :2]
    lab = adata.obs[label_key].astype(str)
    t = adata.obs[t_key].values
    fate, seg_info = terminal_fates(adata)

    cats = sorted(lab.unique())
    pal = dict(zip(cats, plt.get_cmap("tab20")(np.linspace(0, 1, max(len(cats), 2)))))
    fates = sorted(fate.unique())
    fpal = dict(zip(fates, plt.get_cmap("Set2")(np.linspace(0, 1, max(len(fates), 2)))))

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1.0], hspace=0.45, wspace=0.12)

    def blank(ax):
        ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.05)
        for s in ax.spines.values():
            s.set_visible(False)

    ax0 = fig.add_subplot(gs[0, 0])
    for c in cats:
        m = (lab == c).values
        ax0.scatter(emb[m, 0], emb[m, 1], s=pt, c=[pal[c]], lw=0, rasterized=True, label=c)
    ax0.set_title("cell type", fontsize=fs_base, loc="left"); blank(ax0); letter(ax0, "a")

    ax1 = fig.add_subplot(gs[0, 1])
    sc1 = ax1.scatter(emb[:, 0], emb[:, 1], s=pt, c=t, cmap="viridis", lw=0, rasterized=True)
    if draw_graph and "graph" in adata.uns:
        F = np.asarray(adata.uns["graph"]["F"])[:2].T
        B = np.asarray(adata.uns["graph"]["B"])
        for i, j in zip(*np.where(np.triu(B) > 0)):
            ax1.plot(F[[i, j], 0], F[[i, j], 1], c="k", lw=0.6, zorder=3)
    cax = ax1.inset_axes([0.62, 0.04, 0.34, 0.035])
    cb = fig.colorbar(sc1, cax=cax, orientation="horizontal")
    cb.set_label("pseudotime", fontsize=fs_tick, labelpad=1); cb.ax.tick_params(labelsize=fs_tick, pad=1)
    ax1.set_title("pseudotime", fontsize=fs_base, loc="left"); blank(ax1); letter(ax1, "b")

    ax2 = fig.add_subplot(gs[0, 2])
    for f in fates:
        m = (fate == f).values
        ax2.scatter(emb[m, 0], emb[m, 1], s=pt, c=[fpal[f]], lw=0, rasterized=True, label=f)
    ax2.set_title(f"branch ({int(seg_info.terminal.sum())} terminal)", fontsize=fs_base, loc="left")
    blank(ax2); letter(ax2, "c")
    ax2.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False,
               fontsize=fs_tick, markerscale=2.2, handletextpad=0.3, labelspacing=0.25)

    ax3 = fig.add_subplot(gs[1, :])
    order = adata.obs.groupby(label_key, observed=True)[t_key].mean().sort_values().index
    rng = np.random.default_rng(0)
    for i, c in enumerate(order):
        v = t[(lab == str(c)).values]
        ax3.scatter(v, i + rng.uniform(-0.22, 0.22, v.size), s=pt,
                    c=[pal[str(c)]], lw=0, rasterized=True)
        ax3.plot([np.median(v)] * 2, [i - 0.34, i + 0.34], c="0.15", lw=1.2, zorder=4)
    ax3.set_yticks(range(len(order))); ax3.set_yticklabels([str(c) for c in order], fontsize=fs_tick)
    ax3.set_xlabel("pseudotime (bar = median)", fontsize=fs_base)
    ax3.set_title("ordering of cell types along the trajectory", fontsize=fs_base, loc="left")
    ax3.margins(x=0.02, y=0.02); frame(ax3); letter(ax3, "d")

    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=3.5 * k, color=pal[c], label=c)
                        for c in cats],
               loc="upper left", bbox_to_anchor=(0.005, 0.055), ncol=min(len(cats), 5),
               frameon=False, fontsize=fs_tick, handletextpad=0.25, columnspacing=0.9)
    if title:
        fig.suptitle(title, fontsize=round(fs_base * 1.15, 1), x=0.005, ha="left")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    return fig, seg_info
