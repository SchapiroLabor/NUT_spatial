####################
###COZI DOTPLOT###
####################
#Special dotplot that plots cozi scores

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import TwoSlopeNorm


def cozi_dotplot_facet_split(
    res,
    column_name="cell_pair",
    zscore_col="zscore",
    patient_col="patient_ID",
    split_groups=None,
    figsize=None,
    save_path=None,
    size_scale=3,
    max_dot_size=250,
    min_dot_size=8,
    vmin=None,
    vcenter=0,
    vmax=None,
    survival_data=None,
    survival_col="survival_time_months",
    tumor_site_col="tumor_site",
    survival_ylim=None,
    surv_height=0.18,
    show_legend=True,
    site_width_ratios=None,
    title_fontsize=15,
    tick_fontsize=12,
    group_gap=0.06,
    # manually positioned in figure-fraction coords — tune these to taste
    colorbar_pos=(0.92, 0.35, 0.015, 0.35),  # (left, bottom, width, height)
    legend_pos=(0.91, 0.05),                  # (left, bottom) of dot-size legend box
):
    """
    Like cozi_dotplot_facet_by_site but splits rows into visually separated
    groups (one sub-box per myeloid cell type) with whitespace between them.
    No side labels, no separator lines — just space.
    """

    df = res.copy()
    df[patient_col] = df[patient_col].astype(str)
    df[column_name] = df[column_name].astype(str)
    df[zscore_col] = pd.to_numeric(df[zscore_col], errors="coerce")
    df["cond_percent"] = pd.to_numeric(df["cond_percent"], errors="coerce")

    if split_groups is None:
        raise ValueError("split_groups must be provided")
    if survival_data is None:
        raise ValueError("survival_data is required")

    group_labels = list(split_groups.keys())
    group_pairs  = list(split_groups.values())
    n_groups     = len(group_labels)
    all_pairs    = [p for pairs in group_pairs for p in pairs]

    surv_df = survival_data[[patient_col, survival_col, tumor_site_col]].copy()
    surv_df[patient_col] = surv_df[patient_col].astype(str)

    site_to_patients = (
        surv_df.groupby(tumor_site_col)
        .apply(lambda g: g.sort_values(survival_col)[patient_col].tolist())
        .to_dict()
    )
    unique_sites = list(site_to_patients.keys())
    n_sites      = len(unique_sites)

    if site_width_ratios is not None:
        if len(site_width_ratios) != n_sites:
            raise ValueError(f"site_width_ratios must have {n_sites} values")
        panel_widths = list(site_width_ratios)
    else:
        panel_widths = [len(site_to_patients[s]) for s in unique_sites]

    # ── pre-process df ───────────────────────────────────────────────────────
    patients_all = [p for s in unique_sites for p in site_to_patients[s]]
    full_index = pd.MultiIndex.from_product(
        [patients_all, all_pairs], names=[patient_col, column_name]
    )
    df = (
        df.set_index([patient_col, column_name])
        .reindex(full_index)
        .reset_index()
    )
    df["is_na"] = df[zscore_col].isna()
    df["cond_percent_plot"] = df["cond_percent"].fillna(0)
    df["dot_size"] = np.sqrt(df["cond_percent_plot"].clip(lower=0)) * size_scale * 10
    df["dot_size"] = df["dot_size"].clip(lower=min_dot_size, upper=max_dot_size)
    df.loc[df["is_na"], "dot_size"] = min_dot_size

    # ── color scaling ────────────────────────────────────────────────────────
    valid_z = df.loc[~df["is_na"], zscore_col]
    if vmin is None or vmax is None:
        max_abs = np.nanmax(np.abs(valid_z)) if len(valid_z) > 0 else 1
        if not np.isfinite(max_abs) or max_abs == 0:
            max_abs = 1
        vmin, vmax = -max_abs, max_abs
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    # ── figure ───────────────────────────────────────────────────────────────
    group_sizes   = [len(p) for p in group_pairs]
    total_pairs   = sum(group_sizes)

    if figsize is None:
        figsize = (4 * n_sites + 3, total_pairs * 0.55 + n_groups * 0.4 + 3)

    fig = plt.figure(figsize=figsize)

    # Outer grid: (n_groups dotplot rows + 1 survival row) × (n_sites + 1 colorbar)
    # Height ratios: each group proportional to its pair count, then survival
    dot_heights    = [gs / total_pairs for gs in group_sizes]   # relative
    dot_total      = 1 - surv_height
    height_ratios  = [gs * dot_total for gs in dot_heights] + [surv_height]

    outer_gs = GridSpec(
        nrows=n_groups + 1,
        ncols=n_sites,
        width_ratios=panel_widths,
        height_ratios=height_ratios,
        wspace=0.12,
        hspace=group_gap,
        figure=fig,
        left=0.18,    # leave room for y-tick labels
        right=0.88,   # leave room for colorbar + legend on the right
    )

    # ── manually placed colorbar and legend axes ─────────────────────────────
    cax = fig.add_axes(colorbar_pos)   # [left, bottom, width, height]

    lax = fig.add_axes([legend_pos[0], legend_pos[1], 0.12, 0.22])
    lax.axis("off")

    sc = None

    # ── loop over groups (rows) ──────────────────────────────────────────────
    for g_idx, (group_label, pairs) in enumerate(split_groups.items()):
        n_pairs  = len(pairs)
        pair_to_y = {p: i for i, p in enumerate(pairs)}

        group_df = df[df[column_name].isin(pairs)].copy()
        group_df["y"] = group_df[column_name].map(pair_to_y)

        # ── loop over sites (cols) ───────────────────────────────────────────
        for i, site in enumerate(unique_sites):
            site_patients = site_to_patients[site]
            ax = fig.add_subplot(outer_gs[g_idx, i])

            sub_df = group_df[group_df[patient_col].isin(site_patients)].copy()
            local_x = {p: j for j, p in enumerate(site_patients)}
            sub_df["x"] = sub_df[patient_col].map(local_x)

            sc = ax.scatter(
                sub_df.loc[~sub_df["is_na"], "x"],
                sub_df.loc[~sub_df["is_na"], "y"],
                c=sub_df.loc[~sub_df["is_na"], zscore_col],
                s=sub_df.loc[~sub_df["is_na"], "dot_size"],
                cmap="coolwarm",
                norm=norm,
                edgecolor="black",
                linewidth=0.25,
                zorder=3,
            )
            ax.scatter(
                sub_df.loc[sub_df["is_na"], "x"],
                sub_df.loc[sub_df["is_na"], "y"],
                s=sub_df.loc[sub_df["is_na"], "dot_size"],
                color="darkgrey",
                edgecolor="black",
                linewidth=0.25,
                zorder=3,
            )

            # site title only on first group row
            if g_idx == 0:
                ax.set_title(site, fontsize=title_fontsize, pad=6)

            ax.set_xlim(-0.5, len(site_patients) - 0.5)
            ax.set_ylim(-0.5, n_pairs - 0.5)
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.tick_params(axis="x", labelbottom=False)
            ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.5, zorder=0)

            # y-axis labels only on leftmost site
            if i == 0:
                ax.set_yticks(range(n_pairs))
                ax.set_yticklabels(pairs, fontsize=tick_fontsize)
                if g_idx == n_groups // 2:
                    ax.set_ylabel("Cell pair", fontsize=title_fontsize, labelpad=8)
            else:
                ax.set_yticks([])

    for i, site in enumerate(unique_sites):
        site_patients = site_to_patients[site]
        sax = fig.add_subplot(outer_gs[-1, i])

        surv_sub = surv_df[surv_df[patient_col].isin(site_patients)]
        survival_values = (
            surv_sub.set_index(patient_col)
            .reindex(site_patients)[survival_col]
        )
        sax.bar(
            range(len(site_patients)),
            survival_values,
            color="grey",
            edgecolor="black",
            linewidth=0.3,
        )
        sax.set_xticks(range(len(site_patients)))
        sax.set_xticklabels(site_patients, rotation=90, fontsize=tick_fontsize)

        if i == 0:
            sax.set_ylabel("Survival\n(months)", fontsize=10)
        else:
            sax.set_yticks([])
            sax.spines["left"].set_visible(False)

        if survival_ylim is not None:
            sax.set_ylim(survival_ylim)

        sax.spines["top"].set_visible(False)
        sax.spines["right"].set_visible(False)

    #colorbar
    if sc is not None:
        cbar = fig.colorbar(sc, cax=cax)
        cbar.set_label("Neighbor Preference z-score", fontsize=12, labelpad=10)
        cbar.set_ticks([vmin, vcenter, vmax])
        cbar.ax.tick_params(labelsize=11)

    #dot size legend
    if show_legend:
        legend_vals = [5, 25, 50, 75]
        handles = []
        for val in legend_vals:
            size = np.clip(np.sqrt(val) * size_scale * 3, min_dot_size, max_dot_size)
            handles.append(
                lax.scatter([], [], s=size, color="#B30326", edgecolor="black", linewidth=0.25)
            )
        na_handle = lax.scatter([], [], s=min_dot_size, color="darkgrey", edgecolor="black")
        lax.legend(
            handles + [na_handle],
            [f"{v}%" for v in legend_vals] + ["Low cell count"],
            title="Conditional cell\npercentage",
            frameon=False,
            loc="upper left",
            fontsize=11,
            title_fontsize=12,
        )

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    plt.show()
    return fig