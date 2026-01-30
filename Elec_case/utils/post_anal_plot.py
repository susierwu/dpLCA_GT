import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import pickle
import matplotlib.colors as mc
from matplotlib.ticker import FuncFormatter



# ============================================================
# Global SSP display mapping (used in ALL plots)
# ============================================================
SSP_LABEL_MAP = {
    "SSP1-VLLO": "SSP1-19",
    "SSP2-M":    "SSP2-45",
    "SSP5-H":    "SSP5-85",
}

def format_ssp(ssp: str) -> str:
    """Return standardized display name for SSP."""
    return SSP_LABEL_MAP.get(ssp, ssp)




def lighten(color, amount=0.5):
    """Lighten color by mixing with white."""
    try:
        c = mc.cnames[color]
    except KeyError:
        c = color
    r, g, b = mc.to_rgb(c)
    return (1 - amount) + amount * r, (1 - amount) + amount * g, (1 - amount) + amount * b





# step1: slopegraph, absolute rank shift 


def rank_shift_plot_cleanSeparated(
    df_long,
    metric_order,
    metric_labels,
    tech_colors,
    fig_title="Rank shifts across GWP metrics",
    years=(2030, 2050),
    countries=("CN", "US"),
    ssps=("SSP1-VLLO", "SSP2-M", "SSP5-H"),
    lw=3.0,
    sep=0.03,
    anot_off=0.16
):
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    df = df_long.copy()
    df["Rank"] = (
        df.groupby(["Country", "SSP", "Year", "Metric"])["Score"]
          .rank(method="dense", ascending=True)
    )

    dotted_2050 = (0, (1, 3))
    year_ls    = {2030: "solid", 2050: dotted_2050}
    year_alpha = {2030: 1.0, 2050: 0.9}
    zmap = {2030: 5, 2050: 4}

    fig, axes = plt.subplots(
        nrows=len(ssps),
        ncols=len(countries),
        figsize=(12, 12),   # 2nd parameter 10 to 12 -> taller figure
        sharex=True,
        sharey=True,
    )
    plt.subplots_adjust(
        top=0.90,
        bottom=0.18,
        hspace=0.55,         # MORE vertical breathing room
        wspace=0.20,
    )

      
    plt.subplots_adjust(top=0.88, bottom=0.20, hspace=0.35, wspace=0.18)

    x = np.arange(len(metric_order))

    for i, ssp in enumerate(ssps):
        for j, ct in enumerate(countries):
            ax = axes[i, j]
            sub = df[(df["SSP"] == ssp) & (df["Country"] == ct)]

            for tech in ["Hydro", "PV", "Wind"]:
                sub_t = sub[sub["Technology"] == tech]

                for yr in years:
                    tdf = sub_t[sub_t["Year"] == yr]
                    if tdf.empty:
                        continue

                    y_rank = tdf.set_index("Metric")["Rank"].reindex(metric_order)
                    y_score = tdf.set_index("Metric")["Score"].reindex(metric_order)

                    shift     = sep if yr == 2030 else -sep
                    y_shifted = y_rank + shift

                    ax.plot(
                        x,
                        y_shifted,
                        linestyle=year_ls[yr],
                        color=tech_colors[tech],
                        alpha=year_alpha[yr],
                        linewidth=lw,
                        marker="x",
                        markersize=7,
                        zorder=zmap[yr],
                    )

                    # ---- BOXED ANNOTATIONS ----
                    for xi, yi, sc in zip(x, y_shifted, y_score):
                        if pd.isna(sc):
                            continue

                        val = sc * 1000.0
                        y_text = yi + anot_off if yr == 2030 else yi - anot_off
                        va     = "bottom" if yr == 2030 else "top"

                        ax.text(
                            xi,
                            y_text,
                            f"{val:.1f}",
                            ha="center",
                            va=va,
                            fontsize=9.5,
                            color="white",
                            bbox=dict(
                                    boxstyle="round,pad=0.18",      # slightly tighter box
                                    facecolor=tech_colors[tech],
                                    edgecolor="black",
                                    linewidth=0.9,
                                    linestyle="solid" if yr == 2030 else "dashed",
                                    alpha=0.8,                     # lets grid & lines breathe
                                ),
                            zorder=10,
                        )

            #ax.set_title(f"{ct}, {ssp}", fontsize=14, pad=8)
            ax.set_title(f"{ct}, {format_ssp(ssp)}", fontsize=13)
            
            ax.set_xticks(x)

            if i == len(ssps) - 1:
                ax.set_xticklabels(
                    [metric_labels[m] for m in metric_order],
                    rotation=25,
                    ha="right",
                    fontsize=12,
                )
            else:
                ax.set_xticklabels([])

            ax.set_yticks([1, 2, 3])
            ax.set_ylim(0.2, 3.8)
            if j == 0:
                ax.set_ylabel("Rank (1 = lowest GWP)", fontsize=13)

            ax.grid(axis="y", linestyle=":", alpha=0.3)

    # Legend
    handles, labels = [], []
    for tech in ["Hydro", "PV", "Wind"]:
        for yr in years:
            h, = axes[0, 0].plot(
                [],
                [],
                linestyle=year_ls[yr],
                color=tech_colors[tech],
                linewidth=lw,
                marker="x",
                markersize=7,
            )
            handles.append(h)
            labels.append(f"{tech} {yr}")

    fig.legend(
        handles,
        labels,
        title="Technology / Year",
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=6,
        fontsize=12,
        frameon=True,
    )

    fig.suptitle(fig_title, fontsize=16, y=0.97)
    return fig, axes




# step2a: normalized barplots just one year 2030 or 2050, not mixing in one plot
### this stackedbar is not favored by Cecile so not using it 
def plot_stacked_normalized_bars(
    df_full,
    metric_order,
    metric_labels,
    countries=("CN", "US"),
    ssps=("SSP1-VLLO", "SSP2-M", "SSP5-H"),
    technologies=("Hydro", "Wind", "PV"),
    tech_colors=None,
    bar_width=0.23,
):
    if tech_colors is None:
        tech_colors = {"Hydro": "#7B7B7B", "Wind": "#006BA4", "PV": "#FF800E"}

    fig, axes = plt.subplots(
        nrows=len(ssps),
        ncols=len(countries),
        figsize=(18, 15),
        sharey=True,
    )
    plt.subplots_adjust(hspace=0.35, wspace=0.22)

    YPAD = 2  # annotation spacing

    for i, ssp in enumerate(ssps):
        for j, country in enumerate(countries):
            ax = axes[i, j]

            sub = df_full[(df_full["SSP"] == ssp) & (df_full["Country"] == country)]
            x_positions = np.arange(len(metric_order))
            offset = np.linspace(-bar_width, bar_width, len(technologies))

            # ------------------------------
            # Shading for dpLCI metrics
            # ------------------------------
            for k, metric in enumerate(metric_order):
                if metric in [
                    "dyn-FG-LCI + pGWP-fixedCO2",
                    "full dpLCI + pGWP-fixedCO2",
                ]:
                    ax.axvspan(k - 0.5, k + 0.5, color="#CCCCCC", alpha=0.15, zorder=0)

            # ------------------------------
            # Technology bars
            # ------------------------------
            for t_index, tech in enumerate(technologies):
                t_df = sub[sub["Technology"] == tech]

                norm2030, norm2050 = [], []
                raw2030, raw2050 = [], []

                for metric in metric_order:
                    sc2030 = t_df[(t_df["Metric"] == metric) & (t_df["Year"] == 2030)]["Score"]
                    sc2050 = t_df[(t_df["Metric"] == metric) & (t_df["Year"] == 2050)]["Score"]

                    s2030 = sc2030.values[0] if len(sc2030) else np.nan
                    s2050 = sc2050.values[0] if len(sc2050) else np.nan

                    # raw values (*1000)
                    raw2030.append(s2030 * 1000 if not np.isnan(s2030) else np.nan)
                    raw2050.append(s2050 * 1000 if not np.isnan(s2050) else np.nan)

                    # normalization
                    max2030 = sub[(sub["Metric"] == metric) & (sub["Year"] == 2030)]["Score"].max()
                    max2050 = sub[(sub["Metric"] == metric) & (sub["Year"] == 2050)]["Score"].max()

                    n2030 = (s2030 / max2030 * 100) if max2030 else np.nan
                    n2050 = (s2050 / max2050 * 100) if max2050 else np.nan

                    norm2030.append(n2030)
                    norm2050.append(n2050)

                norm2030 = np.array(norm2030)
                norm2050 = np.array(norm2050)
                top_seg = norm2030 - norm2050

                xpos = x_positions + offset[t_index]

                # 2050 base bar
                ax.bar(
                    xpos,
                    norm2050,
                    width=bar_width,
                    color=tech_colors[tech],
                    edgecolor="black",
                    label=f"{tech} 2050" if (i == 0 and j == 0) else None,
                )

                # 2030 top segment
                ax.bar(
                    xpos,
                    top_seg,
                    bottom=norm2050,
                    width=bar_width,
                    color=lighten(tech_colors[tech], 0.40),
                    edgecolor="black",
                    label=f"{tech} 2030" if (i == 0 and j == 0) else None,
                )

                # ------------------------------
                # Annotation placement rules
                # ------------------------------
                for xi, base, top, r2030, r2050 in zip(xpos, norm2050, top_seg, raw2030, raw2050):
                
                    if tech == "Hydro":
                        # --- HYDRO SPECIAL CASE ---
                
                        # 2050 label: always near x-axis (bottom)
                        ax.text(
                            xi,
                            2.0,                     # fixed small height above axis
                            f"{r2050:.1f}",
                            ha="center",
                            va="bottom",
                            fontsize=11,
                            color="white",
                        )
                
                        # 2030 label: always on top of the stacked bar
                        ax.text(
                            xi,
                            base + top + YPAD,
                            f"{r2030:.1f}",
                            ha="center",
                            va="bottom",
                            fontsize=11,
                            color="black",
                        )
                
                    else:
                        # --- WIND & PV: previous smooth logic ---
                        # 2050 annotation
                        if base > 30:
                            # inside bar (white text)
                            ax.text(
                                xi,
                                base * 0.55,
                                f"{r2050:.1f}",
                                ha="center",
                                va="center",
                                fontsize=11,
                                color="white",
                            )
                        else:
                            # above bar
                            ax.text(
                                xi,
                                base + 1.0,
                                f"{r2050:.1f}",
                                ha="center",
                                va="bottom",
                                fontsize=11,
                                color="black",
                            )
                
                        # then ALL 2030 annotation (above stacked segment)
                        ax.text(
                            xi,
                            base + top + YPAD,
                            f"{r2030:.2f}",
                            ha="center",
                            va="bottom",
                            fontsize=10,
                            color="black",
                        )

            # ------------------------------
            # Layout: remove xticks except bottom row
            # ------------------------------
            ax.set_title(f"{country}, {ssp}", fontsize=14, pad=8)
            ax.grid(axis="y", linestyle=":", alpha=0.35)
            ax.set_ylim(0, 125)

            if i == len(ssps) - 1:
                ax.set_xticks(x_positions)
                ax.set_xticklabels(
                    [metric_labels[m] for m in metric_order],
                    rotation=28, ha="right"
                )
            else:
                ax.set_xticks([])
                ax.set_xticklabels([])

    # Legend
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        title="Technology / Year",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=6, fontsize=12,
        frameon=True
    )

    fig.suptitle("Normalized (0–100) bars for each technology among different metrics \n the 2030 GWP score (raw GWP score in black) is stacked on top of the 2050 base bar (raw GWP score in white)", fontsize=18, y=0.95)
    return fig, axes



def plot_bar_oneyear_only(
    df_full,
    year,
    metric_order,
    metric_labels,
    tech_colors,
    countries=("CN", "US"),
    ssps=("SSP1-VLLO", "SSP2-M", "SSP5-H"),
    technologies=("Hydro", "Wind", "PV"),
    bar_width=0.25
):
    """
    PURE ranking per year.
    Normalize panel-by-panel using max score of that panel for that YEAR ONLY.
    Annotate bars using raw GWP × 1000.
    Y-axis fixed to 0–100 (%).
    """

    df_y = df_full[df_full["Year"] == year].copy()

    fig, axes = plt.subplots(
        nrows=len(ssps),
        ncols=len(countries),
        figsize=(14, 10),
        sharex=False,
        sharey=True   # shared since now fixed 0–100
    )

    plt.subplots_adjust(hspace=0.28, wspace=0.18)

    for i, ssp in enumerate(ssps):
        for j, ct in enumerate(countries):

            ax = axes[i, j]
            sub = df_y[(df_y["SSP"] == ssp) & (df_y["Country"] == ct)]

            panel_max = sub.groupby("Metric")["Score"].max().max()
            if pd.isna(panel_max) or panel_max == 0:
                continue

            x = np.arange(len(metric_order))

            for t_i, tech in enumerate(technologies):

                sub_t = sub[sub["Technology"] == tech]
                base_color = tech_colors[tech]

                for k, metric in enumerate(metric_order):

                    df_val = sub_t[sub_t["Metric"] == metric]["Score"]
                    if len(df_val) == 0:
                        continue

                    v = float(df_val.values[0])
                    normalized = v / panel_max * 100
                    raw_annot = v * 1000
                    xpos = k + (t_i - 1) * bar_width

                    ax.bar(
                        xpos,
                        normalized,
                        width=bar_width,
                        color=base_color,
                        edgecolor="black",
                        linewidth=0.7
                    )

                    if normalized > 35:
                        ax.text(
                            xpos, normalized * 0.55,
                            f"{raw_annot:.1f}",
                            ha="center", va="center",
                            fontsize=10, color="white"
                        )
                    else:
                        ax.text(
                            xpos, normalized + 2,
                            f"{raw_annot:.1f}",
                            ha="center", va="bottom",
                            fontsize=10, color="black"
                        )

            # --------------------------------------------------
            # Titles & axes
            # --------------------------------------------------
            ax.set_title(f"{ct}, {format_ssp(ssp)}", fontsize=13)

            #ax.set_ylim(0, 100)
            ### updated adding % 
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.yaxis.set_major_formatter(
                FuncFormatter(lambda x, pos: f"{int(x)}%")
            )
            
            if j == 0:
                ax.set_ylabel("Normalized impact (%)", fontsize=11)

            if i == len(ssps) - 1:
                ax.set_xticks(np.arange(len(metric_order)))
                ax.set_xticklabels(
                    [metric_labels[m] for m in metric_order],
                    rotation=25, ha="right", fontsize=10
                )
            else:
                ax.set_xticks([])

            # Stronger, clearer gridlines
            ax.grid(
                axis="y",
                linestyle="--",
                linewidth=0.9,
                alpha=0.8,
                color="0.35"
            )

    fig.suptitle(
        f"{year} technology GWP scores (normalized to each scenario max)",
        fontsize=16,
        y=0.96
    )

    plt.tight_layout(rect=[0, 0.07, 1, 0.95])

    # --------------------------------------------------
    # Legend
    # --------------------------------------------------
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=tech_colors["Hydro"], ec="black"),
        plt.Rectangle((0, 0), 1, 1, color=tech_colors["Wind"], ec="black"),
        plt.Rectangle((0, 0), 1, 1, color=tech_colors["PV"], ec="black"),
    ]
    labels = ["Hydro", "Wind", "PV"]

    fig.legend(
        handles, labels,
        loc="lower center",
        ncol=3,
        fontsize=12,
        bbox_to_anchor=(0.5, 0)
    )

    return fig, axes

    

# step2: one region two MY as subplots
def plot_bar_oneregion_twoyears(
    df_full,
    years,
    country,
    metric_order,
    metric_labels,
    tech_colors,
    ssps=("SSP1-VLLO", "SSP2-M", "SSP5-H"),
    technologies=("Wind", "Hydro", "PV"),
    bar_width=0.25,
):
    """
    PURE ranking per region.
    Same REGION, two YEARS as columns.
    Normalize panel-by-panel using max score of that (SSP, YEAR) panel.
    Annotate bars using raw GWP × 1000.
    Y-axis fixed to 0–100 (%).
    """

    df_r = df_full[df_full["Country"] == country].copy()

    fig, axes = plt.subplots(
        nrows=len(ssps),
        ncols=len(years),
        figsize=(14, 10),
        sharex=False,
        sharey=True   # shared since fixed to 0–100
    )

    plt.subplots_adjust(hspace=0.28, wspace=0.18)

    for i, ssp in enumerate(ssps):
        for j, year in enumerate(years):

            ax = axes[i, j]
            sub = df_r[(df_r["SSP"] == ssp) & (df_r["Year"] == year)]

            panel_max = sub.groupby("Metric")["Score"].max().max()
            if pd.isna(panel_max) or panel_max == 0:
                continue

            x = np.arange(len(metric_order))

            for t_i, tech in enumerate(technologies):

                sub_t = sub[sub["Technology"] == tech]
                base_color = tech_colors[tech]

                for k, metric in enumerate(metric_order):

                    df_val = sub_t[sub_t["Metric"] == metric]["Score"]
                    if len(df_val) == 0:
                        continue

                    v = float(df_val.values[0])
                    normalized = v / panel_max * 100
                    raw_annot = v * 1000
                    xpos = k + (t_i - 1) * bar_width

                    ax.bar(
                        xpos,
                        normalized,
                        width=bar_width,
                        color=base_color,
                        edgecolor="black",
                        linewidth=0.7
                    )

                    # Always annotate inside bar
                    ax.text(
                        xpos,
                        normalized * 0.55,
                        f"{raw_annot:.1f}",
                        ha="center",
                        va="center",
                        fontsize=12,
                        color="white"
                    )

            # --------------------------------------------------
            # Titles & axes
            # --------------------------------------------------
            ax.set_title(
                f"{country}, {format_ssp(ssp)}, {year}",
                fontsize=13
            )

            #ax.set_ylim(0, 100)
            ### updated adding % 
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.yaxis.set_major_formatter(
                FuncFormatter(lambda x, pos: f"{int(x)}%")
            )
            

            if j == 0:
                ax.set_ylabel("Normalized impact (%)", fontsize=11)

            if i == len(ssps) - 1:
                ax.set_xticks(np.arange(len(metric_order)))
                ax.set_xticklabels(
                    [metric_labels[m] for m in metric_order],
                    rotation=25,
                    ha="right",
                    fontsize=10
                )
            else:
                ax.set_xticks([])

            # Stronger, clearer gridlines
            ax.grid(
                axis="y",
                linestyle="--",
                linewidth=0.9,
                alpha=0.8,
                color="0.35"
            )

    fig.suptitle(
        f"{country} technology GWP scores "
        f"(normalized per year and scenario)", #; labels show absolute scores ×1000
        fontsize=16,
        y=0.96
    )

    plt.tight_layout(rect=[0, 0.07, 1, 0.95])

    # --------------------------------------------------
    # Legend
    # --------------------------------------------------
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=tech_colors[t], ec="black")
        for t in technologies
    ]

    fig.legend(
        handles,
        technologies,
        loc="lower center",
        ncol=len(technologies),
        fontsize=12,
        bbox_to_anchor=(0.5, 0)
    )

    return fig, axes









#step3: pair-wise tech comp DOTS 
def plot_pairwise_gaps_grid_directional(
    df_long,
    country,
    ssps=("SSP1-VLLO", "SSP2-M", "SSP5-H"),
    years=(2030, 2050),
    metric_order=None,
    metric_labels=None,
    pairs=None,
    gap_threshold_pct=5.0,
    figsize=(14, 10),
    size_scale=45,
    pad_ylim_adjust = 90, 
    directional=True,
):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    if metric_labels is None:
        metric_labels = {m: m for m in metric_order}

    if pairs is None:
        raise ValueError("You must provide `pairs`, e.g. (('PV','Wind'), ('Hydro','PV'))")

    # ------------------------------------------------------------------
    # Dynamic, order-safe styling from `pairs`
    # ------------------------------------------------------------------
    marker_cycle = ["o", "s", "D", "^", "v", "P", "X"]
    gray_cycle   = ["0.25", "0.45", "0.65", "0.35", "0.55", "0.75"]

    marker_map = {}
    gray_map   = {}

    for i, (a, b) in enumerate(pairs):
        label = f"{a} vs {b}"
        marker_map[label] = marker_cycle[i % len(marker_cycle)]
        gray_map[label]   = gray_cycle[i % len(gray_cycle)]

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        nrows=len(ssps),
        ncols=len(years),
        figsize=figsize,
        sharey=True,
    )
    plt.subplots_adjust(hspace=0.35, wspace=0.20)

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------
    for i, ssp in enumerate(ssps):
        for j, year in enumerate(years):

            ax = axes[i, j]

            sub = df_long[
                (df_long["Country"] == country) &
                (df_long["SSP"] == ssp) &
                (df_long["Year"] == year)
            ]

            if sub.empty:
                ax.set_axis_off()
                continue

            x = np.arange(len(metric_order))

            for (a, b) in pairs:
                label = f"{a} vs {b}"
                gaps = []

                for m in metric_order:
                    sm = sub[sub["Metric"] == m].set_index("Technology")["Score"]

                    if a in sm.index and b in sm.index:
                        sa, sb = float(sm.loc[a]), float(sm.loc[b])
                        mag = abs(sa - sb) / max(sa, sb) * 100.0

                        if directional:
                            gap = +mag if sa < sb else -mag
                        else:
                            gap = mag
                        gaps.append(gap)
                    else:
                        gaps.append(np.nan)

                gaps = np.array(gaps, dtype=float)

                ax.scatter(
                    x,
                    gaps,
                    s=size_scale * np.sqrt(np.abs(gaps)),
                    marker=marker_map[label],
                    facecolor=gray_map[label],
                    edgecolor="black",
                    linewidth=0.8,
                    alpha=0.9,
                    zorder=5,
                    label=label if (i == 0 and j == 0) else None,
                )

            # ------------------------------------------------------------------
            # Reference lines
            # ------------------------------------------------------------------
            ax.axhline(0, color="black", linewidth=0.9)
            ax.axhline(+gap_threshold_pct, linestyle="--", color="red", alpha=0.7)
            ax.axhline(-gap_threshold_pct, linestyle="--", color="red", alpha=0.7)

            # ------------------------------------------------------------------
            # Formatting
            # ------------------------------------------------------------------
            #ax.set_title(f"{country}, {ssp}, {year}", fontsize=12)
            ax.set_title(
                f"{country}, {format_ssp(ssp)}, {year}",
                fontsize=12
            )
            
            ax.set_xticks(x)

            if i == len(ssps) - 1:
                ax.set_xticklabels(
                    [metric_labels[m] for m in metric_order],
                    rotation=20,
                    ha="right",
                    fontsize=10,
                )
            else:
                ax.set_xticklabels([])

            if j == 0:
                ax.set_ylabel(
                    "Directional pairwise gap (%)\n(+ = first tech lower GWP)",
                    fontsize=11,
                )

            #ax.grid(axis="y", linestyle=":", alpha=0.35)
            # NEW Y-axis ticks and labels (fixed, symmetric, percentage)
            #ax.set_yticks([-60, -35, 0, 30, 60])
            #ax.set_yticklabels(["-60%", "-30%", "0%", "30%", "60%"], fontsize=12)
            ax.set_yticks([-50, 0, 50])
            ax.set_yticklabels(["-50%", "0%", "50%"], fontsize=12)
            
            # Darker, clearer gridlines
            ax.grid(
                axis="y",
                linestyle=":",
                linewidth=0.9,
                color="0.5",
                alpha=0.9
            )


    # ---------------------------------------------------
    # Fix y-limits to account for marker size
    # ---------------------------------------------------
    all_gaps = df_long[
        (df_long["Country"] == country) &
        (df_long["SSP"].isin(ssps)) &
        (df_long["Year"].isin(years))
    ]["Score"]
    
    # Conservative padding (in % points)
    ymin = df_long["Score"].min()
    ymax = df_long["Score"].max()
    pad = pad_ylim_adjust  # % points , adjusting here, to make sure all markers are shown, 90% here becaz the ymax/min is about 75% 
    for ax in axes.flat:
        ax.set_ylim(ymin - pad, ymax + pad)

    # above conservative padding (in % points) too much ylim, caz we know it's already close to 60%
    #for ax in axes.flat:
    #    ax.set_ylim(-60, 60)

    
    # ------------------------------------------------------------------
    # Legend (dynamic, directional, marker-preserving)
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    legend_handles = []
    
    for (a, b) in pairs:
        label = f"{a} vs {b}"
    
        legend_label = (
            f"+ : {a} lower impact than {b} (by Δ %)\n"
            f"− : {a} higher impact than {b} (by Δ %)\n"
            f"Δ = |{a} − {b}| / max({a}, {b}) × 100"
        )
    
        handle = Line2D(
            [0], [0],
            marker=marker_map[label],      # marker still encodes the pair
            linestyle="None",
            markerfacecolor=gray_map[label],
            markeredgecolor="black",
            markeredgewidth=0.8,
            markersize=9,
            label=legend_label,
        )
    
        legend_handles.append(handle)
    
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=len(pairs),
        fontsize=11,
        frameon=True,
        bbox_to_anchor=(0.5, -0.08),
    )


    fig.suptitle(
        f"Pairwise technology gaps across metrics \n red lines show the 5% gap",
        fontsize=15,
        y=0.96,
    )

    return fig, axes