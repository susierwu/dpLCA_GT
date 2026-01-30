import pandas as pd
import matplotlib.pyplot as plt
import os 

def extract_region_block_general(path, region="CAN", sheet="Electricity - generation"):
    df = pd.read_excel(path, sheet_name=sheet, header=None)

    # find all __cats__ rows
    cat_rows = df.index[df[0] == "__cats__"].tolist()

    target_cat_row = None
    for cat in cat_rows:
        region_row = cat - 3
        if region_row >= 0 and df.iloc[region_row, 0] == region:
            target_cat_row = cat
            break

    if target_cat_row is None:
        raise ValueError(f"Region {region} not found")

    header = df.iloc[target_cat_row]

    block = df.iloc[target_cat_row + 1 : target_cat_row + 16].copy()

    # ---- CRITICAL CLEANING FIX ----
    block.columns = (
        header.astype(str)
              .str.strip()
              .str.replace(r"[\x00-\x1F\x7F]", "", regex=True)
    )

    block = block.rename(columns={"__cats__": "Year"})

    # Coerce values to numeric
    block = block.apply(pd.to_numeric, errors="coerce").fillna(0)

    block = block.set_index("Year")

    return block


def aggregate_groups_all_solar(df):
    out = df.copy()
    groups = {
        "Biomass (all)": "biomass",
        "Coal (all)":    "coal",
        "Gas (all)":     "gas",
        "Oil (all)":     "oil",
        "Solar (all)":   "solar", # old one combining all solars
        "Wind (all)":    "wind",
    }
    drops=[]
    for new,prefix in groups.items():
        cols=[c for c in out.columns if isinstance(c,str) and c.lower().startswith(prefix)]
        if cols:
            out[new]=out[cols].sum(axis=1)
            drops+=cols
    out=out.drop(columns=drops)
    tiny=[c for c in out.columns if out[c].abs().max()<1e-6]
    out=out.drop(columns=tiny)
    totals=out.sum(axis=0).sort_values(ascending=False)
    return out[totals.index]


def aggregate_groups(df):
    out = df.copy()
    # ----- fossil & biomass groups -----
    groups = {
        "Biomass (all)": "biomass",
        "Coal (all)":    "coal",
        "Gas (all)":     "gas",
        "Oil (all)":     "oil",
        "Wind (all)":    "wind",
    }
    drops = []
    # simple prefix-based grouping for the above
    for new, prefix in groups.items():
        cols = [c for c in out.columns 
                if isinstance(c, str) and c.lower().startswith(prefix)]
        if cols:
            out[new] = out[cols].sum(axis=1)
            drops += cols

    # ----- SPECIAL SOLAR HANDLING -----
    # Solar CSP separate
    csp_cols = [c for c in out.columns if "Solar CSP" == c]
    # Identify Solar PV technologies (commercial&residential)
    pv_cols = [c for c in out.columns 
               if isinstance(c, str) and c.lower().startswith("solar pv")]
               #if c in ["Solar PV Centralized", "Solar PV Residential"]]

    if pv_cols:
        out["Solar PV"] = out[pv_cols].sum(axis=1) 
        drops += pv_cols
    # DO NOT drop Solar CSP, Only drop solar-PV components
    out = out.drop(columns=drops)
    # ----- remove negligible columns -----
    tiny = [c for c in out.columns if out[c].abs().max() < 1e-6]
    out = out.drop(columns=tiny)
    # ----- by magnitude (largest → smallest) -----
    totals = out.sum(axis=0).sort_values(ascending=False)
    out = out[totals.index]

    return out




def plot_region_mix(region, scenario_paths, iam_scnname = "IMAGE", mode="EJ", up_to_year=2050, palette="Set3", output=False, outdir="plot_outputs"):
    cmap = plt.cm.get_cmap(palette).colors
    fig, axes = plt.subplots(1, len(scenario_paths), figsize=(22, 8), sharey=True)
    legend = None

    # Preferred stacking order (bottom → top)
    preferred_order = [
        "Hydro",
        "Coal (all)",
        "Gas (all)",
        "Oil (all)",
        "Nuclear",
        "Storage, Battery",
        "Biomass (all)",      # << keep this high otherwise biomass is invisible like in previous version 
        "Wind (all)",
        "Solar PV",
        "Solar CSP",
        "Geothermal",
    ]

    for ax, (ssp, path) in zip(axes, scenario_paths.items()):
        df = aggregate_groups(extract_region_block_general(path, region))
        df = df[df.index <= up_to_year]

        if mode == "percent":
            df = df.div(df.sum(axis=1), axis=0) * 100

        # Reorder columns according to preferred_order, dropping any that don’t exist
        cols_ordered = [c for c in preferred_order if c in df.columns]
        df = df[cols_ordered]

        techs = list(df.columns)
        x = df.index.values
        y = df[techs].T.values

        ax.stackplot(x, y, colors=cmap[:len(techs)], alpha=0.9)
        ax.set_title(f"{region} – {ssp}", fontsize = 16)
        ax.set_xlim(df.index.min(), up_to_year)
        ax.grid(alpha=0.3)

        legend = techs

    axes[0].set_ylabel("EJ" if mode == "EJ" else "Share (%)")
    fig.legend(legend, bbox_to_anchor=(1.14, 0.5), loc="center right", title="technology", fontsize = 14)
    plt.suptitle(f"{region}: electricity mix ({mode}) – up to {up_to_year} \n (electricity generation under {iam_scnname})", fontsize=18)
    plt.tight_layout()
    
    if output:
        os.makedirs(outdir, exist_ok=True)
        
        fname = f"electricity_mix_{region}_{mode}_to_{up_to_year}.png"
        fpath = os.path.join(outdir, fname)
        plt.savefig(fpath, dpi=300, bbox_inches="tight")
        print(f"Saved figure to: {fpath}")
    
    plt.show()