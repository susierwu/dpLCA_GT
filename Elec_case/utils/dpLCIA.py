"""
for dynamic prospective LCIA of electricity foregrounds (not based on BW2 LCIA method):
- reads elec_SSPx_YEAR.xlsx inventories,
- finds matching CO2/CH4/N2O dynamic CF Excel files,
- builds dpIRF / dpAGWP / dpGWP kernels (central + ensemble),
- maps BW2 biosphere flows to gas/subtype/sign,
- computes time-resolved LCIA with uncertainty bands,
- offers plotting utilities for:
    (1) each individual BW2 GHG flow,
    (2) aggregate CO2 categories (2a–2d),
    (3) aggregate CH4 categories (3a–3d),
    (4) N2O aggregate,
    (5) all GHGs combined.

Assumptions:
- electricity Excel has columns: ["date", "amount", "flow", "activity", ...]
- GHG Excel naming embeds scenario & year in the filename 
"""

from __future__ import annotations
import os
import re
import bw2data
from bw2data import get_activity
from dataclasses import dataclass
from typing import Dict, Iterable, Literal, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# 1. Utilities: scenario/year parsing & GHG file discovery
# -----------------------------------------------------------------------------
def parse_scenario_year_from_elec(path: str) -> Tuple[str, int]:
    """
    Parse SSP scenario & model year from electricity filename.
    Output filename pattern (example):
        market_for_electricity_hydro_high_voltage_*_SSP1-VLLO_2030.xlsx
        market_for_electricity_hydro_high_voltage_*_SSP5-H_2050.xlsx
    mapping rules:
        SSP1-* -> ssp119
        SSP2-* -> ssp245
        SSP5-* -> ssp585
    Returns
    -------
    scenario : str   (e.g. "ssp119", "ssp245", "ssp585")
    year     : int   (e.g. 2030, 2050)
    """
    fname = os.path.basename(path)
    m = re.search(
        r"_(SSP[125])-[^_]+_(\d{4})\.xlsx$",
        fname,
        flags=re.IGNORECASE,
    )

    if not m:
        raise ValueError(f"Cannot parse SSP scenario/year from '{fname}'")

    ssp_family, year_str = m.groups()

    SSP_MAP = {
        "SSP1": "ssp119",
        "SSP2": "ssp245",
        "SSP5": "ssp585",
    }
    scenario = SSP_MAP[ssp_family.upper()]
    year = int(year_str)
    return scenario, year




def find_ghg_files(
    ghg_dir: str,
    scenario: str,
    year: int,
    gases: Iterable[str] = ("CO2", "CH4", "N2O"),
) -> Dict[str, str]:
    """
    Search ghg_dir for dynamic CF Excel files for each gas that match
    the given scenario & year.

    Relies only on substring matching: gas string (CO2/CH4/N2O), scenario
    (ssp119, ssp2-m, etc.) and year (e.g. '2040') must all appear in the
    filename.

    Returns
    -------
    {gas: filepath}
    """
    files = os.listdir(ghg_dir)
    mapping: Dict[str, str] = {}

    for gas in gases:
        candidates = [
            os.path.join(ghg_dir, f)
            for f in files
            if gas in f
            and scenario in f.lower()
            and str(year) in f
            and f.lower().endswith(".xlsx")
        ]
        if not candidates:
            raise FileNotFoundError(
                f"No GHG Excel file found for gas={gas}, scenario={scenario}, year={year} in {ghg_dir}"
            )
        if len(candidates) > 1:
            # Pick the first, but warn in logs
            print(
                f"[dynamic_ghg_lcia] Multiple matches for {gas}, {scenario}, {year}; using '{candidates[0]}'"
            )
        mapping[gas] = candidates[0]

    return mapping


# -----------------------------------------------------------------------------
# 2. Load kernels (dpIRF / AGWP / GWP / DCF) central + ensemble
# -----------------------------------------------------------------------------

Metric = Literal["rf", "agwp", "gwp", "dcf"]


def _find_sheet(
    xls: pd.ExcelFile,
    metric: Metric,
    gas: str,
    year: int,
    kind: Literal["central", "ensemble"],
) -> str:
    """
    Robustly find the appropriate sheet for a given metric, gas & year.

    - Ensemble sheets start with   '{metric}_wh_ensmb'
    - Central   sheets start with  '{metric}_pointvalue'

    Example names from your CO2 file:
        'rf_wh_ensmb_CO22040'
        'rf_pointvalueCO2_2040'
        'agwp_wh_ensmb_CO22040'
        'agwp_pointvalueCO2_2040'
    """
    if kind == "ensemble":
        prefix = f"{metric}_wh_ensmb"
    else:
        prefix = f"{metric}_pointvalue"

    candidates = [
        s
        for s in xls.sheet_names
        if s.startswith(prefix)
        and gas in s
        and str(year) in s  # catches "CO22040" and "_2040"
    ]
    if not candidates:
        raise ValueError(
            f"No {kind} sheet found for metric={metric}, gas={gas}, year={year} "
            f"in sheets={xls.sheet_names}"
        )
    if len(candidates) > 1:
        print(
            f"[dynamic_ghg_lcia] Multiple {kind} sheets for metric={metric}, gas={gas}, year={year}; using '{candidates[0]}'"
        )
    return candidates[0]



def load_ghg_kernels(
    ghg_files: Dict[str, str],
    metric: Metric,
    max_horizon: int = 100,   # keep only 0..100 by default
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Load central + ensemble kernels for each gas.

    Parameters
    ----------
    ghg_files : {gas: filepath}
        Map from "CO2"/"CH4"/"N2O" to dynamic CF Excel files.
    metric : {"rf", "agwp", "gwp", "dcf"}
        Which family of worksheets to use.
    max_horizon : int
        Maximum horizon to keep from Excel kernels.
        max_horizon=100 means rows for index 0..100, i.e. 101 rows.
    """
    irfs: Dict[str, Dict[str, np.ndarray]] = {}
    n_keep = max_horizon + 1

    for gas, path in ghg_files.items():
        xls = pd.ExcelFile(path)

        fname = os.path.basename(path)
        m = re.search(r"MY(\d{4})", fname)
        if not m:
            raise ValueError(
                f"Cannot infer model year from '{fname}'. Expected pattern '*MYXXXX*'"
            )
        year = int(m.group(1))

        # -------------------
        # Central
        # -------------------
        central_sheet = _find_sheet(
            xls, metric=metric, gas=gas, year=year, kind="central"
        )
        df_central = pd.read_excel(xls, sheet_name=central_sheet)

        gas_col = [c for c in df_central.columns if gas in str(c) or str(c) == gas]
        if not gas_col:
            raise ValueError(f"No '{gas}' column in central sheet '{central_sheet}'")

        central = df_central[gas_col[0]].to_numpy()

        # -------------------
        # Ensemble
        # -------------------
        ens_sheet = _find_sheet(
            xls, metric=metric, gas=gas, year=year, kind="ensemble"
        )
        df_ens = pd.read_excel(xls, sheet_name=ens_sheet)
        ens = df_ens.drop(columns=[df_ens.columns[0]]).to_numpy()

        # -------------------
        # Harmonize horizon
        # -------------------
        if len(central) < n_keep:
            raise ValueError(
                f"{gas} central kernel in '{os.path.basename(path)}' has only "
                f"{len(central)} rows, but needs at least {n_keep} rows "
                f"for horizon 0..{max_horizon}."
            )
        if ens.shape[0] < n_keep:
            raise ValueError(
                f"{gas} ensemble kernel in '{os.path.basename(path)}' has only "
                f"{ens.shape[0]} rows, but needs at least {n_keep} rows "
                f"for horizon 0..{max_horizon}."
            )

        if len(central) > n_keep or ens.shape[0] > n_keep:
            print(
                f"[dynamic_ghg_lcia] Truncating {gas} kernels from "
                f"{len(central)} / {ens.shape[0]} rows to {n_keep} rows "
                f"(index 0..{max_horizon}) from file '{os.path.basename(path)}'"
            )

        central = central[:n_keep]
        ens = ens[:n_keep, :]

        irfs[gas] = {
            "central": central,
            "ensemble": ens,
        }

    return irfs
    

# -----------------------------------------------------------------------------
# 3. BW2 flow metadata: gas/subtype/sign classification
# -----------------------------------------------------------------------------

@dataclass
class FlowInfo:
    flow_id: int
    name: str
    gas: Optional[str]         # "CO2", "CH4", "N2O" or None for non-GHG
    subtype: Optional[str]     # aggregation key, see below
    sign: int                  # +1 or −1 applied to kernel


def classify_flow_from_name(name: str) -> FlowInfo:
    """
    classifier for CO₂ / CH₄ / N₂O that follows ghg category rules:
    CO2:
      2.a fossil positive
      2.b biogenic (pos/neg except resource correction/in air)
      2.c negative other: 
            - Carbon dioxide, in air*
            - Carbon dioxide, non-fossil, resource correction
      2.d total (handled at aggregation step)
    CH4: fossil / non-fossil / biomass
    N2O: single category
    All matching is case-insensitive.
    """
    n = name.strip().lower()

    # ------------------------------------------------------------------
    # CO₂ family
    # ------------------------------------------------------------------
    if "carbon dioxide" in n:

        # ---------- Category 2.c (NEGATIVE OTHER: air + resource correction) ----------
        # Matches:
        #   "carbon dioxide, in air"
        #   "carbon dioxide, in air, long-term"
        #   "carbon dioxide, non-fossil, resource correction"
        if n.startswith("carbon dioxide, in air"):
            return FlowInfo(
                flow_id=-1,
                name=name,
                gas="CO2",
                subtype="co2_other_negative",   # category 2.c
                sign=-1,
            )

        if n == "carbon dioxide, non-fossil, resource correction":
            return FlowInfo(
                flow_id=-1,
                name=name,
                gas="CO2",
                subtype="co2_other_negative",   # category 2.c
                sign=-1,
            )

        # ---------- Category 2.b (BIOGENIC CO2) ----------
        # Includes:
        #   "carbon dioxide, non-fossil"
        #   "carbon dioxide, from soil or biomass stock" (positive)
        #   "carbon dioxide, to soil or biomass stock" (negative)
        if "non-fossil" in n:
            # skip resource correction case (already handled above)
            if "resource correction" not in n:
                return FlowInfo(
                    flow_id=-1,
                    name=name,
                    gas="CO2",
                    subtype="co2_biogenic",   # category 2.b
                    sign=+1,
                )

        if "from soil or biomass stock" in n:
            return FlowInfo(
                flow_id=-1,
                name=name,
                gas="CO2",
                subtype="co2_biogenic",
                sign=+1,
            )

        if "to soil or biomass stock" in n:
            return FlowInfo(
                flow_id=-1,
                name=name,
                gas="CO2",
                subtype="co2_biogenic",
                sign=-1,
            )

        # ---------- Category 2.a (FOSSIL + POSITIVE CO2) ----------
        if "fossil" in n and "non-fossil" not in n:
            return FlowInfo(
                flow_id=-1,
                name=name,
                gas="CO2",
                subtype="co2_fossil_positive",
                sign=+1,
            )

        # Fallback: treat any remaining CO2 as fossil-positive
        return FlowInfo(
            flow_id=-1,
            name=name,
            gas="CO2",
            subtype="co2_fossil_positive",
            sign=+1,
        )

    # ------------------------------------------------------------------
    # CH₄ family
    # ------------------------------------------------------------------
    if "methane" in n:
        gas = "CH4"
        if "fossil" in n and "non-fossil" not in n:
            sub = "ch4_fossil"
        elif "non-fossil" in n:
            sub = "ch4_non_fossil"
        elif "from soil or biomass stock" in n:
            sub = "ch4_biomass"
        else:
            sub = "ch4_other"
        return FlowInfo(-1, name, gas, sub, +1)

    # ------------------------------------------------------------------
    # N₂O
    # ------------------------------------------------------------------
    if "dinitrogen monoxide" in n or "nitrous oxide" in n:
        return FlowInfo(-1, name, "N2O", "n2o_all", +1)

    # Non-GHG
    return FlowInfo(-1, name, None, None, +1)




def build_flow_info_mapping(flow_ids: Iterable[int]) -> Dict[int, FlowInfo]:
    """
    Use BW2 to look up flow names and classify them.
    Must call this in a Brightway2 context.
    Parameters
    ----------
    flow_ids : iterable of biosphere flow IDs (ints)
    Returns
    -------
    {flow_id: FlowInfo}
    """
    try:
        from bw2data import get_activity
    except ImportError:
        raise ImportError("bw2data is required to build flow information mapping")

    mapping: Dict[int, FlowInfo] = {}
    for fid in flow_ids:
        act = get_activity(int(fid))
        name = act["name"]
        tmpl = classify_flow_from_name(name)
        info = FlowInfo(
            flow_id=int(fid),
            name=name,
            gas=tmpl.gas,
            subtype=tmpl.subtype,
            sign=tmpl.sign,
        )
        mapping[int(fid)] = info

    return mapping


# -----------------------------------------------------------------------------
# 4. dp-LCIA computation (central + uncertainty) by category, adding H parameter
# -----------------------------------------------------------------------------
def compute_dynamic_lcia_by_category(
    elec_path: str,
    ghg_dir: str,
    metric: Metric = "rf",
    max_horizon: int = 100,
    ): # -> Tuple[pd.DatetimeIndex, pd.DataFrame, pd.DataFrame, Dict[int, FlowInfo]]:

    # ---- 1. Read electricity inventory ----
    # Ensure read "flow" as str (converted before saving to excel) 
    elec = pd.read_excel(elec_path, dtype={"flow": str})
    # Ensure dates are parsed
    elec["date"] = pd.to_datetime(elec["date"], errors="coerce")
    # Extract year (needed for pivot)
    elec["year"] = elec["date"].dt.year
    # convert back str(flow) to int64 so that bw2data.get_activity() works
    elec["flow"] = elec["flow"].astype("int64")

    # New pivot table (year × flow)
    pivot = elec.pivot_table(
        index="year",
        columns="flow",
        values="amount",
        aggfunc="sum"
    ).fillna(0)

    flow_ids = pivot.columns.astype(int).tolist()

    # ---- 2. Scenario/year + GHG kernels ----
    scenario, year = parse_scenario_year_from_elec(elec_path)
    ghg_files = find_ghg_files(ghg_dir=ghg_dir, scenario=scenario, year=year)
    irfs = load_ghg_kernels(ghg_files, metric=metric, max_horizon=max_horizon)

    # ---- 3. Flow metadata ----
    flow_info = build_flow_info_mapping(flow_ids)

    # ---- 4. Evaluation horizon ----
    first_year = int(pivot.index[0])
    last_year  = int(pivot.index[-1])

    #eval_years = pd.date_range(
    #    f"{first_year}-01-01",
    #    f"{last_year + 100}-01-01",
    #    freq="YS",
    #)
    # now with updated H parameter, do not hard-code H = 100 
    eval_years = pd.date_range(
        f"{first_year}-01-01",
        f"{last_year + max_horizon}-01-01",
        freq="YS",
    )
    
    n_eval = len(eval_years)

    # ---- 5. Dynamic convolution per subtype ----
    subtypes = sorted({
        fi.subtype
        for fi in flow_info.values()
        if fi.gas in irfs and fi.subtype is not None
    })

    central = pd.DataFrame(0.0, index=eval_years, columns=subtypes)
    band_lower = pd.DataFrame(0.0, index=eval_years, columns=subtypes)
    band_upper = pd.DataFrame(0.0, index=eval_years, columns=subtypes)

    # length and ensemble count
    ens_example = next(iter(irfs.values()))["ensemble"]
    n_steps, n_ens = ens_example.shape

    for subtype in subtypes:
        cols = [fid for fid, fi in flow_info.items()
                if fi.subtype == subtype and fi.gas in irfs]

        scores_central = np.zeros(n_eval)
        scores_ens = np.zeros((n_eval, n_ens))

        for fid in cols:
            fi = flow_info[fid]
            kernels = irfs[fi.gas]

            kernel_c = fi.sign * kernels["central"]
            kernel_e = fi.sign * kernels["ensemble"]

            series = pivot[fid]

            for date, amt in series.items():
                if amt == 0:
                    continue

                e_year = int(date)                 # correct for pivot
                start_idx = e_year - first_year

                scores_central[start_idx:start_idx+n_steps] += amt * kernel_c
                scores_ens[start_idx:start_idx+n_steps, :] += amt * kernel_e

        central[subtype] = scores_central
        band_lower[subtype] = np.percentile(scores_ens, 5, axis=1)
        band_upper[subtype] = np.percentile(scores_ens, 95, axis=1)

    # ---- 6. Aggregate categories (robust version) ----

    # CO2
    central["co2_total"] = (
        central.get("co2_fossil_positive", 0) +
        central.get("co2_biogenic", 0) +
        central.get("co2_other_negative", 0)
    )
    band_lower["co2_total"] = (
        band_lower.get("co2_fossil_positive", 0) +
        band_lower.get("co2_biogenic", 0) +
        band_lower.get("co2_other_negative", 0)
    )
    band_upper["co2_total"] = (
        band_upper.get("co2_fossil_positive", 0) +
        band_upper.get("co2_biogenic", 0) +
        band_upper.get("co2_other_negative", 0)
    )

    # CH4
    central["ch4_total"] = (
        central.get("ch4_fossil", 0) +
        central.get("ch4_non_fossil", 0) +
        central.get("ch4_biomass", 0)
    )
    band_lower["ch4_total"] = (
        band_lower.get("ch4_fossil", 0) +
        band_lower.get("ch4_non_fossil", 0) +
        band_lower.get("ch4_biomass", 0)
    )
    band_upper["ch4_total"] = (
        band_upper.get("ch4_fossil", 0) +
        band_upper.get("ch4_non_fossil", 0) +
        band_upper.get("ch4_biomass", 0)
    )

    # N2O
    central["n2o_total"] = central.get("n2o_all", 0)
    band_lower["n2o_total"] = band_lower.get("n2o_all", 0)
    band_upper["n2o_total"] = band_upper.get("n2o_all", 0)

    # All GHGs
    central["all_ghg"] = (
        central["co2_total"] +
        central["ch4_total"] +
        central["n2o_total"]
    )
    band_lower["all_ghg"] = (
        band_lower["co2_total"] +
        band_lower["ch4_total"] +
        band_lower["n2o_total"]
    )
    band_upper["all_ghg"] = (
        band_upper["co2_total"] +
        band_upper["ch4_total"] +
        band_upper["n2o_total"]
    )

    # ---- 7. Pack uncertainty bands ----
    band = pd.concat({"lower": band_lower, "upper": band_upper}, axis=1)

    #return eval_years, central, band, flow_info
    
    return {
        "eval_years": eval_years,
        "central": central,
        "band": band,
        "flow_info": flow_info,
        "SSP": scenario,
        "ModelYear": year 
    }


# -----------------------------------------------------------------------------
# 5. Plotting helpers
# -----------------------------------------------------------------------------
def plot_with_band(
    years: pd.DatetimeIndex,
    central: pd.Series,
    lower: pd.Series,
    upper: pd.Series,
    title: str,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    ax.fill_between(years, lower, upper, alpha=0.3)
    ax.plot(years, central, marker="o")
    ax.set_title(title)
    ax.set_ylabel("Dynamic impact (a.u.)")
    ax.grid(True)
    return ax


def plot_gas_categories(
    eval_years: pd.DatetimeIndex,
    central: pd.DataFrame,
    band: pd.DataFrame,
    gas: Literal["CO2", "CH4", "N2O", "ALL"],
):
    """
    Convenience plotting wrapper implementing your category scheme:

    CO2:
        2.a co2_fossil_positive
        2.b co2_biogenic
        2.c co2_other_negative
        2.d co2_total
    CH4:
        3.a ch4_fossil
        3.b ch4_non_fossil
        3.c ch4_biomass
        3.d ch4_total
    N2O:
        n2o_total
    ALL:
        all_ghg
    """
    if gas == "CO2":
        cats = [
            ("co2_fossil_positive", "CO₂ fossil & positive"),
            ("co2_biogenic", "CO₂ biogenic (±)"),
            ("co2_other_negative", "CO₂ negative tech (-)"),
            ("co2_total", "CO₂ total"),
        ]
    elif gas == "CH4":
        cats = [
            ("ch4_fossil", "CH₄ fossil"),
            ("ch4_non_fossil", "CH₄ non-fossil"),
            ("ch4_biomass", "CH₄ biomass"),
            ("ch4_total", "CH₄ total"),
        ]
    elif gas == "N2O":
        cats = [("n2o_total", "N₂O total")]
    elif gas == "ALL":
        cats = [("all_ghg", "All major GHGs")]
    else:
        raise ValueError(gas)

    n = len(cats)
    fig, axes = plt.subplots(n, 1, figsize=(10, 4 * n), sharex=True)
    if n == 1:
        axes = [axes]

    for ax, (col, title) in zip(axes, cats):
        plot_with_band(
            eval_years,
            central[col],
            band["lower"][col],
            band["upper"][col],
            title,
            ax=ax,
        )

    axes[-1].set_xlabel("Year")
    fig.tight_layout()
    return fig, axes




def compute_dynamic_lcia_for_scenarios(
    elec_files: Dict[Tuple[str, int], str],
    ghg_dir: str,
    metric: Metric = "rf",
    max_horizon: int = 100
):
    """
    Compute dynamic LCIA for multiple (SSP, Year) electricity files.
    Parameters
    ----------
    elec_files : dict
        Keys: (ssp, year)
        Values: path to electricity Excel for that SSP-year.
    
    Returns
    -------
    results : dict : keys of the dictionary are tuples (ssp, year), value is a dict, example below
         {
            ("SSP1-19", 2030): {
                "eval_years": DatetimeIndex([...]),
                "central": DataFrame(...),
                "band": DataFrame(...),
                "flow_info": {25040179: FlowInfo(...), ...},
            },
            ("SSP5-85", 2030): {
                ...
            },
        }
    """
    results = {}

    for (ssp, year), path in elec_files.items():
        out = compute_dynamic_lcia_by_category(
            elec_path=path,
            ghg_dir=ghg_dir,
            metric=metric,
            max_horizon=max_horizon,
        )
        results[(ssp, year)] = out

    return results





def _ssp_to_color(ssp_name: str):
    """
    Robustly map any SSP name (SSP1, SSP1-19, ssp119, SSP1_VLLO, ...)
    to the correct scenario family color.

    Returns
    -------
    hex color string
    """
    s = ssp_name.upper().replace(" ", "").replace("_", "").replace("-", "")
    # SSP1 family → SSP1-19 deep blue
    if s.startswith("SSP1"):
        return "#006BA4"
    # SSP2 family → SSP2-45 orange
    if s.startswith("SSP2"):
        return "#FF800E"
    # SSP5 family → SSP5-85 gray
    if s.startswith("SSP5"):
        return "gray"
    # fallback — black (unknown SSP)
    return "black"





def plot_multi_scenario_bands(
    results: dict,
    metric: str = "rf", 
    gas_to_plot: str = "CO2",
    subcat: str = "co2_total",
    title: str | None = None,
    figsize=(12, 6),
    #save_path: str | None = None,
    #dpi: int = 300,
):
    """
    Plot dpLCIA uncertainty bands (central + 5–95%) across multiple SSP/Year scenarios for a single GHG category.
    Parameters
    ----------
    results : dict
        Output of compute_dynamic_lcia_for_scenarios().
        Keys: (ssp_name, model_year)
        Values: {
            "eval_years": DatetimeIndex,
            "central": DataFrame,
            "band": DataFrame,
            "flow_info": dict,
        }
    gas_to_plot : {"CO2", "CH4", "N2O", "ALL"}
        Gas category to extract.
    subcat : str
        Subcategory column (e.g., "co2_total").
        If gas_to_plot == "ALL", the plot title will instead show:
        "all CO2, CH4, N2O" and subcat is ignored.
    title : str or None
        Custom title; otherwise auto-generated.
    figsize : tuple
        Figure size.
    """
       
    _METRIC_INFO = {
        "rf": (
            "instantaneous radiative forcing",
            "W m⁻² kg⁻¹",
        ),
        "agwp": (
            "AGWP (cumulative radiative forcing)",
            "W m⁻² yr kg⁻¹",
        ),
        "gwp": (
            "GWP100-equivalent metric",
            "kg CO₂-eq",
        ),
        "dcf": (
            "dynamic characterization factor",
            "a.u.",
        ),
    }
    
    # ------------------------------------------
    # Print available options as reminder
    # ------------------------------------------
    print("Gas options:  'CO2', 'CH4', 'N2O', 'ALL'")
    print("metric options: 'rf', 'agwp'")
    print("Subcategory options only for a single GHG (not for ALL): 'co2_total', 'ch4_total', 'n2o_total', \n co2_other_negative, co2_biogenic, co2_fossil_positive, ch4_fossil, ch4_non_fossil, ch4_biomass")

    # -------------------------------
    # Metric info lookup
    # -------------------------------
    metric_key = metric.lower()
    if metric_key not in _METRIC_INFO:
        raise ValueError(f"Unknown metric '{metric}' — must be one of {list(_METRIC_INFO.keys())}")
    metric_fullname, metric_unit = _METRIC_INFO[metric_key]

    # ------------------------------------------
    # Build display name for gas/subcat
    # ------------------------------------------
    if gas_to_plot.upper() == "ALL":
        disp_cat = "all CO₂, CH₄, N₂O"
        subcat_used = "all_ghg"   # override
    else:
        disp_cat = f"{gas_to_plot.upper()} – {subcat}"
        subcat_used = subcat

    # ------------------------------------------
    # Build default title
    # ------------------------------------------
    if title is None:
        title = f"dpLCIA ({metric_fullname}) multi-scenario – {disp_cat}"

    # -----------------------------------------
    # Plot
    # ------------------------------------------
    plt.figure(figsize=figsize)

    for (ssp, model_year), out in results.items():

        eval_years = out["eval_years"]
        central_df = out["central"]
        band_df = out["band"]

        if subcat_used not in central_df.columns:
            print(f"WARNING: '{subcat_used}' not found in scenario {ssp} {model_year}")
            continue

        central = central_df[subcat_used]
        lower = band_df["lower"][subcat_used]
        upper = band_df["upper"][subcat_used]

        color = _ssp_to_color(ssp)
        label = f"{ssp} – {model_year}"

        # band
        plt.fill_between(eval_years, lower, upper, color=color, alpha=0.18)

        # central curve
        plt.plot(eval_years, central, lw=2.2, color=color, label=label)

    plt.grid(True)
    plt.xlabel("Year", fontsize =12 )
    plt.ylabel(f"{metric} [{metric_unit}]", fontsize =14)
    plt.title(title, fontsize =16)
    plt.legend(fontsize =16)
    plt.tight_layout()
    plt.show()





def plot_multi_scenario_split_each_ghg_totalscore(
    results: dict,
    metric: str = "rf",
    gases_to_plot: list = None,     #  only totals
    ssp_color_map: dict = None,
    plot_uncertainty: bool = True,  
    figsize: tuple = (14, 8),
    x_limit: int | None = None,
    x_end_year: int | None = None,
    left_pad_years: int = 5,  # with the new x_limit, to have some margin on the left side 
    save_path: str | None = None,
    dpi: int = 300,
    
):
    """
    One combined dpLCIA plot comparing multiple SSP scenarios.
    DEFAULT behavior:
        Plot ONLY the four major totals:
            all_ghg, co2_total, ch4_total, n2o_total

    If gases_to_plot is provided:
        Plot exactly those categories (totals OR subcategories).
    """

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    # ===========================================================
    # 1. Collect available columns (from first results entry)
    # ===========================================================
    sample_out = next(iter(results.values()))
    available_cols = list(sample_out["central"].columns)

    # ===========================================================
    # 2. Define label and style registry (for totals + subcats)
    # ===========================================================
    style_registry = {
        # Totals
        "all_ghg": ("All GHGs", "-"),
        "co2_total": ("CO₂ total", "-."), #--
        "ch4_total": ("CH₄ total", ":"), #-.
        "n2o_total": ("N₂O total", "--"),

        # Subcategories
        "co2_fossil_positive": ("CO₂ fossil", "--"),
        "co2_biogenic": ("CO₂ biogenic", ":"),
        "co2_other_negative": ("CO₂ negative tech", "-."),

        "ch4_fossil": ("CH₄ fossil", "--"),
        "ch4_non_fossil": ("CH₄ non-fossil", ":"),
        "ch4_biomass": ("CH₄ biomass", "-."),
    }

    # ===========================================================
    # 3. DEFAULT GASES TO PLOT (major totals ONLY)
    # ===========================================================
    default_totals = ["all_ghg", "co2_total", "ch4_total", "n2o_total"]

    if gases_to_plot is None:
        gases_to_plot = [g for g in default_totals if g in available_cols]

    # Validate chosen gases exist
    for g in gases_to_plot:
        if g not in available_cols:
            raise ValueError(
                f"Requested gas '{g}' not found.\n"
                f"Available: {available_cols}"
            )

    # Build ordered list
    gas_order = [(g, style_registry[g][0], style_registry[g][1])
                 for g in gases_to_plot]

    # ===========================================================
    # 4. Metric metadata
    # ===========================================================
    metric_map = {
        "rf": ("instantaneous radiative forcing", "W m⁻² kg⁻¹"),
        "agwp": ("AGWP (cumulative radiative forcing)", "W m⁻² yr kg⁻¹"),
        "gwp": ("GWP100-equivalent weighting", "kg CO₂-eq / kg"),
        "dcf": ("dynamic characterization factor", "a·kg⁻¹"),
    }
    metric_fullname, metric_unit = metric_map.get(metric.lower(), ("metric", ""))

    # ===========================================================
    # 5. SSP color map
    # ===========================================================
    if ssp_color_map is None:
        ssp_color_map = {
            "SSP1-19": "#006BA4",
            "SSP2-45": "#FF800E",
            "SSP5-85": "gray",
        }

    plt.figure(figsize=figsize)

    # ===========================================================
    # 6a. Plot curves and (optional) uncertainty bands
    # ===========================================================
    for gas_key, gas_label, style in gas_order:
        for (ssp, year), out in results.items():
    
            eval_years = out["eval_years"]
            central = out["central"][gas_key]
    
            color = ssp_color_map.get(ssp, "black")
    
            # --- uncertainty band (optional) ---
            if plot_uncertainty:
                lower = out["band"]["lower"][gas_key]
                upper = out["band"]["upper"][gas_key]
    
                plt.fill_between(
                    eval_years,
                    lower,
                    upper,
                    color=color,
                    alpha=0.12,
                    zorder=1,
                )
    
            # --- central trajectory ---
            plt.plot(
                eval_years,
                central,
                linestyle=style,
                linewidth=2,
                color=color,
                zorder=2,
            )

    # ===========================================================
    # 6b. CO2 > CH4 cumulative AGWP crossover (per scenario) 
    # ADD a  horizontal line (axhline) in the same SSP color, where co2_total > ch4_total
    # ===========================================================
    if metric == "agwp": 
        ax = plt.gca()
        for (ssp, model_year), out in results.items():
            central = out["central"]
            if not {"co2_total", "ch4_total"}.issubset(central.columns):
                continue
    
            diff = central["co2_total"] - central["ch4_total"]
            crossover_mask = diff > 0
            # no crossover within horizon
            if not crossover_mask.any():
                continue

            # first year where CO2 surpasses CH4
            t_cross = crossover_mask.idxmax()   # Timestamp
            cross_year = int(pd.Timestamp(t_cross).year)
            color = ssp_color_map.get(ssp, "black")
            # --- vertical line at crossover year ---
            ax.axvline(
                x=t_cross,
                color=color,
                linestyle=":",
                linewidth=2.0,
                alpha=0.9,
                zorder=0,
            )
            # --- annotate the year near the top of the plot ---
            # place text using axis fraction for y so it stays visible regardless of scale
            ax.text(
                t_cross,
                0.98,  # near top
                f"{ssp}: CO₂ > CH₄ ({cross_year})",
                color=color,
                fontsize=10,
                rotation=90,
                va="top",
                ha="right",
                transform=ax.get_xaxis_transform(),  # x in data, y in axes fraction
            )
    
    # ===========================================================
    # 7. Legend 1 — SSP scenarios (top right)
    # ===========================================================
    ssp_patches = [
        Patch(
            facecolor=ssp_color_map.get(ssp, "black"),
            edgecolor="none",
            alpha=0.6,  # clearer than band
            label=f"{ssp} – {year}"
        )
        for (ssp, year) in results.keys()
    ]

    legend1 = plt.legend(
        handles=ssp_patches,
        title="Scenarios",
        loc="upper right",
        bbox_to_anchor=(1.0, 1.00),
        fontsize=11,
        title_fontsize=11,
        frameon=True,
        borderpad=0.3,
    )
    plt.gca().add_artist(legend1)

    # ===========================================================
    # 8. Legend 2 — Gas categories (directly under legend1)
    # ===========================================================
    gas_lines = [
        Line2D([0], [0], color="black", linestyle=sty, linewidth=2, label=lbl)
        for (_, lbl, sty) in gas_order
    ]

    plt.legend(
        handles=gas_lines,
        title="GHG",
        loc="upper right",
        bbox_to_anchor=(1.0, 0.84),  # slightly below legend1
        fontsize=11,
        title_fontsize=11,
        frameon=True,
        borderpad=0.3,
        handlelength=2.0,
    )

    # ===========================================================
    # 9. Final axes, titles
    # ===========================================================
    plt.grid(True)
    plt.xlabel("Year", fontsize=12)
    plt.ylabel(f"{metric} [{metric_unit}]", fontsize=14)

    #desc = ", ".join(gases_to_plot)
    plt.title(
        f"dpLCIA ({metric_fullname}) ", #— {desc}
        fontsize=16,
    )


    # new with the x-limit so no sudden drop of the impacts along x-axis    # apply x-range cutoff here
    ax = plt.gca()
    
    sample_out = next(iter(results.values()))
    first_plot_year = pd.Timestamp(sample_out["eval_years"][0]).year
    
    if x_limit is not None and x_end_year is not None:
        raise ValueError("Use either x_limit or x_end_year, not both.")
    
    if x_limit is not None:
        ax.set_xlim(
            pd.Timestamp(f"{first_plot_year - left_pad_years}-01-01"),
            pd.Timestamp(f"{first_plot_year + x_limit}-01-01"),
        )
    
    if x_end_year is not None:
        ax.set_xlim(
            pd.Timestamp(f"{first_plot_year - left_pad_years}-01-01"),
            pd.Timestamp(f"{x_end_year}-01-01"),
        )
    
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")

    plt.show()