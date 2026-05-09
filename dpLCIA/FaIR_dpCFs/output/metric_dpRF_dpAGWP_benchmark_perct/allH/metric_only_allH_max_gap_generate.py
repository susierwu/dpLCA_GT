from pathlib import Path
import pickle

import nbformat as nbf
import numpy as np
import pandas as pd


BASE = Path("/Users/susierwu/dpLCA_main/dpLCIA/FaIR_dpCFs/output/metrics_AGWP_wh_ModC_CC_T")
OUTDIR = Path("/Users/susierwu/Documents/Codex/2026-05-08/files-mentioned-by-the-user-results")

SCENARIOS = [("SSP1-1.9", "ssp119"), ("SSP2-4.5", "ssp245"), ("SSP5-8.5", "ssp585")]
MODEL_YEARS = [2030, 2040, 2050]
GASES = ["CO2", "CH4", "N2O"]
METRICS = ["dpRF", "dpAGWP"]

METRIC_SPECS_FINAL = {
    ("CO2", "dpRF"): "rf",
    ("CO2", "dpAGWP"): "agwp",
    ("CH4", "dpRF"): "rf_final",
    ("CH4", "dpAGWP"): "agwp_final",
    ("N2O", "dpRF"): "rf_final",
    ("N2O", "dpAGWP"): "agwp_final",
}


def flow_label(gas):
    return {"CO2": "CO2 total", "CH4": "CH4 total", "N2O": "N2O total"}[gas]


def load_obj(scenario_code, model_year):
    path = BASE / f"metrics_by_pIRF_tstep1ALL_{scenario_code}_fair_start1750MY{model_year}.pkl"
    with path.open("rb") as f:
        return pickle.load(f)


def horizon_series(obj, gas, metric_key):
    arr = np.asarray(obj["metrics"][gas][metric_key], dtype=float)
    values = np.nanmean(arr, axis=1) if arr.ndim == 2 else arr
    horizons = np.asarray(obj["meta"]["H"], dtype=float)
    return horizons, values


def build_all_h_table():
    rows = []
    for model_year in MODEL_YEARS:
        cache = {(label, code): load_obj(code, model_year) for label, code in SCENARIOS}
        for metric in METRICS:
            for gas in GASES:
                metric_key = METRIC_SPECS_FINAL[(gas, metric)]
                series = {}
                horizons = None
                for label, code in SCENARIOS:
                    horizons, values = horizon_series(cache[(label, code)], gas, metric_key)
                    series[label] = values
                for i, horizon in enumerate(horizons):
                    s1 = float(series["SSP1-1.9"][i])
                    s2 = float(series["SSP2-4.5"][i])
                    s5 = float(series["SSP5-8.5"][i])
                    rows.append(
                        {
                            "Metric": metric,
                            "Metric source": metric_key,
                            "ModelYear": model_year,
                            "Horizon H": int(horizon),
                            "Flow category": flow_label(gas),
                            "SSP1-1.9 value": s1,
                            "SSP2-4.5 value": s2,
                            "SSP5-8.5 value": s5,
                            "SSP2-4.5 - SSP1-1.9": s2 - s1,
                            "SSP5-8.5 - SSP1-1.9": s5 - s1,
                            "Difference SSP2-4.5 vs SSP1-1.9": (s2 - s1) / s1 * 100 if s1 != 0 else np.nan,
                            "Difference SSP5-8.5 vs SSP1-1.9": (s5 - s1) / s1 * 100 if s1 != 0 else np.nan,
                        }
                    )
    return pd.DataFrame(rows)


def max_gap_table(df, gap_type):
    rows = []
    comparisons = [
        ("SSP2-4.5 vs SSP1-1.9", "SSP2-4.5 - SSP1-1.9", "Difference SSP2-4.5 vs SSP1-1.9"),
        ("SSP5-8.5 vs SSP1-1.9", "SSP5-8.5 - SSP1-1.9", "Difference SSP5-8.5 vs SSP1-1.9"),
    ]
    group_cols = ["Metric", "ModelYear", "Flow category", "Metric source"]
    for group_vals, sub in df.groupby(group_cols, dropna=False):
        for comp_label, diff_col, pct_col in comparisons:
            if gap_type == "absolute_difference":
                idx = sub[diff_col].abs().idxmax()
            elif gap_type == "absolute_percent":
                idx = sub[pct_col].replace([np.inf, -np.inf], np.nan).abs().idxmax()
            else:
                raise ValueError(gap_type)
            r = df.loc[idx]
            rows.append(
                {
                    "Gap type": gap_type,
                    "Comparison": comp_label,
                    "Metric": r["Metric"],
                    "ModelYear": int(r["ModelYear"]),
                    "Flow category": r["Flow category"],
                    "Metric source": r["Metric source"],
                    "Horizon of maximum gap": int(r["Horizon H"]),
                    "SSP1-1.9 value": r["SSP1-1.9 value"],
                    "SSP2-4.5 value": r["SSP2-4.5 value"],
                    "SSP5-8.5 value": r["SSP5-8.5 value"],
                    "Absolute difference at H": r[diff_col],
                    "Percent difference at H": r[pct_col],
                    "Abs SSP1 denominator at H": abs(r["SSP1-1.9 value"]),
                }
            )
    return pd.DataFrame(rows)


def global_max_table(df, gap_type):
    rows = []
    comparisons = [
        ("SSP2-4.5 vs SSP1-1.9", "SSP2-4.5 - SSP1-1.9", "Difference SSP2-4.5 vs SSP1-1.9"),
        ("SSP5-8.5 vs SSP1-1.9", "SSP5-8.5 - SSP1-1.9", "Difference SSP5-8.5 vs SSP1-1.9"),
    ]
    for (metric, flow), sub in df.groupby(["Metric", "Flow category"], dropna=False):
        for comp_label, diff_col, pct_col in comparisons:
            if gap_type == "absolute_difference":
                idx = sub[diff_col].abs().idxmax()
            else:
                idx = sub[pct_col].replace([np.inf, -np.inf], np.nan).abs().idxmax()
            r = df.loc[idx]
            rows.append(
                {
                    "Gap type": gap_type,
                    "Comparison": comp_label,
                    "Metric": r["Metric"],
                    "Flow category": r["Flow category"],
                    "ModelYear": int(r["ModelYear"]),
                    "Metric source": r["Metric source"],
                    "Horizon of maximum gap": int(r["Horizon H"]),
                    "SSP1-1.9 value": r["SSP1-1.9 value"],
                    "SSP2-4.5 value": r["SSP2-4.5 value"],
                    "SSP5-8.5 value": r["SSP5-8.5 value"],
                    "Absolute difference at H": r[diff_col],
                    "Percent difference at H": r[pct_col],
                    "Abs SSP1 denominator at H": abs(r["SSP1-1.9 value"]),
                }
            )
    return pd.DataFrame(rows)


def format_for_reading(df):
    out = df.copy()
    sci_cols = [
        "SSP1-1.9 value",
        "SSP2-4.5 value",
        "SSP5-8.5 value",
        "Absolute difference at H",
        "Abs SSP1 denominator at H",
    ]
    for col in sci_cols:
        if col in out.columns:
            out[col] = out[col].map(lambda x: f"{x:.3e}")
    if "Percent difference at H" in out.columns:
        out["Percent difference at H"] = out["Percent difference at H"].map(lambda x: f"{x:+.1f}%")
    return out


def write_notebook(path):
    code = Path(__file__).read_text()
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Metric-only all-horizon maximum gap benchmark\n\n"
            "This notebook loops through every metric horizon in the Module C pickle outputs "
            "(`H = 0...200`) and identifies the H where the SSP2-4.5 and SSP5-8.5 gaps "
            "relative to SSP1-1.9 are largest. It reports both maximum absolute percent "
            "gap and maximum absolute metric-unit gap, because dpRF can cross or approach zero."
        ),
        nbf.v4.new_code_cell(code),
    ]
    nbf.write(nb, path)


def main():
    all_h = build_all_h_table()
    max_pct = max_gap_table(all_h, "absolute_percent")
    max_abs = max_gap_table(all_h, "absolute_difference")
    global_pct = global_max_table(all_h, "absolute_percent")
    global_abs = global_max_table(all_h, "absolute_difference")

    xlsx_path = OUTDIR / "metric_only_allH_max_gap_ssp119_ssp245_ssp585.xlsx"
    csv_all = OUTDIR / "metric_only_allH_values_ssp119_ssp245_ssp585.csv"
    csv_global_pct = OUTDIR / "metric_only_allH_global_max_percent_gap.csv"
    csv_global_abs = OUTDIR / "metric_only_allH_global_max_absolute_gap.csv"
    nb_path = OUTDIR / "metric_only_allH_max_gap_ssp119_ssp245_ssp585.ipynb"

    all_h.to_csv(csv_all, index=False)
    global_pct.to_csv(csv_global_pct, index=False)
    global_abs.to_csv(csv_global_abs, index=False)
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        format_for_reading(global_pct).to_excel(writer, sheet_name="global_max_percent", index=False)
        format_for_reading(global_abs).to_excel(writer, sheet_name="global_max_absdiff", index=False)
        format_for_reading(max_pct).to_excel(writer, sheet_name="by_MY_max_percent", index=False)
        format_for_reading(max_abs).to_excel(writer, sheet_name="by_MY_max_absdiff", index=False)
        all_h.to_excel(writer, sheet_name="all_H_numeric", index=False)
        pd.DataFrame(
            {
                "Item": ["Input folder", "Horizon range", "Metric convention", "Percent formula"],
                "Value": [
                    str(BASE),
                    "All H values available in pickles: 0...200",
                    "CO2 uses rf/agwp; CH4 and N2O use rf_final/agwp_final",
                    "(SSPx - SSP1-1.9) / SSP1-1.9 * 100",
                ],
            }
        ).to_excel(writer, sheet_name="notes", index=False)
    write_notebook(nb_path)

    print(f"WROTE {xlsx_path}")
    print(f"WROTE {csv_all}")
    print(f"WROTE {csv_global_pct}")
    print(f"WROTE {csv_global_abs}")
    print(f"WROTE {nb_path}")
    print("\nGLOBAL MAX ABSOLUTE PERCENT GAP")
    print(format_for_reading(global_pct).to_string(index=False))
    print("\nGLOBAL MAX ABSOLUTE METRIC-UNIT GAP")
    print(format_for_reading(global_abs).to_string(index=False))


if __name__ == "__main__":
    main()
