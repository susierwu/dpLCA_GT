# Elec_case

Electricity case study workspace for the current manuscript.

This folder now has two parallel tracks:

- **Primary track:** `dpLCI_v2026` using updated metrics linked to **Module C** (prospective forcing + prospective decay from FaIR).
- **Legacy track:** `dynLCI_legacy-ModB_metric` using **IPCC-style analytical impulse response functions** (Module B pathway).

## Folder navigation

| Folder | Main purpose | Status | How to use it |
|---|---|---|---|
| `dpLCI_v2026` | Main electricity workflow for dynamic prospective LCI + dpLCIA in the revised manuscript | **Primary** | For the case study |
| `dpLCI_v2026/dpLCIA_hydro` | Hydro-focused dpLCIA workflow (date shift, compute/save, and RF/AGWP plotting/decomposition) | **Primary** | Use for hydro deep-dive analysis |
| `dpLCI_v2026/plot` | Manuscript and SI plotting/analysis notebooks (rank shift, split-MY, gap plots) | **Primary** | Run after generating case outputs |
| `dynLCI_legacy-ModB_metric` | Legacy electricity workflow built around old IPCC-style analytical IRF assumptions | Legacy/reference | Use only for back-comparison with old approach |
| `set_up` | Environment/database/method setup (BW25, ecoinvent/premise, method import/export) | **Primary** Shared dependency | Start here and run before case study notebooks |
| `elec_share` | Electricity-mix/share preprocessing and support inputs | Shared dependency | Use when scenario electricity shares need refresh/check |
| `utils` | Utilities for case-study processing and helper functions | Shared dependency | Imported by notebooks in this folder |


## Recommended run order for reproduction

1. **Setup**
   - Use notebooks/scripts in `set_up` to prepare BW25 project, background DB, and metric methods.
2. **Generate dynamic inventories and impacts**
   - Run main technology notebooks in `dpLCI_v2026`:
     - `dynElec_hydro.ipynb`
     - `dynElec_wind.ipynb`
     - `dynElec_pv.ipynb`
3. **Hydro metric deep dive (if needed)**
   - Run `dpLCI_v2026/dpLCIA_hydro` sequence:
     - `1_dpLCIA_dateshift.ipynb`
     - `2_calc_save_dpLCIA.ipynb`
     - `3_plot_rf_agwp_lines.ipynb`
     - `4_plot_rf_single-flow-decomposed.ipynb`
4. **Post-analysis and figures**
   - Run notebooks in `dpLCI_v2026/plot` for rank-shift, gap, and split-MY analyses.
5. **Legacy comparison (optional)**
   - Only then run `dynLCI_legacy-ModB_metric` notebooks for old-vs-new comparison.



## Outputs and reproducibility

- Large generated outputs (including some `.pkl`, `.nc`, or figure-heavy artifacts) may be excluded from Git tracking.
- If required files are missing, regenerate by rerunning notebooks in order above.
- Legacy and primary tracks should not be mixed in the same result table without clearly labeling metric assumptions.

