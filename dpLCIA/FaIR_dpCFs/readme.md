# How we use FaIR in this project

This folder stores notebooks, modules, and outputs used to derive dynamic prospective climate metrics.

## Output folders and which one to use

Use this section first when you only need to know which files to consume.

| Output path | Workflow type | Main metrics in files | Name in excel/pickle files (Table S1.2 style) | Typical use |
|---|---|---|---|---|
| `output/metrics_AGWP_wh_ModC_CC_T/` | **Updated pathway** (prospective forcing + prospective decay, Module C) | RF, AGWP, AGTP, iAGTP, GWP, and CC-adjusted variants for non-CO2 gases | `rf`, `agwp`, `agtp`, `iagtp`, `gwp`, `rf_cc`, `agwp_cc`, `rf_final`, `agwp_final` | **Primary outputs for current manuscript interpretation** |
| `output/metrics_2026-02-23/` | Legacy pathway (prospective forcing + analytical IPCC-style decay, Module B) | RF/AGWP/GWP-oriented outputs from old approach | legacy naming in `agwp_dcf_gwp100_tstep1...xlsx` files | Benchmark/comparison only |
| `output/pIRF/*compare*/` and `output/*comp*` notebooks | Diagnostics and benchmarking | Comparison plots/tables (including cross-study checks) | n/a | Supporting analysis / SI diagnostics |

### Quick interpretation of acronym

| Name in files | Meaning |
|---|---|
| `rf` | Instantaneous radiative forcing |
| `agwp` | Absolute Global Warming Potential (cumulative RF) |
| `agtp` | Absolute Global Temperature Change Potential |
| `iagtp` | Integrated AGTP |
| `gwp` | Global Warming Potential |
| `rf_cc` (non-CO2) | Carbon-cycle feedback adjustment to RF |
| `agwp_cc` (non-CO2) | Carbon-cycle feedback adjustment to AGWP |
| `rf_final` (non-CO2) | Final RF = direct RF + carbon-cycle feedback RF |
| `agwp_final` (non-CO2) | Final AGWP = direct AGWP + carbon-cycle feedback AGWP |


## What changed in the 2026 revision

- We explicitly separate two pathways (as shown in the flowchart figure):
  - **Older pathway**: prospective forcing + analytical IPCC-style decay (Module A + B)
  - **Updated pathway**: prospective forcing + **prospective decay factor** from FaIR pulse reruns (Module A.2 + C)
- The updated pathway is the core physics used for the new metric interpretation in the manuscript.

## Method map to manuscript equations (eq. 4-11)

| Manuscript equation block | Physical meaning | Where implemented in this folder |
|---|---|---|
| Eq. 4 / 4a | Scenario- and ModelYear-dependent marginal ERF from FaIR background state | `utils/majorghg_modA.py` (forcing-side extraction), then consumed by module B/C |
| Eq. 5 | Pulse-perturbed minus baseline concentration response (`dC`) | `utils/majorghg_modA_pulse_rerun-fair.py` |
| Eq. 6 | Prospective decay factor from normalized `dC` | `utils/majorghg_modA_pulse_rerun-fair.py` (pIRF) |
| Eq. 7 | Absolute forcing kernel combining ERF amplitude and decay kernel | `utils/majorghg_modC_replaceB.py` |
| Eq. 8-9 | Time integration to absolute metrics (dpRF, dpAGWP) in discrete annual implementation | `utils/majorghg_modC_replaceB.py` |
| Eq. 10 | Relative metric with fixed CO2 reference (2019/IPCC-style denominator) | Computed downstream in fixed-reference workflow (`../AGWPCO2_fixed_IPCCAR6`) |
| Eq. 11 | Relative metric with dynamic prospective CO2 reference (same SSP/MY background as numerator) | Computed downstream in dynamic-reference workflow (`../dp_approach_ModBvsC`) |



## Overall flowchart

![FaIR Module A/B/C function flowchart](utils/FaIR_ModuleABC_function_flowchart.jpg)



## Module architecture used here

| Module | Role | Status |
|---|---|---|
| Module A.1 (`majorghg_modA.py`) | Extract FaIR-based forcing-side quantities and ERF-per-unit-perturbation | Legacy-compatible and still useful for comparison |
| Module A.2 (`majorghg_modA_pulse_rerun-fair.py`) | Rerun FaIR with pulse perturbation and derive prospective decay factors (pIRF) | **Key updated component** |
| Module B (`majorghg_modB_whRF.py`) | Analytical IPCC-style decay-based metric calculations | Legacy/reference pathway |
| Module C (`majorghg_modC_replaceB.py`) | Replace analytical decay with FaIR-derived prospective decay; compute prospective metrics numerically | **Primary updated pathway** |


### `utils/` folder 

| File | What it is for |
|---|---|
| `majorghg_modA.py` | Forcing-side extractor for major GHGs (CO2/CH4/N2O); computes scenario/MY-dependent ERF terms from FaIR state variables. |
| `majorghg_modA_pulse_rerun-fair.py` | Pulse-rerun engine to compute `dC` and normalized prospective decay (pIRF); foundation for the updated Module C pathway. |
| `majorghg_modB_whRF.py` | Legacy analytical metric engine (IPCC-style decay forms); computes RF/AGWP/GWP with optional carbon-cycle adjustment terms. |
| `majorghg_modC_replaceB.py` | Updated numerical metric engine using FaIR-derived pIRF + ERF; computes RF/AGWP/AGTP/iAGTP and optional carbon-cycle feedback adjustments. |
| `minorghg_modAB_combined.py` | Minor-GHG metric calculation utilities, combining parameter extraction and analytical metric evaluation with FaIR species data. |
| `majorghg_extra_analy.py` | Post-processing/diagnostic helper class for collecting and plotting rolling-year metric outputs across SSP and ModelYear slices. |
| `premise_gwp_name_map_TBD.txt` | Working name-mapping scaffold for minor GHG species between FaIR naming and premise/BW2 naming conventions. |
| `FaIR_ModuleABC_function_flowchart.jpg` | Current module flow diagram (A/B/C architecture). |
| `ModuleAB_function_flowchart.jpg` | Older flowchart reflecting pre-Module-C workflow. |
| `pIRF_test/` | Development/testing scripts for pIRF and module variants during method iteration. |



### implementation guide (concept -> code)

| Concept | Symbol | Implementation pointer |
|---|---|---|
| Scenario | SSP | Class init `scn=...` in module A/B/C |
| Calendar year in FaIR | `y` | FaIR `timebounds/timepoints` and `year_index` handling |
| ModelYear | `MY` | Selected via `year_index` (background state selector) |
| Metric time | `tau` | `H` arrays (`0..H_max`) used for RF/AGWP trajectories |
| Assessment horizon | `H` | `H_max` and `ts_per_year` in module classes |
| Forcing-side ERF trajectory | `dERFi` | Module A forcing extraction |
| Forcing amplitude | `dERFi(MY|SSP)` | selected at `erf_diff_t[year_index]` |
| Analytical decay option | IPCC analytical IRF | Module B (`majorghg_modB_whRF.py`) |
| Prospective decay option | Prospective decay factor | Module A.2 (`majorghg_modA_pulse_rerun-fair.py`) |
| Absolute metrics | dpRF / dpAGWP | Module C (`majorghg_modC_replaceB.py`) |
| Carbon-cycle feedback adjustment | CC adjustment | Module B/C adjustment routines |


#### Prospective decay 

The prospective decay curve is obtained by rerunning FaIR with a unit pulse perturbation and normalizing by the initial post-pulse concentration jump:

Detailed note and derivation, see:`prospective-decay-rerunFaIR_perturbation.md`


## Notes

- Scenario set commonly used: `ssp119`, `ssp126`, `ssp245`, `ssp434`, `ssp460`, `ssp585`.
- Under strong forcing backgrounds, the normalized prospective decay can exceed 1 at later metric times for CO2 because decay is state-dependent (not fixed analytical IRF behavior).
