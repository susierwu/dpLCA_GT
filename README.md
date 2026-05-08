# dpLCA_GT

Repository for dynamic/prospective LCA workflow development.

## Scope update in 2026 revision

- **Primary project**: `Elec_case` (new electricity case study)
- **Methods library**: `dpLCIA/FaIR_dpCFs`  
- **Legacy side project**: `GT_case` (garbage truck case)



## Practical environment note

- The repository does not yet include a single pinned environment file (`requirements.txt`/`environment.yml`) at root.
- If you want strict reproducibility, create environment spec after confirming package versions in your working environment shown below:




## Folder map

| Folder | Role in repo | Priority |
|---|---|---|
| `Elec_case` | **Main case study for current manuscript** (electricity-focused workflows, notebooks, scenario analyses) | High |
| `dpLCIA` | Dynamic LCIA method and metric generation framework (FaIR + LCIA preparation) | High |
| `utils` | Shared helper functions used earlier for the old case study calculation | Low |
| `GT_case` | Garbage-truck case-study materials kept for reference/side analyses | Low |
| `image` | Static image assets (legacy diagram retained as archival reference only) | Low |



## Recommended navigation

| Step | Folder | What is inside | What to run / read first |
|---|---|---|---|
| 1 | `dpLCIA/FaIR_dpCFs` | FaIR-based metric generation (RF/AGWP/GWP-related dynamic CF workflow) | Start here if rerunning FaIR (we used v2.1.2) or just need new metric output files |
| 2 | `dpLCIA/LCIA` | LCIA method preparation under brightway2/25 and export pipeline | Use after new metric generation |
| 3 | `Elec_case/set_up` | Environment and setup helpers for electricity workflow | Start here when setting up a fresh run |
| 4 | `Elec_case/dpLCI_v2026` | Main electricity-case notebooks/data pipeline for dynamic LCI and downstream analysis | Entry point for manuscript reproduction |
| 5 | `Elec_case/elec_share` | Electricity mix/share inputs and scenario support files | Use when checking SSP assumptions and electricity trajectories |
| 6 | `Elec_case/utils` | Electricity-case-specific helper functions/scripts | Refer when notebooks call local utility modules |
| 7 | `Elec_case/dpLCI_v2026/dpLCIA_hydro` | Hydro-focused dpRF/dpAGWP calculations, saved outputs, and plotting notebooks | Use for hydro benchmark tables/figures in manuscript |
| 8 | `utils` | Shared cross-project utilities | Read when imports point to top-level `utils` |





### Dependencies and versions

Below are the core dependencies needed to run the current workflow.

| Package / data | Role in workflow | Version |
|---|---|---|
| `FaIR` | Climate model used to generate scenario-dependent forcing and prospective decay kernels | `v2.1.2` |
| `bw2data`, `bw2io`, `bw2calc` | Core Brightway libraries used across setup, LCIA method building, and calculations | `v4.5.3 / v0.9.11 / v2.2.2` |
| `brightway25` | Compatibility meta-package used in setup notebooks (`import brightway25`) | `v1.1.0` |
| `bw_timex` | Dynamic timeline-based LCA calculations in electricity legacy notebooks and related workflows | `v0.3.3` |
| `premise` | Prospective IAM-linked background database generation (ecoinvent transformations) | `v2.3.2` |
| `ecoinvent` database files | Source LCI database used in BW projects | `ecoinvent 3.9.1` for legacy case study `ecoinvent-3.11` for the new case study |
| `numpy`, `pandas`, `xarray` | Numerical arrays, tabular handling, and multidimensional data processing | `v1.26.4 / v2.3.3 / v2025.11.0` |
| `matplotlib`, `seaborn` | Plotting and manuscript-figure generation | `v3.10.7 / v0.13.2` |
| `openpyxl` | Excel read/write in utility and analysis scripts | `v3.1.5` |
| `scipy` | Scientific helpers in selected metric notebooks | `v1.16.3` |



## Data and reproducibility notes

### .pkl / .nc files 


Many `.pkl` outputs are excluded from the GitHub repository for two reasons:

1. File size constraints (commonly larger than GitHub limits, e.g., `>25 MB`).
2. Some pickle files contain calculated outputs tied to **ecoinvent-derived** intermediate/final results.

For reproducibility, users should regenerate these outputs by running the notebooks in sequence.

heavy files (for example large `.nc` outputs in some workflows) may also be absent from the repo and should be regenerated locally unless separately provided.
