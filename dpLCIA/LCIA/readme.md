# LCIA (Brightway implementation layer)

This folder is the implementation layer that converts generated dynamic prospective metric outputs into Brightway-compatible LCIA methods and export artifacts.

## Scope and caution

- This folder is primarily for **BW2 method construction/export workflows**.
- For the current manuscript, prioritize **v2026** workflows.
- Most notebooks without `v2026` in the name are legacy and tied to older metric logic.

## Folder navigation

| Path | Main purpose | Status | How to use it |
|---|---|---|---|
| `data/v2026` | NetCDF inputs used by v2026 LCIA creation workflows | Primary | Use as source data for v2026 notebooks |
| `2_createBW2_A_add_dpIRF_majorGHG.ipynb` | Build BW2 LCIA methods for the A-side (major-GHG IRF-related method step) | Legacy-to-transition | Use only if your pipeline explicitly needs this step |
| `3_v2026_createBW2_B_add_pGWP100_whminorGHG.ipynb` | Main v2026 method-building notebook for prospective GWP100 (incl. minor GHG integration) | Primary | One of the key notebooks to run |
| `6a_v2026_export_pGWP100.ipynb` | Export v2026-built methods to BW2 datapackage artifacts | Primary | Run after method creation |
| `6b_test_import_pGWP100.ipynb` | Validate datapackage import into BW project | Primary QA | Use as post-export check |
| `fixed_AGWPCO2_createBW2/` | Legacy fixed-CO2 BW2 creation path from earlier workflow versions | Legacy | Keep for historical comparison only |
| `minorGHG/` | Supporting notebooks for minor-GHG metric preparation | Supporting | Use when updating/rebuilding minor-GHG components |
| `WIP/` | Work-in-progress notebooks and experiments | Development | Do not treat as stable production path |
| `output_plot/` | Plotting artifacts related to LCIA-side checks | Supporting | Optional diagnostics only |
| `premise_gwp/` | Premise-linked GWP support files and mappings | Supporting | Needed by selected method-assembly steps |

## Recommended run order (v2026)

1. Confirm required v2026 `.nc` inputs are present in `data/v2026`.
2. Run `3_v2026_createBW2_B_add_pGWP100_whminorGHG.ipynb`.
3. Run `6a_v2026_export_pGWP100.ipynb` to export datapackage(s).
4. Run `6b_test_import_pGWP100.ipynb` to verify import and method availability.

## Relationship to other folders

| Upstream/downstream | Folder | Relation |
|---|---|---|
| Upstream | `../FaIR_dpCFs` | Generates metric inputs consumed by LCIA method construction |
| Parallel analysis | `../AGWPCO2_fixed_IPCCAR6` and `../dp_approach_ModBvsC` | Produces denominator/reference-side analyses that inform metric interpretation |
| Downstream packaging | `../pGWP100_bw2package` | Receives exported BW2-compatible datapackage artifacts |

## Legacy note

- Older notebooks in this folder may reflect the previous pathway that coupled FaIR forcing with IPCC-style analytical impulse response assumptions.
- For current manuscript reproduction, use the **v2026** notebooks and data only unless you are explicitly doing old-vs-new method comparison.
