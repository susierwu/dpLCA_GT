# dpLCIA

This folder contains the dynamic prospective climate emission metrics pipeline and post processing used in the manuscript.

## Folder navigation

| Folder | Main purpose | How to use it |
|---|---|---|
| `FaIR_dpCFs` | Core metric generation from FaIR (scenario-dependent forcing + decay; module A/B/C workflow) | Start here to regenerate dpRF/dpAGWP and dpGWP inputs |
| `LCIA` | Implement the generated metrics in Brightway methods/datapackages | Use after metric generation for BW2 integration |
| `AGWPCO2_fixed_IPCCAR6` | Fixed-CO2 denominator pathway (`dpGWPfixed-CO2`) | Use for IPCC-anchored relative-metric calculations and comparison |
| `dp_approach_ModBvsC` | Dynamic prospective CO2 reference pathway (`dpGWPdp-CO2`, Eq. 11 denominator logic) | Use for scenario/MY-dependent denominator analyses under dynamic-reference formulation |
| `fixed_approach_ModBvcC` | Fixed CO2 reference pathway analyses (`dpGWPfixed-CO2`, Eq. 10 denominator logic) and old/new implementation comparison | Use for fixed-reference relative-metric calculations, sensitivity checks, and diagnostics |
| `Fig3-4_metrics` | Analysis and plotting notebooks/figures for manuscript metric interpretation | Primarily figure-generation and manuscript analysis |
| `pGWP100_bw2package` | Exported BW2 packages of GWP factors | Final distribution artifacts for BW2 usage |


## Recommended workflow order

1. Generate forcing/decay/metric outputs in `FaIR_dpCFs`.
2. Prepare and validate LCIA methods in `LCIA`.
3. Run fixed-vs-dynamic reference analyses (`AGWPCO2_fixed_IPCCAR6`, `dp_approach_ModBvsC`, `fixed_approach_ModBvcC`).
4. Produce manuscript-oriented analysis figures in `Fig3-4_metrics`.
5. Export reusable BW2 packages from `pGWP100_bw2package`.

## Manuscript mapping (renumbered figures)

The folder names were created during drafting. In the newest manuscript revision, earlier “Figure 3–4” analyses are now renumbered in the main text.

| Analysis content | Main folder | Manuscript figure mapping note |
|---|---|---|
| Absolute metric behavior (forcing/decay/AGWP trajectories) | `Fig3-4_metrics` + `FaIR_dpCFs/output` | Earlier Figure 3–4 workflow, now aligned to the renumbered Figure 4–5 sequence |
| Relative metric behavior (GWP variants; fixed-CO2 vs dp-CO2) | `AGWPCO2_fixed_IPCCAR6`, `dp_approach_ModBvsC`, `fixed_approach_ModBvcC` | Supports the updated relative-metric figures in the renumbered manuscript |
| Brightway implementation and case-study deployment | `LCIA` + downstream case folders | Method implementation layer; not just plotting |
