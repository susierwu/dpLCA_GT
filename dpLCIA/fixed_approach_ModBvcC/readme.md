## fixed GWP/AGWP/RF metric calculation: 


now we can run all metric calculation in one "*_unified*" notebook, though AGWP / GWP are kept as two notebook with only the "metric_value" changed during configuration, so we can see rendered plot outputs in each notebook. 




when EXPORT_OUTPUTS == TRUE, 
the output folder `output_dpGWP_fixedCO2_modC` contains data (.nc) needed for later prepare BW2 LCIA method, which we'll copy to the "../LCIA/data" folder  
we DO NOT export data with full 1001 ensemble, unless it's needed for BW2 calculation 

for now, we don't have yet exported all figure output 


### important note, 
#### the fixedCO2 vs. dpCO2 output has different nrows:: if we include all 6 SSPs, 3 MYs, then under fixedCO2 approach, we only have 1818 rows output, but under dpCO2 approach, we have 3618 rows output, this is because in the new runs of FaIR, we set len = 200 years, however, under fixed CO2 approach, we need the fixedCO2 AGWP from `AGWPCO2_fixed_IPCCAR6` folder where we only run 100 years.  so to make length consistent (H=100) for BW2 LCIA of both dp vs. fixedCO2 approach, use the "../dpCO2_unified_metric_GWP_ModBvsC_saveFig_DTlen101_forBW2" for export ModuleC GWP for dpCO2, caz there we updated `export_metric_outputs()` which trim to the configured horizon before saving 