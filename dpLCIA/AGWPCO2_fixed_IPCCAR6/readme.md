now in the 2026/Feb update, we need to rerun following notebooks: 

- `2_dpGWP_CO2CH4N2O_per_fixedAGWPCO2` 


no need to rerun 1a/1b for metric output, they're just testing files (though they're also updated and tested in 1a_v2026 / 1b_v2026)
- `1a_dpGWP_CH4N2O_per_fixedAGWPCO2.ipynb` : the CH4/N2O dpAGWP metric folder now is with new data
- `1b_dpGWP_CO2_per_fixedAGWPCO2` same as above change

though, keep using the AGWP[CO2] unchanged as in `0_agwp_CO2_fixed_SSP.ipynb`, caz it's per IPCC not changing d/p



then the new output (only pRE, with IPCC impulse response function) is: 
- CF_v2026_GWP0_100_perSSP_MY_majorghgs.nc
- CF_v2026_GWP1_100_perSSP_MY_majorghgs.nc
- dpGWP_v2026_CO2CH4N2O_output_fixed_AGWPCO2.xlsx



or, new output ( pRE +  prospeoctive decay factor) to be used in further dpLCIA preparation in BW2 is: 
- CF_v2026_use4BW2_whModC_GWP0_100_perSSP_MY_majorghgs.nc
- CF_v2026_use4BW2_whModC_GWP1_100_perSSP_MY_majorghgs.nc
- dpGWP_v2026_use4BW2_whModC_CO2CH4N2O_output_fixed_AGWPCO2.xlsx




Benchmark two metric approach:

1. pRE with the IPCC impulse response function (as the old module B) 
2. pRE with new prospeoctive decay factor to replace the IPCC impulse response function (as new module C for metric calculation) 

benchmark notebook: `3_v2026_dpGWP_whpIRF_CO2CH4N2O_per_fixedAGWPCO2`


we also had two seperate new folders `AGWPCO2_dp_approach` and `AGWPCO2_fixed_approach` for unified metric calculation (gwp / agwp / irf) to replace old three notebooks  under "LCIA/1_CFanal_CRF.ipynb | 1_CFanal_GWP.ipynb | 1_CFanal_IRF.ipynb"  