
###  to navigate the folders: 

The dpIRF and dpCRF (aka dpAGWP) for CO2, CH4, N2O are computed in `FaIR_dpCFs` folder, then the new metrics are prepared as BW2 LCIA methods under the `LCIA` folder, for the dpIRF and the dpGWP.  

Two approaches to compute the relative metric - dpGWP

1.  dpGWPfixed-CO2 approach: fixed reference substance approach: fixing the AGWP[CO2] as in the IPCC_AR6, in the `AGWPCO2_fixed_IPCCAR6` folder
2.  dpGWPdy-CO2 approach: dynamic prospective reference substance approach: the denominator AGWP[CO2] itself is dynamically and prospectively changing as the nominator, in the `FaIR_dpCFs` folder 

The dpGWP100 Characterization factors from the two approach for computing dpGWP are further exported as BW2 datapackage, in the `pGWP100_bw2package` folder 
 


