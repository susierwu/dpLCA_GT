import re
import os
import sys
import pandas as pd
import numpy as np
import xarray as xr
from itertools import zip_longest
import random
import bw2data, bw2io
import bw2calc


def get_my_pGWP100(ssp = 'SSP585', my = 2030, metric = 'pGWP'): 
    """ to get a matched dpLCIA for each SSP[x], v_year[t] """
    sp = re.findall(r'\d+', ssp)[0] 
    bw2m_list = [m for m in bw2data.methods if sp in str(m) and  metric in str(m) and str(my) in str(m)]
    if len(bw2m_list) == 1: 
        bw2m_list = bw2m_list[0]
    else: 
        print("more than one pGWP")
    return (bw2m_list)


def get_f2_fuel_elec_kwh_diesel_kg(powertrain, Y, v_lifetime, v_yearlyFU, SP, year_inP, v_pneed, to_print_consmpt = False  ): 
    """ 
    to get annual fuel consumption from consumption parameter file 
    Inputs: 
    -------
    Y: each running year, supplied during calculation
    v_yearlyFU: yearly running mileage, default 36000km 
    
    Return: 
    -------
    f_consumperyear: annual consumption data
    """
    if powertrain == "BEV": 
        regex_pattern = re.compile("elec_consump_per")  
        fuel_consumperkm = 2.3 # eelc kwh per km a static default value for HDV-GT
        f_unit = "kWh"
    elif powertrain == "ICEVd": 
        regex_pattern = re.compile("diesel_consump_per")
        fuel_consumperkm = 0.42 # diesel kg per km a static default value for HDV-GT
        f_unit = "kg"
    else:
        pass
        print("to add hybrid fuel")

    for key in v_pneed:
        if regex_pattern.search(key): 
            f_consumperkm = v_pneed[key] 
            break
        else:
            f_consumperkm = fuel_consumperkm 
            
    f_consumperkm2 = f_consumperkm[0] if isinstance(f_consumperkm, list) else f_consumperkm  
    if isinstance(v_yearlyFU, list): 
        tt = Y - year_inP - 1 # [tt] is the list_index if v_yearlyFU is a list, -1 if it starts running from Y[1], not Y[0]
        f_consumperyear = v_yearlyFU[tt]  *  f_consumperkm2
    else: 
        f_consumperyear = v_yearlyFU  * f_consumperkm2

    
    if to_print_consmpt == True:
        print(f"  ⛽️ 🔌 for v_year{year_inP}, running_year {Y}, SSP {SP}, the annual fuel_consump is {f_consumperyear} {f_unit} ")   
    else:
        pass

    return(f_consumperyear)




def get_plcia_wh_v_mod_icev_f2 (SSP_list, run_Y_list,tech_list, Year_list, GWP_list, v_module = "F2_fuel_diesel", lciacol_names = ['IPCC2021_default', "premise_GWP",  'pGWP100' ]):
    """  
    when there's no TRL engaged, for F2 with annual dynamic inventory, run this function to get final lifecycle impacts
    Inputs: 
    ---------
    lciacol_names: impact category names, e.g., for GWP, three impacts hard-coded
    
    Returns: 
    ---------
    final_df: lifetime impact, used for preparing final ds in function: get_plcia_ds()
    dlcia: annual impact, using dynamic inventory
    """  
    f_SSP_list = [item for sublist in SSP_list for item in sublist]
    f_run_Y_list = [item for sublist in run_Y_list  for item in sublist]
    f_tech_list = [item for sublist in tech_list for item in sublist]
    f_year_list = [item for sublist in Year_list for item in sublist]
    f_gwp_list =  [item for sublist in GWP_list  for item in sublist]
    #lciacol_names = ['IPCC2021_default', "premise_GWP",  'pGWP100' ] #list(range(0,len(methods_touse) ) )
    plcia = pd.DataFrame( f_gwp_list, columns = lciacol_names) 
    plcia["SSP"] = f_SSP_list
    plcia["v_year"] = f_year_list
    plcia["tech_list"] = f_tech_list
    plcia["run_year"] = f_run_Y_list 
    plcia["v_module"] = [v_module] * len(f_SSP_list)
    plcia.set_index([ 'SSP', 'v_module', 'v_year', 'tech_list', 'run_year'], inplace=True) 
    final_df = plcia.groupby(level=['SSP', 'v_module', 'v_year',  'tech_list']).sum()
    return(final_df, plcia)


def get_plcia_wh_v_mod_bev_grid (SSP_list, run_Y_list,tech_list, Year_list, GWP_list, v_module = "F2_fuel_elec", lciacol_names = ['IPCC2021_default', "premise_GWP",  'pGWP100' ]):
    """  
    for electricity, using grid as default, no PV (as it involved TRL for diffrent running years) for final aggregation  
    Inputs: 
    ---------
    lciacol_names: impact category names, e.g., for GWP, three impacts hard-coded
    Returns: 
    ---------
    final_df: lifetime impact, from plcia_filter, with only grid tech 
    plcia_filter: annual impacts, only grid tech 
    plcia: annual impact, using dynamic inventory, unfiltered, for each v_year, there're some running years with selected tech not avaiable 
    """  
    f_SSP_list = [item for sublist in SSP_list for item in sublist]
    f_run_Y_list = [item for sublist in run_Y_list  for item in sublist]
    f_tech_list = [item for sublist in tech_list for item in sublist]
    f_year_list = [item for sublist in Year_list for item in sublist]
    f_gwp_list =  [item for sublist in GWP_list  for item in sublist]

    plcia = pd.DataFrame( f_gwp_list, columns = lciacol_names) 
    plcia["SSP"] = f_SSP_list
    plcia["v_year"] = f_year_list
    plcia["tech_list"] = f_tech_list
    plcia["run_year"] = f_run_Y_list 
    plcia["v_module"] = [v_module] * len(f_SSP_list)
    # filter only grid elec, for final lifetime impacts sum 
    plcia_filter = plcia[plcia['tech_list'] == 'grid']
    plcia_filter.set_index([ 'SSP', 'v_module', 'v_year', 'tech_list', 'run_year'], inplace=True) 
    final_df = plcia_filter.groupby(level=['SSP', 'v_module', 'v_year',  'tech_list']).sum()    
    return(final_df,  plcia_filter, plcia )



def get_plcia_wh_V1AV1B_BEV (SSP_list, year_list, tech_list, GWP0_list, GWP1_list, GWP2_list,
                        default_tech = "NMC622", v_module = "V1AB_init", lciacol_names = ['IPCC2021_default', "premise_GWP",  'pGWP100' ]):      
    plcia = pd.DataFrame(zip(GWP0_list, GWP1_list, GWP2_list), columns = lciacol_names) 
    plcia["SSP"] = SSP_list
    plcia["v_year"] = year_list
    plcia["tech_list"] = tech_list
    plcia["v_module"] = [v_module] * len(SSP_list)
    plcia_default_tech = plcia[plcia['tech_list'] == default_tech]
    plcia_default_tech.set_index([ 'SSP', 'v_module', 'v_year', 'tech_list'], inplace=True) 
    plcia_default_tech = plcia_default_tech.drop_duplicates()
    # .set_index after preparing DF for default_tech
    plcia.set_index([ 'SSP', 'v_module', 'v_year', 'tech_list'], inplace=True) 
    plcia_pertech = plcia.drop_duplicates()
    return(plcia_default_tech, plcia_pertech)


def get_plcia_wh_V1A_ICEV (SSP_list, year_list, tech_list, GWP0_list, GWP1_list, GWP2_list,
                          v_module = "V1AB_init", lciacol_names = ['IPCC2021_default', "premise_GWP",  'pGWP100' ]):      
    # keep the v_module name same as BEV for "V1AB_init", for final benchmark if needed 
    plcia = pd.DataFrame(zip(GWP0_list, GWP1_list, GWP2_list), columns = lciacol_names) 
    plcia["SSP"] = SSP_list
    plcia["v_year"] = year_list
    plcia["tech_list"] = tech_list
    plcia["v_module"] = [v_module] * len(SSP_list)
    plcia.set_index(['SSP', 'v_module', 'v_year', 'tech_list'], inplace=True) 
    return(plcia)



def get_plcia_ds (final_df): 
    plcia_melt = pd.melt(final_df.reset_index(), id_vars=final_df.index.names, value_vars= final_df.columns )
    plcia_melt.set_index(['SSP', 'v_module', 'v_year', 'tech_list', 'variable'], inplace=True) 
    plcia_ds = plcia_melt.to_xarray()
    plcia_ds = plcia_ds.rename({'variable': 'Impact_category'})
    return(plcia_ds)


