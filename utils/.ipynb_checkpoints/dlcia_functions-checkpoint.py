import re
import os
import sys
from lca_algebraic import initProject
import pandas as pd
import numpy as np
import xarray as xr
import lca_algebraic as agb
from lca_algebraic import *
from lca_algebraic.stats import * 
import bw2data, bw2io
import bw2calc
from premise import *
import premise as prm
from itertools import zip_longest
import random



def get_my_dpLCIA(ssp = 'SSP585',my = 2030, metric = 'dpIRF'): 
    """ to get a matched dpLCIA for each SSP[x], v_year[t] """
    sp = re.findall(r'\d+', ssp)[0] 
    bw2m_list = [m for m in bw2data.methods if sp in str(m) and  metric in str(m)  and str(my) in str(m)] 
    yy = random.randint(0, 100) 
    print(f"dpLCIA methods for SSP {sp}, MY{my} is, for example for year {yy},: {bw2m_list[yy]}, total len: {len(bw2m_list)}  ")
    return (bw2m_list)



def get_dlcia_rawtable_a (SSP_list,run_Y_list,tech_list, Year_list, GWP_list):
    """ once calc. done, below _a / _b / _c functions to cleanup data to prepare final ds"""  
    f_SSP_list = [item for sublist in SSP_list for item in sublist]
    f_run_Y_list = [item for sublist in run_Y_list  for item in sublist]
    f_tech_list = [item for sublist in tech_list for item in sublist]
    f_year_list = [item for sublist in Year_list for item in sublist]
    f_gwp_list =  [item for sublist in GWP_list  for item in sublist]
    
    """ 101 lcia score for each run_year, used as df.col, needs to be improved, linked to Time Horizon of LCIA, could be 101, 501 ... """ 
    #print(f"len of SSP/v_year/tech comb is: {len(f_SSP_list)}" )
    """ prepare final pandas DF , columns are yearly dlcia score, multiindex with SSP, v_year, tech_list, run_year """ 
    lciacol_names = list(range(0,101))
    dlcia = pd.DataFrame( f_gwp_list, columns = lciacol_names) 
    dlcia["SSP"] = f_SSP_list
    dlcia["v_year"] = f_year_list
    dlcia["tech_list"] = f_tech_list
    dlcia["run_year"] = f_run_Y_list
    dlcia.set_index(['SSP', 'v_year', 'tech_list', 'run_year'], inplace=True)
    
    return(dlcia)


def get_dlcia_ds_b (dlcia): 
    # Melt the DataFrame while keeping the MultiIndex intact
    dlcia_melt = pd.melt(dlcia.reset_index(), id_vars=dlcia.index.names, value_vars= dlcia.columns )
    dlcia_melt.set_index(['SSP', 'v_year', 'tech_list', 'run_year', 'variable'], inplace=True)
    # convert DF to xarrary.ds 
    dlcia_ds = dlcia_melt.to_xarray()
    dlcia_ds['variable'] = dlcia_ds['variable'].astype(int)
    #dgwp_ds
    dlcia_ds = dlcia_ds.rename({'variable': 'lcia_year'})
    return(dlcia_ds)



def get_dlcia_final_ds_empty_c (dlcia_ds,  v_module_touse = None, TH=100, v_lifetime = 20, b_replace_y = 10): 
    dpLCA_type = ['dpIRF', 'dpCRF', 'dpGWP']
    print (" enter one of the following for v_module to arrange dpLCIA results:  ['V1A_V1B_init'] for init_veh both ICEV&BEV,  \
         ['V1B']/ for BEV bat replacement, ['V2_oper'] for ICEV on-road emissions,  ['F2_fuel_diesel'] or ['F2_fuel_elec'] for fuel ") 
    SSP = dlcia_ds['SSP'].values
    v_year = dlcia_ds['v_year'].values
    tech_list = dlcia_ds['tech_list'].values  
    v_module = v_module_touse 
    ensems = np.arange(1001) #this is not used as we don't have uncertainty for dpLCIA yet, only ensemble = 0 has the final score stored 
    
    """ initial manuf. starts at YEAR[0], 100 is the LCIA_TH here, then operationa module start from Y[1] """
    """ battery replace year is special, it starts at Y[10], which is also the b_lifetime for simplicity """ 
    if 'V1A_V1B' in str(v_module): 
        lcia_year_tt = np.arange(0, TH+1, 1) 
    elif 'F2' in str(v_module) or 'V2' in str(v_module): 
        lcia_year_tt = np.arange(1, TH+1+v_lifetime, 1)  
    elif 'V1B' in str(v_module) and 'V1A' not in str(v_module):  
        lcia_year_tt = np.arange(b_replace_y, b_replace_y+TH+1)
    else: 
        print("check if you've entered the right v_module name")
    #print(lcia_year_tt)
     
    # Create a new xarray dataset with dimensions
    cord = {'SSP': SSP, 'v_year': v_year, 'v_module': v_module, 'tech_list':tech_list, 'lcia_year': lcia_year_tt, 'dpLCA_IC': dpLCA_type, 'ensemble': ensems}
    finalds_empty = xr.Dataset(coords=cord)
   
    fakedd = np.full((len(SSP), len(v_year), len(v_module), len(tech_list), len(lcia_year_tt), len(dpLCA_type), len(ensems) ), np.nan)    
    dplcia = xr.DataArray(fakedd, dims=('SSP', 'v_year', 'v_module', 'tech_list', 'lcia_year', 'dpLCA_IC', 'ensemble'), coords=cord )

    """ Why ds can't assign data variable with diff. values, I wanna have three seperate variable (w/o) dpLCA_IC dim,  """  
    #fakedd = np.full((len(SSP), len(v_year), len(v_module), len(tech_list), len(lcia_year_tt),  len(ensems) ), np.nan)
    #dpCRF = xr.DataArray(fakedd, dims=('SSP', 'v_year', 'v_module', 'tech_list', 'lcia_year', 'ensemble'), coords=cord )
    #dpGWP = xr.DataArray(fakedd, dims=('SSP', 'v_year', 'v_module', 'tech_list', 'lcia_year', 'ensemble'), coords=cord )
    
    # Assign the data variable to the dataset
    finalds_empty = finalds_empty.assign(dplcia = dplcia )
    print("empty final ds prepared with fake data")
    
    return(finalds_empty)






""" below functions to get consumption data on diesel / electricity """ 
############# if wanna check TRL together when getting annual consumption for each running_y, using the cloese neighbor year from TRL ds #########
def get_f2_fuel_elec_kwh_whTRL ( Y, v_lifetime , v_yearlyFU , SP , year_inP, v_pneed, data  ): 
    """ input: 
         v_yearlyFU could be either a float (same yearly FU), or a list e.g., [12000] * v_lifetime 
         SP -> ssp_to_get_inP , year_inP -> year_to_get_inP , v_pneed  will be supplied during calc 
         TRL data needed for yr_toget -> to filter through TRL so that each running[y] has a matched closest neighbor year to fetch tech 
         this won't affect the annualy consumption though, however, we assert there's some avail. tech for run_yeear[t]
    """ 
    yr_toget = min(data.v_year.values.tolist(), key=lambda x: abs(x - Y))
    print(f"vehicle running year {Y}, neighboring year to get tech {yr_toget}" )
    #ff = #data2.sel(powertrain = "BEV", Module = "F2_fuel_elec", v_type = "poubelle", SSP = ssp_to_get_inP, year = yr_toget)["tech_avail"].values
        
    ff = data.sel(Module = "F2_fuel_elec", SSP = SP, v_year = yr_toget, v_type = "garbage_truck", powertrain = "BEV") 
    ff = ff.where(ff["tech_avail"] == 1, drop = True)
    assert len(ff.Tech) != 0, f"We don't have any available tech to run your BEV for selected year {self.v_year} and {ssp}"
    assert len(v_pneed) != 0, f"didn't read in any foreground data, check the P_file, do not capitalize SSP"
        
    regex_pattern = re.compile("elec_consump_per") # it could be per km or per hour
    for key in v_pneed:
        if regex_pattern.search(key): 
            elec_consumperkm = v_pneed[key] 
            #print(f" 😳😳😳 for v_year{year_to_get_inP}, running_year {Y}, SSP {ssp_to_get_inP}, the elec_consumpperkm is {elec_consumperkm}")  
            break
        else:
            elec_consumperkm = 2.3
            
    elec_consumperkm2 = elec_consumperkm[0] if isinstance(elec_consumperkm, list) else elec_consumperkm  
    if isinstance(v_yearlyFU, list): 
        tt = Y - year_inP - 1 # [tt] is the list_index if v_yearlyFU is a list, -1 if as starts running from Y[1], not Y[0]
        elec_consumperyear = v_yearlyFU[tt]  * elec_consumperkm2
    else: 
        elec_consumperyear = v_yearlyFU  * elec_consumperkm2
    #print(f"  🔌🔌 for v_year{year_inP}, running_year {Y}, SSP {SP}, the elec_consumpperkm is {elec_consumperyear} annually ")  
    return(elec_consumperyear)


##### can combine it with get_f2_fuel_diesel_kg, adding a new input parameter for re.compile(" "), and change default consumperkm #### 
def get_f2_fuel_elec_kwh (Y, v_lifetime , v_yearlyFU , SP , year_inP, v_pneed ): 

    regex_pattern = re.compile("elec_consump_per") # it could be per km or per hour
    for key in v_pneed:
        if regex_pattern.search(key): 
            elec_consumperkm = v_pneed[key] 
            break
        else:
            elec_consumperkm = 2.3
            
    elec_consumperkm2 = elec_consumperkm[0] if isinstance(elec_consumperkm, list) else elec_consumperkm  
    if isinstance(v_yearlyFU, list): 
        tt = Y - year_inP - 1 # [tt] is the list_index if v_yearlyFU is a list, -1 if as starts running from Y[1], not Y[0]
        elec_consumperyear = v_yearlyFU[tt]  * elec_consumperkm2
    else: 
        elec_consumperyear = v_yearlyFU  * elec_consumperkm2
        
    #print(f" 🔌🔌 for v_year{year_inP}, running_year {Y}, SSP {SP}, the elec_consumpperkm is {elec_consumperyear} annually ")  
    return(elec_consumperyear)
    


def get_f2_fuel_diesel_kg(Y, v_lifetime , v_yearlyFU  , SP , year_inP, v_pneed, to_print_consmpt = False  ): 
    """ v_yearlyFU could be either a float (same yearly FU), or a list e.g., [12000] * v_lifetime   """
    # SP -> ssp_to_get_inP , year_inP -> year_to_get_inP ,   v_pneed  will be supplied during calc
    regex_pattern = re.compile("diesel_consump_per") # it could be per km or per hour
    for key in v_pneed:
        if regex_pattern.search(key): 
            d_consumperkm = v_pneed[key] 
            break
        else:
            d_consumperkm = 0.42 # kg per km a static default value
    d_consumperkm2 = d_consumperkm[0] if isinstance(d_consumperkm, list) else d_consumperkm  
    if isinstance(v_yearlyFU, list): 
        tt = Y - year_inP - 1 # [tt] is the list_index if v_yearlyFU is a list, -1 if it starts running from Y[1], not Y[0]
        d_consumperyear = v_yearlyFU[tt]  * d_consumperkm2
    else: 
        d_consumperyear = v_yearlyFU  * d_consumperkm2

    if to_print_consmpt == True:
        print(f"  ⛽️⛽️ for v_year{year_inP}, running_year {Y}, SSP {SP}, the annual diesel_consump is {d_consumperyear}")  
    else:
        pass
    return(d_consumperyear)

