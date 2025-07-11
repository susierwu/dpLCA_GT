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

# like in dLCA, calculate the V2 exhaust emission directly w/o creating more UPRs in the database
# for CO2, premise_gwp and new pGWP share the same UUID, but not for IPCC 2021, 
default_ipcc = ('IPCC 2021 no LT', 'climate change no LT', 'global warming potential (GWP100) no LT')
premise_gwp = ('IPCC 2021', 'climate change', 'GWP 100a, incl. H and bio CO2')

map_majorghg_cfs = { "CO2": ["d6235194-e4e6-4548-bfa3-ac095131aef4", 1],
             "CH4": ["7f0ba7c9-341e-413d-80f6-8753727d65d1", 0], 
             "N2O": ["473826ae-125a-4b02-8c8e-c84322491d80", 0] 
          }


def get_my_pGWP100(ssp = 'SSP585', my = 2030, metric = 'pGWP'): 
    """ to get a matched dpLCIA for each SSP[x], v_year[t] """
    sp = re.findall(r'\d+', ssp)[0] 
    bw2m_list = [m for m in bw2data.methods if sp in str(m) and  metric in str(m)  and str(my) in str(m)]
    if len(bw2m_list) == 1: 
        bw2m_list = bw2m_list[0]
    else: 
        print("more than one pGWP")
    return (bw2m_list)


def load_bw2mm_map_cf_for_CFdict (bw2_mm, map_majorghg_cfs) : 
    """ 
    Inputs: 
    --------
    bw2_mm: the bw2-method for LCIA, e.g., IPCC default, pGWP100  
    map_majorghg_cfs: the empty dict with 0 CFs
    Return: 
    --------
    map_majorghg_cfs: CF_dict with cf values assigned, used for calc_co2e_majorghg_v2()
    """
    #bw2_mm: above bw2m_list for SSP[x], MY[t]
    cfs = bw2data.Method(bw2_mm).load()
    for gas in map_majorghg_cfs.keys():
        # for CO2, premise_gwp and new pGWP share the same UUID, but not for IPCC 2021
        if gas == "CO2": 
            cf_value = 1 
        else: 
            for i in cfs: 
                if i[0][1] == map_majorghg_cfs.get( gas )[0]: 
                    if isinstance(i[1], (int, float)):
                        cf_value = i[1] 
                    elif isinstance(i[1].item(), (int, float)) :
                        cf_value = i[1].item()
                    else:   #isinstance(i[1].item(), dict): # in case dict with uncertainty
                        try: 
                            cf_value = i[1].get("amount")
                        except: 
                            pass
                            print(f"no CF value found for gas {gas}")
                #print(f"for {gas}, cf turple is {i}, cf_value is {cf_value}" )
        # assign new cf_value to the dict, 2nd element of the value
        map_majorghg_cfs[gas][1] = cf_value
    return map_majorghg_cfs



def calc_co2e_majorghg_v2 (lci_df, cf_dict): 
    """ 
    Inputs: 
    --------
    lci_df: dataframe for certain v_year, annual emis, will sum up for final lifecycle emission, define lci_df during calculation for each ssp/v_year lci_df = LCI_dt_dict.get('vy2030')
    cf_dict: map_majorghg_cfs, the CF dict assigned with cf_values f
    Return: 
    --------
    total_co2e: final lifecycle LCA score for v_year[t], SSP[x] vehicle  
    """
    # summing all annual emission from lci_df for each gas 
    sum_lci_df = pd.DataFrame(lci_df.sum(axis=0))
    
    lca_list = [] 
    for gas in sum_lci_df.index: 
        if gas != 'year': 
            cf_gas = cf_dict.get(gas)
            cf_gas_value = cf_gas[1]
            emission = sum_lci_df.loc[gas].values[0]
            print(f"for {gas}, lifecycle emission kg are: {emission}, using CF of {cf_gas} ") 
            lca_score = emission * cf_gas_value
            lca_list.append(lca_score)
    total_co2e = sum(lca_list)
    return(total_co2e)

