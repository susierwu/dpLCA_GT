import os
import pandas as pd
import numpy as np
from numpy import nan
import pickle
from itertools import product 
from IPython.display import display

PROJECT_NAME = "iveo_v1"    
#eidb = "ecoinvent 3.9.1"
BIO_DB = "biosphere3"
PDB_NAME = "iveo_Parameterized_v1"
#truck_DB = "truck"
YEARS = [2040, 2045, 2050]
YEARS_STR = list(str(year) for year in YEARS)
#PARAMS_FILE = "data/p_scn_test.xlsx"
hypothesis = "point value"
 
ssp_remind_Pname_map = {
                        "SSP1-PkBudg500" : "ssp119", 
                        #"SSP1-PkBudg1150": "ssp126", 
                        "SSP2-Base": "ssp245", 
                        "SSP5-Base": "ssp585"
                }

#ssp_scn  = ["_remind_SSP1-PkBudg500", "_remind_SSP1-PkBudg1150", "_remind_SSP2-Base" , "_remind_SSP5-Base" ]  
ssp_scn  = ["_remind_SSP1-PkBudg500", "_remind_SSP2-Base" , "_remind_SSP5-Base" ]  
ssp_yr = [2030, 2040, 2050]
bg_ei_ = "ecoinvent_cutoff_3.9"
SSP_premise = "SSP1_base"    #input "ssp" if using premise_SSP database, otherwise, put whatever name e.g., "ei"


print("For dpLCIA, match dpCFs for each SSP[x], Year[t] by calling functions from dlcia_functions.py  " )
print("Mapping premise_remind_DB to SSPx defined here: ", ssp_remind_Pname_map)
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")



#initProject(PROJECT_NAME) # no need to iniiProject if already bw2data/bw.projects.set_current('iveo_v1')
def loadYearlyParams_multisheet(PARAMS_FILE = "data/p_scn_test.xlsx", s_name = ['V1A_g_truck', 'V1B_bat_LSB_perkWh'], hypothesis = 'point value', 
                                years = [2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060], 
                                SSP = "ssp119" ):
    """Load yearly parameters from excel, multiple worksheet, return a dict of {param_name -> [values]} """
    res_all = []
    for sheet in s_name: 
        df = pd.read_excel(PARAMS_FILE,  sheet_name = sheet, header=[0, 1])
        """ screen out only one SSPx """
        df = df[df['General info']["SSP"] == SSP]
        # note: BW parameters only linked with '_' allowed, need to rename the input excel , and I deleted <,> in the excel, only space
        df.index = df["General info"].name.apply(lambda x: x.replace(" ", "_") if type(x) == str else x)  #.replace(' ', '_').replace(',', '_')
        res = df[hypothesis][years].T.to_dict("list")
        if nan in res:
            del res[nan]
        res_all.append(res)

    ## final output is a dict, not list
    res_dict = {}
    for d in res_all:
        res_dict.update(d)
    return res_dict


def loadYearlyParams(PARAMS_FILE = "data/p_scn_test.xlsx", s_name = None, hypothesis = 'point value', years = [2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060], SSP = "ssp119" ):
    """Load yearly parameters from excel, return a dict of {param_name -> [values]} """
    df = pd.read_excel(PARAMS_FILE,  sheet_name = s_name, header=[0, 1])
    """ screen out only one SSPx """
    df = df[df['General info']["SSP"] == SSP]
    # note: BW parameters only linked with '_' allowed, need to rename the input excel , and I deleted <,> in the excel, only space
    df.index = df["General info"].name.apply(lambda x: x.replace(" ", "_") if type(x) == str else x)  #.replace(' ', '_').replace(',', '_')
    res = df[hypothesis][years].T.to_dict("list")
    if nan in res:
        del res[nan]
    return res



def loadYearlyParamswhsize(PARAMS_FILE = "data/p_scn_test.xlsx", s_name = 'LSB_battery', hypothesis = 'point value', years = [2020, 2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060], SSP = "ssp119", size = None):
    """Load yearly parameters from excel, return a dict of {param_name -> [values]} """
    df = pd.read_excel(PARAMS_FILE,  sheet_name = s_name, header=[0, 1])
    """ screen out only one SSPx and one size group """
    df = df[df['General info']["SSP"] == SSP]
    df = df[df['General info']["size"] == size]
    # note: BW parameters only linked with '_' allowed, need to rename the input excel , and I deleted <,> in the excel, only space
    df.index = df["General info"].name.apply(lambda x: x.replace(" ", "_") if type(x) == str else x)  #.replace(' ', '_').replace(',', '_')
    res = df[hypothesis][years].T.to_dict("list")
    if nan in res:
        del res[nan]
    return res

