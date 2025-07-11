import os
import xarray as xr
import sys
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
import ipywidgets as widgets
from IPython.display import display


def user_define_bat_tech(): 
    """ 
    Function to display two dropdowns and handle user selection.
    Return:
    -------
    bat1 and bat2  as global variable
    """ 
    # use label1/2  instead of description within dropdown to make sure not cutoff
    label1 = widgets.Label(value='Select initial battery technology:')
    label2 = widgets.Label(value='Select battery replacement technology:')
    
    # Define the first dropdown list
    dropdown1 = widgets.Dropdown(
        options=['LFP', 'LSB', 'LTO', 'NCA', 'NMC622'],
        value='NMC622',
        #description='Select initial bat tech:', #use label1/2 above instead otherwise cutoff
        layout=widgets.Layout(width='236px' )
    )
    
    # Define the second dropdown list
    dropdown2 = widgets.Dropdown(
        options=['LFP', 'LSB', 'LTO', 'NCA', 'NMC622'],
        value='LFP',
        #description='Select bat replacement tech:', #use label1/2 above instead otherwise cutoff
        layout=widgets.Layout(width='236px' )
    )
    # Function to handle the selection for the first dropdown
    def on_change1(change):
        global bat1
        bat1 = change['new']
        print(f'Select initial battery tech: {bat1}')
    
    # Function to handle the selection for the second dropdown
    def on_change2(change):
        global bat2
        bat2 = change['new']
        print(f'Select battery replacement tech: {bat2}')
    
    # Attach the functions to the dropdowns
    dropdown1.observe(on_change1, names='value')
    dropdown2.observe(on_change2, names='value')
    # Display the dropdowns
    display(widgets.VBox([label1, dropdown1, label2, dropdown2]))

    return dropdown1, dropdown2 



#at1 = dropdown1.value
#bat2 = dropdown2.value

def get_bev_df_per_input_tech(bev, bat1, bat2): 
    for i in range(len(bev)):  
        if (bev[i][0].get('V1AB_init')) ==  bat1 and (bev[i][0].get('V1B')) ==  bat2:
            print(f'index_number {i} from BEV dict for selected bat_tech') 
            bev_ind = i
    mybev = bev[bev_ind][1]
    return mybev


def get_icev_map_bev_drop_nan(mybev, icev):
    mybev2 = mybev.dropna()
    icev2 = icev[(icev['SSP'].isin(list(np.unique(mybev2["SSP"])))) & (icev['v_year'].isin(list(np.unique(mybev2["v_year"])))) ] 
    #print(len(icev2), len(mybev2))
    assert len(icev2) == len(mybev2), "BEV and ICEV DataFrames has different length."

    # adding a column for powertrain
    mybev2['Powertrain'] = 'BEV'
    icev2['Powertrain'] = 'ICEV'
    # rename icev dataframe V1A to V1AB to make it a common col_name for plotting
    icev2 = icev2.rename(columns={'V1A': 'V1AB' }) 
    
    return  mybev2, icev2


def prepare_concat_df_toplot(mybev2, icev2):
    df = pd.concat([mybev2, icev2], ignore_index=True)
    df = df[df['IC'] == 'pGWP100']
    
    df_list = []
    # for final plotting, always nine plots regardless of whether or not has avaiable BEV 
    for sp in ["ssp119","ssp245","ssp585"]: 
        for vy in [2030,2040,2050]: 
            dd = df [ (df['SSP'] == sp)  & (df['v_year'] == vy)   ]
            df_list.append(dd)
    
    return df_list 
