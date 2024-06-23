import xarray as xr
import numpy as np
import pandas as pd
import os


print("before running TRL module, you should have empty vehicle dataset prepared already, with TRL assigned, this module is designed to get a TRL dataset with three TRL data variable assigned (TRL score (0-9) prepared already in TRL_1B, TRL_applicability (0 or 1) and TRL_availability (0 or 1)) after proper vehicle mapping, which will be used when building up a future vehicle model. Finally, you can visualize the TRL table")


# index_label for the raw DF (containing [0]type , [1]size,  [2]powertrain, [3]year, [4]SSP)
index_vtype = 0
index_size = 1
index_pw = 2
index_yr = 3
index_ssp = 4
index_mod = 5
index_tech = 6

""" 
in this module, we define several mapping functions, which can be used if more vehicle types/ sizes/ powertrains are to be investigated, 
for the garbage truck, we only need to map powertrain to its technology 
""" 

def load_v_data(data_folder, file_name, s_name, index_col ):
        full_path = os.path.join(data_folder, file_name)
        data = pd.read_excel(full_path, sheet_name =s_name , index_col=index_col)
        return data 


def get_avail_size_for_V (data, map_dict): 
    """
     0. define a mapping dict, for each type of vehicle, what are applicable powertrain
     returned: filtered data, a new pd.DF with reduced rows
    """
    print(f" example vehicle type <-> size mapping dict under folder [data/mapping]  ")
    assert all(isinstance(item, dict) for item in map_dict), "mapping_vehicle - to size should be a list of {dictionary}"
    
    indices_to_delete = []
    for index, row in data.iterrows():
        key, value = index[index_vtype], index[index_size]
        xd = {key: [value]}
        
        if any(xd == d for d in map_dict):
            pass  
        else:
            indices_to_delete.append(index)
    df = data.drop(indices_to_delete)
    return df


def get_avail_pw_for_V (data, map_dict): 
    """
     0. define a mapping dict, for each type of vehicle, what are applicable powertrain
     returned: filtered data, a new pd.DF with reduced rows
    """
    print(f" example vehicle type <-> powertrain mapping dict under folder [data/mapping]  ")
    assert all(isinstance(item, dict) for item in map_dict), "mapping_vehicle - to powertrain should be a list of {dictionary}"
    
    indices_to_delete = []
    for index, row in data.iterrows():
        key, value = index[index_vtype], index[index_pw]
        xd = {key: [value]}
        if any(xd == d for d in map_dict):
            pass  
        else:
            indices_to_delete.append(index)
    df = data.drop(indices_to_delete)
    return df

    
def get_avail_tech_for_V (data, map_dict): 
    """
     0. define a mapping dict, for each type of powertrain, what are applicable v_modules
     1. prepare the technology_applicablity score (depending on mapping powertrain and module, e.g., V1B_bat won't be applicable to ICEV
     2. only when applicablity == 1, then get availability score per TRL, if TRL > 6 for SSPx,Yt, then availability = 1
     returned: prepared a new pd.DF with two new col for technology_applicablity and technology_availability
    """
   
    assert all(isinstance(item, dict) for item in map_dict), "mapping_tech - to - module should be a list of {dictionary}"
    applicable = []
    availbale = []
    for index, row  in data.iterrows():
        key, value = index[2], index[5]  #w/o size is index[4],  w/h size is index[5] for module 
        xd = {key: [value]}
        if any(xd == d for d in map_dict):
            applicable.append(1)
        else:
            applicable.append(0)
    # assign availability score according to applicability and TRL 
    for a , t in zip(applicable, data["TRL"]): 
        if a == 0: 
            availbale.append(0)
        else: 
            if t>6: 
                availbale.append(1)
            else: 
                availbale.append(0)
    data["tech_appl"] = applicable
    data["tech_avail"] = availbale
    return data


class prepare_veh_tech:
    r""" 
     vd: the data_array prepared from either: 1) raw data w/o tech_appli or tech_avail score, or 2) data prepared with get_avail_tech_for_V()
     return: a prepared data array, vehicle_type,_module,_tech,for SSPx&v_year ?
     """
    def __init__(self, v_data = None):
        #v_data = pd.read_excel("data/xr_sample_4TRL.xlsx", sheet_name="test",  index_col=[0,1,2,3,4,5])
        vd = v_data.to_xarray()
        self.v_data = v_data
        self.vd = vd 
        self.ssp_scn =  vd["SSP"]
        self.vtype = vd["v_type"]
        self.powertrains =  vd["powertrain"]
        self.year = vd["v_year"]
        #self.module = vd["Module"]
        self.tech = vd["Tech"]
        #self.parameters =  

    def get_V_P_avail_comb_dict(self):
        """
        purpose: to truncate the prepared vd dataarray into a vd_dict (see below) according to available module/tech for each vehicle / powertrain
        no input needed, call self.vd once it's prepared by get_avail_tech_for_V()
        return: 
        a vd_dict in the format of: two layer keys: 
            out_key is vehicle type: vehicle[i], inner_key is available powertrain for each vehicle[i], powertrain[j]; 
            value for each inner_key is then a dataarray for all available module/tech. combination for v[i],p[j] for all SS[x]P/Year[t]
        to call a specific dataarray: using vd_dict["out_key:v_type"]["inner_key:powertrain_type"], e.g., vd_dict['g_truck']['BEV']
        """              
        vd = self.vd
        outkey_list = self.vtype.values.tolist()
        inkey_list = self.powertrains.values.tolist()
        
        vd_dict_xr, vd_dict_list, vd_dict_list2_right_dimupto4th  = {}, {}, {}   #vd_dict_list has the dim not combined, use vd_dict_list2_right_dimupto4th
        dd36 = self.v_data
        """
                attention: the final combination len is not: len(vvdd.v_year) * len(vvdd.SSP) * len(vvdd.Module) * len(vvdd.Tech) or
                len(vvdd.v_year) * len(vvdd.SSP) * len(vvdd.Tech)>, as modules are parellel 
                I use the raw self.v_data pandas then filter through pandas index
        """
        for outer_key in outkey_list:
            vd_dict_xr[outer_key] = {}     # Create an empty inner dictionary for each outer key
            vd_dict_list[outer_key] = {} 
            vd_dict_list2_right_dimupto4th[outer_key] = {} 
            for inner_key in inkey_list:
                vdd = vd.sel(v_type = outer_key, powertrain = inner_key)
                vvdd = vdd.where(vdd['tech_avail'] == 1, drop=True)
                comb = len(vvdd.v_year) * len(vvdd.SSP)  * len(vvdd.Tech)   #* len(vvdd.Module)
                print(f"available tech. for vehicle {outer_key} and powertrain {inner_key} for all years and SSPx before_mod/Tech combination: {comb}")
                vd_dict_xr[outer_key][inner_key] = vvdd
                
                """
                dd36 to get the right tech. combination, note comb_N is wrong
                """
                # format: query_string = f"(letter == '{condition[0]}' and number == {condition[1]})"
                pd_list , pd_list2_right_dimupto4th = [] , []
                for y in vvdd.v_year.values:
                    for sp in vvdd.SSP.values: 
                        for m in vvdd.Module.values: 
                            d63 = dd36.query(f" v_type == '{outer_key}' and (Module == '{m}' and v_year == {y} and powertrain == '{inner_key}' and SSP == '{sp}') ")
                            d36 = d63.index.tolist()
                            pd_list2_right_dimupto4th.append(d36[0]) 
                            for i in range(len(d36)):
                                pd_list.append(d36[i])
                                
                vd_dict_list[outer_key][inner_key] = pd_list
                vd_dict_list2_right_dimupto4th[outer_key][inner_key] = pd_list2_right_dimupto4th
    
        return vd_dict_xr, vd_dict_list, vd_dict_list2_right_dimupto4th
   

    
    def get_tech_TRL (self, vd_tech, to_display = True ):
        """
        purpose: to display the TRL table for a specific input tech, vd_tech
        input: 
            vd: the input data_array prepared
            vd_tech: the technology to assess 
            to_display: default true is to display the TRL table prepared for the selected vd_tech 
        return: one TRL table for a specific input tech
        """
        vd = self.vd 
        TRL_df = pd.DataFrame(index=vd["v_year"].values, columns=vd["SSP"].values)
        for index, row in TRL_df.iterrows():
            for column in TRL_df.columns:
                ar = vd["TRL"].sel(v_year = index, SSP = column, Tech = vd_tech).values
                # because there could be multiple coordinates not put in the sel() above, let's get all non-nan values, 
                # and check if they've been assigned with same TRL (you should have one TRL for each technology under SSP_x and Year_t)
                ar1 = ar[~np.isnan(ar)]
                #print(ar1)
                if not np.all(ar1 == ar1[0]): 
                    raise ValueError(f" The technologuy {vd_tech} does not have same TRL for v_year_ {index} and SSP_ {column}. ") 
                else: 
                    TRL_df.at[index, column] = int(ar1[0])
        # Highlight TRL >=7 
        def highlight_TRL(value):
            return 'background-color: lightgreen' if value > 6 else ''
            
        # Apply the styling function using df.style.apply
        if (to_display):
            styled_df = TRL_df.style.applymap(highlight_TRL)
            print(f"for the technology {vd_tech} the TRL table is: ")
            display(styled_df)
        return TRL_df

    
    def get_tech_TRL_dict(self, to_display = True):
        """
        purpose: build on get_tech_TRL(), to get a final dict, {key is all available tech, value is their TRL table}
        """
        vd = self.vd
        key_list = vd["Tech"].values.tolist()
        my_dict = {key: None for key in key_list}
        for key in my_dict:
            TRL_df = pd.DataFrame(index=vd["v_year"].values, columns=vd["SSP"].values)
            for index, row in TRL_df.iterrows():
                for column in TRL_df.columns:
                    ar = vd["TRL"].sel(v_year = index, SSP = column, Tech = key).values
                    ar1 = ar[~np.isnan(ar)]
                    #print(ar1)
                    if not np.all(ar1 == ar1[0]): 
                        raise ValueError(f" The technologuy {key} does not have same TRL for v_year_ {index} and SSP_ {column}. ") 
                    else: 
                        TRL_df.at[index, column] = int(ar1[0])
            # Assign the TRL table to the corresponding key
            my_dict[key] = TRL_df
            
        # Highlight TRL >=7 
        def highlight_TRL(value):
            return 'background-color: lightgreen' if value > 6 else ''
            
        for key in my_dict:
            df = my_dict[key]
            if (to_display):
                styled_df = df.style.applymap(highlight_TRL)
                print(f"for the technology {key} the TRL table is: ")
                display(styled_df)
                
        return my_dict
