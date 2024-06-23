import os
from lca_algebraic import initProject
import pandas as pd
import numpy as np
import lca_algebraic as agb
from lca_algebraic import *
from lca_algebraic.stats import * 
import bw2data, bw2io
import bw2calc
from premise import *
from premise.geomap import Geomap

from utils import *
from utils.getAct import *
from utils.utils import *


class switch_manuf_elec:
    print("not used in this GT study, utils.switchElec can update an UPRs with the switch elec parameter, allow users to switch manuf. electricity grid per their defined locations - market group or grid-specific")
    r""" 
    purpose is to update existing UPRs with the switch elec parameter, depending on user-defined locations for elec.
    to improve:
        change the location of changed activity
        allow different electricity voltage (currently only medium for manufacturing), linked to getAct.findSpecificElecwhSSP(voltage = 'medium')
        allow different electricity tech e.g., PV (currently only ecoinvent grid electricity)
    """
    def __init__(
        self,
        marketgroup = True,
        loc_to_search_marketgroup = ["CA", "ENTSO-E", "US", "CN", "GLO"],
        loc_to_search_nonmarketgroup = ["CA-QC",  "CA-BC", "CA-ON", "CA-AB",  #Canada doesnt play major role in battery manuf though
                                "JP", "KR", "SE", "FR", "DE",         #Japan, Korea, Sweden, France, Germany major manuf counties, w/o aggregrated elec market group 
                                "US-WECC", "US-SERC", "US-MRO", "US-NPCC",
                                "CN-CCG", "CN-CSG", "CN-SWG", 
                               ],
        PDB_NAME = "iveo_Parameterized_v1",
        bg_db = "ecoinvent 3.9.1",
        ssp_scn = ["SSP1_base" ], # "none_EI", "ssp126", "ssp245", "ssp370", "ssp434", "ssp460", "ssp534-over", "ssp585"
        vgroup = 'g_truck',
        vpowertrain = "BEV",
        manuf_y = 2025,              #2030,2040,2050,2060
        v_module = "1B_bat"      # 1A_veh, 
        #adding to parameter name to allow selecting different locations for different vehicle module e.g., battery vs. vehicle produced in different locations
    ):
        """Initialise"""
        self.PDB_NAME = PDB_NAME
        self.vgroup = vgroup
        self.vpowertrain = vpowertrain
        self.v_module = v_module
        self.manuf_y = manuf_y
        self.bg_db = bg_db
        self.ssp_scn = ssp_scn
        self.marketgroup = marketgroup
        self.loc_to_search_marketgroup = loc_to_search_marketgroup
        self.loc_to_search_nonmarketgroup = loc_to_search_nonmarketgroup

    def get_elec_mix(self):
        db_name = self.bg_db
        mg = self.marketgroup
        
        if self.marketgroup == True: 
            loc_to_search = self.loc_to_search_marketgroup
        elif self.marketgroup == False: 
            loc_to_search = self.loc_to_search_nonmarketgroup 
        
        manuf_elec_list = []
        for l in loc_to_search: 
            elec_i = getAct.findSpecificElecwhSSP(voltage = "medium", loc = l, marketgroup = mg, db_name = db_name)
            manuf_elec_list.append(elec_i)
        print("using background database '%s', for manufacturing year '%s', SSP scenario '%s', \
              following electricity_switch will be applied " % (self.bg_db, self.manuf_y, self.ssp_scn ) )
        for m in manuf_elec_list: print(m)
        return manuf_elec_list
        
    def get_elec_dict_foragb(self, manuf_elec_list):
        loc_value = [x.as_dict()["location"].replace(" ", "_").replace('-', '_') for x in manuf_elec_list]
        print(loc_value)
        elec_dict = {}
        for key, value in zip(loc_value, manuf_elec_list):
            elec_dict[key] = value
        print(elec_dict)
        return elec_dict

    def get_elec_mix_agb_Param(self, manuf_elec_list, elec_dict): 
        # Use a regular expression to match commas, spaces, square brackets 
        """
        once the elec_dict is prepared, get the electricity EnumParam and Switch parameters
        input: 
            manuf_elec_list: output from the above get_elec_mix() function according to user-defined locations (either market group or specific grid)
            elec_dict: output from the above get_elec_dict_foragb() function 
        return: 
            1.newEnumParam Parameter: Manuf_elec_mix_..., use Manuf_elec_mix_... .name if wanna run calculation for specified elec
            2.newSwitchAct Parameter: Manuf_elec_mix_..., this will be an input for switch_on_elec_forAct() function to switch-on the changing elec.
        """
        pattern = re.compile(r'[,\s\[\]\'\"\-\+]')
        loc_value = [x.as_dict()["location"].replace(" ", "_").replace('-', '_') for x in manuf_elec_list]
        X_Manufacturing_electricity_mix = agb.newEnumParam(
            name = 'Manuf_elec_mix_' + self.vgroup + self.vpowertrain + self.v_module + pattern.sub('_', str(self.ssp_scn)) + pattern.sub('_', str(self.manuf_y)),   
            values = loc_value,
            default = loc_value[0],
            #label = 'Electricty mix',
            description = "Switch on electricty mix used for X manufacture",
            group = self.vgroup)
            
        elec_switch_act = agb.newSwitchAct(dbname = self.PDB_NAME,
                                   name = 'Manuf_elec_mix_' + self.vgroup + self.vpowertrain + self.v_module + pattern.sub('_', str(self.ssp_scn)) + pattern.sub('_', str(self.manuf_y)), 
                                   paramDef = X_Manufacturing_electricity_mix, 
                                   acts_dict = elec_dict)
            
        print("use this EnumParam name for electricity:", X_Manufacturing_electricity_mix.name) 
        return X_Manufacturing_electricity_mix, elec_switch_act

        
    def switch_on_elec_copiedAct(self, bw2act_raw, elec_switch_act ): 
        """
        this is similar to switch_on_elec_forAct() function, except that the raw bw2act (e.g.,from a database) won't be changed directly but copied first, then changes made on the copied bw2act, make sure to link this copied bw2act to the product system.
        input: 
        bw2act_raw is the background UPR: #agb.findActivity('assembly operation, for lorry', db_name = 'ei391_SSP2_image_RCP26_2050')
        #elec_name: string input, the electricity to update, depending on what's used in the bw2act_raw, e.g, 'market group for electricity, low voltage'
        #elec_amount: either using the same electricity amount in the original bw2act_raw, to supply your own float number or parameter
        """
        bw2act_copy = copyActivity(
                self.PDB_NAME,
                bw2act_raw,
                bw2act_raw.as_dict()["name"] + '_to change elec')
        #print("this is the exchanges in the raw dataset")
        #for exc in bw2act_copy.exchanges(): print(exc)
        for exc in bw2act_raw.listExchanges():
            if 'electricity' in exc[0]:
                print(f"there is {exc} input in the original process {bw2act_raw}, will update it.")
                elec_name = exc[0]  #exc[1].as_dict()["name"] #first element of tuple 
                elec_amount = exc[2]                 #third element of tuple is the amount, use the same electricity input amount 

        bw2act_copy.updateExchanges({
                # update electricity : 'market group for electricity, low voltage' is the UPR used in original UPR, so overwirte it 
                # Exchange: 0.100805369 kilowatt hour 'market group for electricity, low voltage' (kilowatt hour, WEU, None) to 'assembly operation, for lorry, to change elec' (kilogram, RER, None)>
                elec_name: dict(amount= elec_amount , input=elec_switch_act),
            })
        #print("this is the exchanges in the new dataset with electricity updated")
        #for exc in bw2act_copy.exchanges(): print(exc)
        agb.printAct(bw2act_raw, bw2act_copy)
        return bw2act_copy
    
    
    def switch_on_elec_forAct(self, bw2act, elec_loc, elec_switch_act, elec_amount): 
        """
        once the elec_switch_act is prepared, add in the new electricity switch for a selected activity
        if adding new elec_exchanges is needed, calling getAct.findSpecificElecwhSSP(), do not enter name, enter voltage instead
        input: 
            bw2act: bw2activity to change the electricity 
            elec_loc: input the location, need to be one of the initially supplied in self.loc_to_search_nonmarketgroup or self.loc_to_search_marketgroup
            elec_swtich_act: elec switch parameter, 2nd output from the above get_elec_mix_agb_Param() function 
            elec_amount : could be a fixed value, or floatparameter itself
        return: 
            bw2activity with updated to-be-switch-on electricity 
        """
        if any('electricity' in exc[0] for exc in bw2act.listExchanges()):  #exc[0] is first element of tuple, string type name of the act
            print("here is exchanges before changing the electricity. ")
            agb.printAct(bw2act)
       
            if self.marketgroup == True: 
                if elec_loc not in self.loc_to_search_marketgroup:
                    raise ValueError(f" {elec_loc} is not a valid location. As you wanna market group, valid elements are: {self.loc_to_search_marketgroup}")
                            
                for x in bw2act.listExchanges(): 
                    if "market group for electricity" in str(x):  # aligned 
                        print("you wanna electricity market group and you have the market group for electricity in your UPR, aligned, updating it")
                        bw2act.updateExchanges({
                            'market group for electricity, medium voltage': dict(amount = elec_amount, input=elec_switch_act),
                        })
                        agb.printAct(bw2act)
                    elif "market for electricity" in str(x):     # not aligned 
                        print("you wanna electricity market group but you have the market for electricity in your UPR, \
                                delete the market for electricity, then add market group for electricity and switch on elec")
                        bw2act.deleteExchanges(name = 'market for electricity, medium voltage')
                        bw2act.addExchanges({
                              getAct.findSpecificElecwhSSP(voltage = 'medium', loc = elec_loc, marketgroup=True, db_name = self.bg_db): dict(amount = elec_amount, input=elec_switch_act),
                            })
                        agb.printAct(bw2act)
                    else:
                        pass 

            
            elif self.marketgroup == False:
                if elec_loc not in self.loc_to_search_nonmarketgroup:
                    raise ValueError(f" {elec_loc} is not a valid location. As you wanna specific grid elec not market group, valid elements are: {self.loc_to_search_nonmarketgroup}")
                for x in bw2act.listExchanges(): 
                    if "market for electricity" in str(x):  # aligned 
                        print("you wanna market for electricity and you have that in your UPR, aligned, updating it")
                        bw2act.updateExchanges({
                            'market for electricity, medium voltage': dict(amount = elec_amount, input=elec_switch_act),
                        })
                        agb.printAct(bw2act)
                    elif "market group for electricity" in str(x):     # not aligned 
                        print("you wanna market for electricity but you have the market group in your UPR, \
                                delete the market group for electricity, then add market for electricity and switch on elec")
                        bw2act.deleteExchanges(name = 'market group for electricity, medium voltage')
                        bw2act.addExchanges({
                              getAct.findSpecificElecwhSSP(voltage = 'medium', loc = elec_loc,  marketgroup=False, db_name = self.bg_db): dict(amount = elec_amount, input=elec_switch_act),
                            })
                        agb.printAct(bw2act)
            else: 
                print("failed to switch on electricity")

        else:
            print(f"there is no electricity input in the exchanges of {bw2act}. We'll add an electricity exchange. ")
            bw2act.addExchanges({
                    getAct.findSpecificElecwhSSP(voltage = 'medium', loc = elec_loc, marketgroup=self.marketgroup, db_name = self.bg_db): dict(amount = elec_amount, input=elec_switch_act),
            })
            agb.printAct(bw2act)
        return bw2act