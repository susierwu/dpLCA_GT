import bw2data, bw2io
import bw2calc
import numpy as np
import pandas as pd
import xarray as xr


class assign_dpIRF:
    """  
    to assign dIRF/CRF values to majorGHG per SSP[x]
    input: 
        cf_inputds: .nc data array prepared from FaIR output metrics
    """
    def __init__(self,
        cf_inputds, 
        ssp = '119',
        fairMY = 2030, 
        metric = 'IRF',
        TH = 101, #101 if GWP100, 501 if GWP500
        majorGHG_names = ['Carbon dioxide, non-fossil', 'Carbon dioxide, in air', 'Carbon dioxide, to soil or biomass stock', 
             'Carbon dioxide, from soil or biomass stock', 'Carbon dioxide, fossil', 'Carbon dioxide, non-fossil, resource correction', 
             'Methane, from soil or biomass stock', 'Methane, fossil', 'Methane, non-fossil' , 
             'Dinitrogen monoxide'],
        ):
        
        """Initialise"""
        self.cf_inputds = cf_inputds
        self.ssp = ssp
        self.TH = TH
        self.fairMY = fairMY
        self.metric = metric
        self.majorGHG_names = majorGHG_names


    def prep_empty_C (self ):
        """  
        to prepare empty_C datafrme, 
        output: 
            empty_C dataframe to be assigned with value, input "C" for assign_majorghg_dCC()
        """
        indexrow = []
        for i in range(0, self.TH):   #IRF needs to start from Year[0]
            x36 = self.metric + str(i)
            indexrow.append(x36)
        newcol = self.majorGHG_names
        empty_C = pd.DataFrame(0, index=indexrow,columns = newcol)
        return empty_C
           

    
    def assign_majorghg_dCC (self, C ):
        """  
        to prepare final_C dataframe with CF assigned, 
        input: 
            C: empty_C from prep_empty_C()
            the raw cf_inputds 
        output: 
            C dataframe with CF assigned
        """
 
        metric = self.metric
        cf_rawds = self.cf_inputds
        cf_touse = cf_rawds.sel(SSP=self.ssp, ModelYear=self.fairMY)
        co2_ ,ch4_,n2o_ = 'CO2_' + metric , 'CH4_' + metric ,  'N2O_' + metric 
        assert len(cf_touse[ch4_].values ) == len(C["Methane, fossil"]) , "the raw CF_ds has different numbers of CF as in the C_df to be assigned"

        # same as  GWP
        for col in C.columns:
            if "Carbon dioxide" in col: 
                if "Carbon dioxide, in air" in col:
                    C[col] = cf_touse[co2_].values * (-1)
                elif "Carbon dioxide, non-fossil, resource correction" in col:
                    C[col] = cf_touse[co2_].values * (-1)
                elif "Carbon dioxide, non-fossil" in col or 'Carbon dioxide, fossil' in col:   
                    C[col] = cf_touse[co2_].values
                elif "Carbon dioxide, to soil or biomass stock" in col :   
                    C[col] = cf_touse[co2_].values * (-1) 
                elif "Carbon dioxide, from soil or biomass stock" in col :   
                    C[col] = cf_touse[co2_].values  
            elif "Methane, fossil" in col or "Methane, non-fossil" in col or "Methane, from soil or biomass stock" in col:
                C[col] = cf_touse[ch4_].values 
            elif "Dinitrogen monoxide" in col: 
                C[col] = cf_touse[n2o_].values 
            else: 
                #print(f" for GHG {col}, assigned with static GWP100 already " ) 
                pass
        return C


    def prep_data_for_bw2method (self, C, mybio):
        """ 
        prepare final data for final bw2.method  
        note that one stressor mapped to multiple compartment, with the same CF, this holds true for climate change 
        input: 
            C: the fully_prepared C DF with all GHG's GWP100 assigned, from assign_majorghg_dCC()
            mybio: bw2 bio database 
        output: 
            data to be used in prep_final_dCC_bw2method(), it's a list with len of self.TH,  type(data[0][0]) is a turple for each BW2 flow
        """ 
        print(f"start preparing data to be assigned as new bw2data.methods, for SSP {self.ssp} and fair_MY{self.fairMY}")
        data = []  #for IRF:: 101 list of tuple for SSP[x]_MY[t], 101 new_method 
        for ind in C.index: # to loop through all years 
            cf_t = []
            for stressor in C.columns:
                all_efs = [i for i in mybio if i['name'] == stressor]
                """ you'll have multiple compartment for one emission """ 
                for ef in all_efs: 
                    stressor_code = ef['code']
                    #print(f" for {stressor}, for Year{ind}, you have code: {stressor_code} ")
                    """ we still relying on 'biosphere3', so name it as it is for LCIA """ 
                    cf_i = (('biosphere3', stressor_code), C.loc[ind, stressor])  
                    cf_t.append(cf_i)
            data.append(cf_t)
        return data

    
    def prep_final_dCC_bw2method (self, data):
        """ 
        write BW2method from data 
        input: 
            data from prep_data_for_bw2method()
        output: 
            null, write to BW2 project directly 
        """ 
        ssp, MY = str('SSP' + self.ssp) , str('MY' + str(self.fairMY))   
        m_name = 'IPCC 2021 - dpCFs' + ssp + '_' + MY + ' - year'
        print(f" start creating new method, with name: {m_name}")
        
        for i in range( 0, len(data) ):  
            name = (m_name + str(i) , 'climate change', 'dp' + self.metric + str(i) )
            #name = (m_name + str(i+1) , 'climate change', 'd' + self.metric + str(i+1) )
            new_method = bw2data.Method(name)
            new_method.register()
            if self.metric == 'IRF': 
                new_method.metadata["unit"] = 'W m-2 kg-1'
            elif self.metric == 'CRF': 
                new_method.metadata["unit"] = 'W m–2 yr kg–1'
            else: 
                print("metric is either IRF or CRF")
            new_method.write(data[i])
        print("finishing preparing new methods")






class assign_dpGWP:
    """  
    to assign dGWP values to majorGHG per SSP[x], while using the static GWP100 for other GHG per permise_gwp: 
    https://github.com/polca/premise_gwp/tree/main/premise_gwp/data -> lcia_gwp2021_100a_w_bio.xlsx
    """
    def __init__(self,
        premise_gwp100_inputdf,
        cf_inputds, 
        ssp = '119',
        fairMY = 2030, 
        metric = 'GWP', 
        TH = 101, #101 if GWP100, 501 if GWP500
        GWP100_only = True, # this first version, we don't have any dpLCIA for minorGHGs, so only using premise_gwp100 values
        majorGHG_names = ['Carbon dioxide, non-fossil', 'Carbon dioxide, in air', 'Carbon dioxide, to soil or biomass stock', 
             'Carbon dioxide, from soil or biomass stock', 'Carbon dioxide, fossil', 'Carbon dioxide, non-fossil, resource correction', 
             'Methane, from soil or biomass stock', 'Methane, fossil', 'Methane, non-fossil' , 
             'Dinitrogen monoxide'],
        ):
        
        """Initialise"""
        self.premise_gwp100_inputdf = premise_gwp100_inputdf
        self.cf_inputds = cf_inputds
        #assert isinstance(self.premise_gwp100_inputdf, pd.DataFrame), "premise_gwp100_df is not a Pandas DataFrame"
        self.ssp = ssp
        self.fairMY = fairMY
        self.metric = metric
        self.TH = TH
        self.majorGHG_names = majorGHG_names
        self.GWP100_only = GWP100_only
        
 
    def get_minorand_allGHG(self ): 
        prem_gwp_df = self.premise_gwp100_inputdf
        assert isinstance(prem_gwp_df, pd.DataFrame), "premise_gwp100_df is not a Pandas DataFrame"

        prem_gwp_df_allghg = list(prem_gwp_df['name'].values) 
        minorGHG_names  = list(filter(lambda x: x not in  self.majorGHG_names, prem_gwp_df_allghg))
        minorGHG_namesuniq = list(np.unique(minorGHG_names))
        allghg_names =  np.concatenate(( np.unique(minorGHG_names) , self.majorGHG_names), axis=0) 
        # check all GHG flows are complete 
        assert len(np.unique(prem_gwp_df['name'].values)) == len(allghg_names), "the GHG_flows are not complete"
        return  minorGHG_namesuniq, allghg_names

    def prep_empty_C (self , allghg_names ):
        indexrow = []
        for i in range(1 , self.TH):
            x36 = 'GWP' + str(i)
            indexrow.append(x36)
        newcol = allghg_names
        empty_C = pd.DataFrame(np.nan, index=indexrow,columns = newcol)
        return empty_C

    def assign_minorghg_to_C_TBD_whdpLCIA_minorGHG (self, C, minorGHG_namesuniq): 
        """ 
        input: C: the empty_C DF prepared from self.prep_empty_C(), 
               minorGHG_namesuniq: the 1st returned value from get_minorand_allGHG(), 
               prem_gwp_df: the raw premise_GWP DF 
        """ 
        prem_gwp_df = self.premise_gwp100_inputdf
        for gas in C.columns: 
            if gas in minorGHG_namesuniq:
                gwp100_gas = prem_gwp_df[prem_gwp_df["name"] == gas]['amount'].values
                #print(f' for GHG {gas} using static GWP100 from prem_gwp_df, GWP100 is  {gwp100_gas}' )
                if all(x == gwp100_gas[0] for x in gwp100_gas):
                    gwp100_gas_value = gwp100_gas[0]
                else:
                    print("diff GWP for diff compartment")
                C[gas] = [gwp100_gas_value] * len(C[gas].values) 
            else: 
                pass 
        return C


    def assign_minorghg_to_C_GWP100 (self, C, minorGHG_namesuniq): 
        """
        get static GWP100 from premise_gwp
        """
        prem_gwp_df = self.premise_gwp100_inputdf
        if self.GWP100_only == True:
            for ind in C.index: 
                if ind == "GWP100":
                    for gas in C.columns: 
                        if gas in minorGHG_namesuniq:
                            gwp100_gas = prem_gwp_df[prem_gwp_df["name"] == gas]['amount'].values
                            #print(f' for GHG {gas} using static GWP100 from prem_gwp_df, GWP100 is  {gwp100_gas}' )
                            if all(x == gwp100_gas[0] for x in gwp100_gas):
                                gwp100_gas_value = gwp100_gas[0]
                            else:
                                print("diff GWP for diff compartment")
                            C.loc[ind, gas] = gwp100_gas_value 
                        else:
                            pass 
        else:
            print("dpLCIA for minor GHGs are to be prepared")
        return C


    def assign_majorghg_dCC (self, C, metric = 'GWP'):
        """ 
        input: C: the partially filled C from  assign_minorghg_to_C_GWP100()
        """ 
        cf_rawds = self.cf_inputds
        cf_touse = cf_rawds.sel(SSP=self.ssp, ModelYear=self.fairMY)
        co2_ ,ch4_,n2o_ = 'CO2_' + metric , 'CH4_' + metric ,  'N2O_' + metric 
        assert len(cf_touse[ch4_].values ) == len(C["Methane, fossil"]) , "the raw CF_ds has different numbers of CF as in the C_df to be assigned"
    
        for col in C.columns:
            if "Carbon dioxide" in col: 
                if "Carbon dioxide, in air" in col:
                    C[col] = cf_touse[co2_].values * (-1)
                elif "Carbon dioxide, non-fossil, resource correction" in col:
                    C[col] = cf_touse[co2_].values * (-1)
                elif "Carbon dioxide, non-fossil" in col or 'Carbon dioxide, fossil' in col:   
                    C[col] = cf_touse[co2_].values
                elif "Carbon dioxide, to soil or biomass stock" in col :   
                    C[col] = cf_touse[co2_].values * (-1) 
                elif "Carbon dioxide, from soil or biomass stock" in col :   
                    C[col] = cf_touse[co2_].values 
    
            elif "Methane, fossil" in col or "Methane, non-fossil" in col or "Methane, from soil or biomass stock" in col:
                C[col] = cf_touse[ch4_].values 
    
            elif "Dinitrogen monoxide" in col: 
                C[col] = cf_touse[n2o_].values 
            
            else: 
                #print(f" for GHG {col}, assigned with static GWP100 already " ) 
                pass
        return C

    
    def prep_data_for_bw2method (self, C, mybio):
        """ 
        prepare final data for final bw2.method new_method2
        note that one stressor mapped to multiple compartment, and we'll use same CF
        input: C: the fully_prepared C DF with all GHG's GWP100 assigned, from assign_majorghg_dCC()
               mybio: bw2 bio database 
        """ 
        print(f"start preparing data to be assigned as new bw2data.methods, for SSP {self.ssp} and fair_MY{self.fairMY}")
        data = []  #100 list of tuple for SSP[x]_MY[t], 100 new_method
        
        if self.GWP100_only != True: 
            print("dpLCIA for minor GHGs are to be prepared")
        else: 
            for ind in C.index: # to loop through all years  
                if ind == 'GWP100': 
                    cf_t = []
                    for stressor in C.columns:
                        all_efs = [i for i in mybio if i['name'] == stressor]
                        """ you'll have multiple compartment for one emission """ 
                        for ef in all_efs: 
                            stressor_code = ef['code']
                            #print(f" for {stressor}, for Year{ind}, you have code: {stressor_code} ")
                            """ we still relying on 'biosphere3', so name it as it is for LCIA """ 
                            cf_i = (('biosphere3', stressor_code), C.loc[ind, stressor])  
                            cf_t.append(cf_i)
                    data.append(cf_t)
        return data

    def prep_final_dCC_bw2method (self, data):
        ssp, MY = str('SSP' + self.ssp) , str('MY' + str(self.fairMY))   
        m_name = 'IPCC 2021 - dpCFs' + ssp + '_' + MY + ' - year'
        print(f" start creating new method, with name: {m_name}")
        for i in range( 0, len(data) ):  # pGWP only has one if self.GWP100_only == True, prepare for pGWP100
            if self.GWP100_only != True:
                name = (m_name + str(i) , 'climate change', 'dpGWP' + str(i) )
            elif self.GWP100_only == True:
                name = (m_name + "100" , 'climate change', 'pGWP100' )
            new_method = bw2data.Method(name)
            new_method.register()
            new_method.metadata["unit"] = 'kg CO2-Eq'
            new_method.write(data[i])
        print("finishing preparing new methods")
        
