import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

class rolling_Y_analysis:
    """ once fair is run and output RF / AGWP / GWP is saved to metric folder, fetch from each file for each gas: 
        if rolling_Y, only get IRF, as CRF / GWP will be calculated from this new IRF fetched for each MY[Y[0]]
    """
    print( " for SSP scn, only enter number e.g., 119, 585 without SSP \n for metric_toget, fetch only IRF if rolling_Y approach \n MYs should match MY defined for FaIR"  )
    def __init__(
        self,
        len_cf = 2,     # for rolling-Y, H_max = 1, (v2.0.1 notebook), len is 2, if H_max = 100, then len_cf = 101 
        ghggass = ["CO2",  "CH4", "N2O"],
        scn = ["119", "245", "585"],          # anal / plot three SSP
        metric_toget = 'IRF', 
        MYs = 2000 + np.arange(25,100+25,1),  # corresponding to v2.0.1 setup 
        cf_folder = None, 
        f_name_prefix = "agwp_dcf_gwp1_tstep1" # this is the excel prefix for rolling-Y approach, filename is hard-coded below  
        
    ):
        """Initialise"""
        self.len_cf = len_cf
        self.ghggass = ghggass
        self.scn = scn
        self.MYs = MYs
        self.metric_toget = metric_toget
        self.cf_folder = cf_folder
        self.f_name_prefix = f_name_prefix
        
        
    """ 3 functions to plot RE / alpha / concentration, can combine with more input_excel & if plotting label not hardcoded""" 
    def plot_gas_RE_whichgas (self, re_df, whichgas): 
        """ to plot re,  
        Input: 
        ----- 
        re_df:   stored under 'output/RollY' 
        whichgas: define which gas to plot CO2 , CH4 , N2O, then hardcode label
        """
        if whichgas == "CO2": 
            col_from_re = 'co2_erf_1ppm'
            tt = whichgas+'_RE:W/m2/ppm'
        elif whichgas == "CH4": 
            col_from_re = 'ch4_erf_1ppb'
            tt = whichgas+'_RE:W/m2/ppb'
        elif whichgas == "N2O": 
            col_from_re = 'n2o_erf_1ppb' 
            tt = whichgas+'_RE:W/m2/ppb'
        
        grouped = re_df.groupby("SSP")
        fig, ax = plt.subplots(figsize=(8,6))
        g_plot = lambda x:x.plot(x = "MY", y = col_from_re, ax=ax, label=x.name, 
                                 title = (tt + '  rolling-Y per MY: from ' +  str(self.MYs[0] ) + ' to ' + str(self.MYs[-1]) ))
        grouped.apply(g_plot)
        plt.xlim(self.MYs[0] - 5 ,  self.MYs[-1] + 5 )
        plt.show()
        
    
    def plot_gas_alpha_whichgas (self, alf_df, whichgas): 
        if whichgas == "CH4": 
            col_from_alf = 'alpha_ch4'
        elif whichgas == "N2O": 
            col_from_alf = 'alpha_n2o' 
        tt = whichgas+'_lifetime alpha '
        grouped = alf_df.groupby("SSP")
        fig, ax = plt.subplots(figsize=(8,6))
        g_plot = lambda x:x.plot(x = "MY", y = col_from_alf, ax=ax, label=x.name, 
                                 title = (tt + '  rolling-Y per MY: from ' +  str(self.MYs[0] ) + ' to ' + str(self.MYs[-1]) ))
        grouped.apply(g_plot)
        plt.xlim(self.MYs[0] - 5 ,  self.MYs[-1] + 5 )
        plt.show()
        
        
    def plot_gas_concentration_whichgas (self, concentration_df, whichgas): 
        if whichgas == "CO2": 
            col_from_df = 'co2_concentration'
            tt = whichgas+'_concentration: ppm'
        elif whichgas == "CH4": 
            col_from_df = 'ch4_concentration'
            tt = whichgas+'_concentration: ppb'
        elif whichgas == "N2O": 
            col_from_df = 'n2o_concentration' 
            tt = whichgas+'_concentration: ppb'
        
        grouped = concentration_df.groupby("SSP")
        fig, ax = plt.subplots(figsize=(8,6))
        g_plot = lambda x:x.plot(x = "MY", y = col_from_df, ax=ax, label=x.name, 
                                 title = (tt + '  rolling-Y per MY: from ' +  str(self.MYs[0] ) + ' to ' + str(self.MYs[-1]) ))
        grouped.apply(g_plot)
        plt.xlim(self.MYs[0] - 5 ,  self.MYs[-1] + 5 )
        plt.show()
        
        
        
    # get calculated metrics from ModuleB
    def get_gas_metric(self ):  
        """ 
        Returns: 
        ------- 
        cf_col_gas_allgas: list for three self.ghggass, cf_col_gas_allgas[i] has all the point value for metric read from excel   
        cf_col_gas: the last N2O output, used as input data shape for preparing DF in prep_CF_pds()
        """ 
        
        if self.metric_toget == 'IRF':
            metric_persheet =  'rf_pointvalue' 
        elif self.metric_toget == 'CRF':
            metric_persheet =  'agwp_pointvalue' 
        elif self.metric_toget == 'GWP': 
            metric_persheet =  'gwp_pointvalue' 
            
        cf_col_gas_allgas = []
        for ghggas in self.ghggass: 
            cf_col_gas = []
            for sp in self.scn:
                for MY in self.MYs:
                    filename = self.f_name_prefix  + ghggas + "_ssp" + sp + "_fair_start2000MY" + str(MY) + ".xlsx"
                    cf_file = os.path.join(self.cf_folder, filename)
                    print(f"read point value for {ghggas} SSP {sp}, model year {MY}, metric: {self.metric_toget}")
                    cf_df = pd.read_excel(cf_file, sheet_name = metric_persheet + ghggas+ "_" + str(MY) ) 
                    cf_pointvalue = cf_df[ghggas].values
                    cf_col_gas.append(cf_pointvalue)
                    
            cf_col_gas_allgas.append(cf_col_gas)
        
        return(cf_col_gas_allgas, cf_col_gas) 
    
    
    
    def prep_CF_pds (self , cf_col_gas_allgas, cf_col_gas ): 
        """ 
        Returns: 
        ------- 
        cf_df: characterization factor (for self.metric) for each gas
        """ 
        cf_col_gas_f = [item for sublist in cf_col_gas for item in sublist]
        # prepare other columns 
        len_cf = self.len_cf
        MYs = self.MYs  
        ssps = self.scn  
        sp_col, my_col = [], []

        for sp in ssps: 
            repeat_sp = [sp] * len_cf * len(MYs)
            sp_col.append(repeat_sp)
        sp_col_f = my_col_f = [item for sublist in sp_col for item in sublist]

        for my in MYs: 
            repeat_my = [my] * len_cf  
            my_col.append(repeat_my)

        my_col = my_col * len(ssps)
        my_col_f = [item for sublist in my_col for item in sublist] 
        year_col =  [np.arange(len_cf)] *len(ssps) * len(MYs)
        year_col_f =  [item for sublist in year_col for item in sublist]
        assert len(year_col_f) == len(sp_col_f) == len(my_col_f) == len(cf_col_gas_f)
        
        #co2_irf_all, ch4_irf_all, n2o_irf_all = cf_col_gas_allgas[0], cf_col_gas_allgas[1], cf_col_gas_allgas[2] 
        gas_all_f_list = []
        # gas_irf_all_f_list[i] is the ith gas from ghggass, flattened array for all SSPs, MYt in order 
        for i in range(len(cf_col_gas_allgas)): 
            gas_all_f = [item for sublist in cf_col_gas_allgas[i] for item in sublist]
            gas_all_f_list.append(gas_all_f)
        
        # col_names for the gas, to be combined with other col_names SSP MY Year
        col_gas_name = [] 
        for gsname in self.ghggass: 
            col_gas_name.append (gsname + "_" + self.metric_toget ) 
        
        col_names =  ["SSP", "ModelYear", "Year"] + col_gas_name
        print(col_names)
        cf_df = pd.DataFrame(zip(sp_col_f, my_col_f, year_col_f, gas_all_f_list[0], gas_all_f_list[1], gas_all_f_list[2]) , 
                     columns = col_names  ) 
        
        return cf_df
    
    
    
    def CF_pds_filterY0 (self ,cf_df): 
        """ 
        Input: 
        ------- 
        cf_df: prepared dataframe on the selected metric from prep_CF_pds()
        Returns: 
        ------- 
        cf_df_perMY: run this for IRF only, for rolling-Y, only get Y[0] output for IRF, for CRF/AGWP, need calculate from IRF 
        """ 
        if self.metric_toget == 'IRF':
            filter_y = 0 
        else:
            filter_y = 1 
        cf_df_perMY = cf_df[cf_df["Year"] == filter_y ]
        return cf_df_perMY
