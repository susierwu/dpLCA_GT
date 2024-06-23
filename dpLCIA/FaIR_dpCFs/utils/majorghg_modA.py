import os 
import numpy as np
import pandas as pd 

class majorghg_get_f:
    """ in order to calculate climate metrics, this first module A runs 3 major functions:
    1. get majorghg fair output parameters 
    2. get majorghg ERF according to meinshausen2020(): src/fair/forcing/ghg.py
        - computer hardware (1.4 GHz Quad-Core Intel Core i5|8 GB 2133 MHz) reports memory erros running fair.meinshausen2020()
    3. calculate majorghg ERF per unit emission pertubation ([W m-2 ppm-1] CO2, [W m-2 ppb-1] CH4 and N2O), 
    which are used as input parameters in Module B calculating metrics 
    ########## there must be more efficient way to handling these functions through xarray DataArray ##########
    """
    
    def __init__(
        self,
        f,
        #majorghg = ["CO2",  "CH4", "N2O"],
        scn,
        H_max=100,
        ts_per_year = 1,
        scenarios = ["ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp460", "ssp534-over", "ssp585"],
        year_index = 269,  #1750 + 269 = year2019
        #y_range=[2020,2030,2040,2050,2060],
        #M_ATMOS = 5.1352E18,
        #M_AIR = 28.97E-3,
        #M_CO2 = 44.01E-3,
        #M_C = 12.0E-3,
        #M_CH4 = 16.043E-3,  
    ):
        """Initialise"""
        self.f = f
        self.scn = scn
        self.H_max = H_max
        self.year_index = year_index 
        self.ts_per_year = ts_per_year
        self.H = np.linspace(0, self.H_max, self.H_max * self.ts_per_year + 1) #parameter for co2_analytical()
        self.co2_f = self.call_f_from_fair_gas()[0]
        self.ch4_f = self.call_f_from_fair_gas()[1]
        self.n2o_f = self.call_f_from_fair_gas()[2] 
        self.n2o_c = self.call_f_from_fair_gas()[2]['N2O_C']


    # step A.1: call from fair, get parameter needed for climate metric calculation 
    def call_f_from_fair_gas(self):
        """
        call external fair class, its output parameters are input parameters for metric calculation 
        gas_c vs. gas_c_plus1 is used to calculate pulse emission for final metrics
        Inputs:
        -------
        Returns: 
        -------
        dict for {co2_f}, {ch4_f}, {n2o_f}:
            _C: concentration, _C_plus1: concentration + 1ppm/ppb;  _RF: reference baseline_concentration  _FS: forcing_scale 
        """
        # for each gas, prepare empty dict, with keys, value to be added 
        co2_key = ['CO2_C', 'CO2_C_plus1', 'CO2_RF', "CO2_FS", "N2O_C"]
        ch4_key = ['CH4_C', 'CH4_C_plus1', 'CH4_RF', "CH4_FS", "N2O_C", "CH4_lifetime"]
        n2o_key = ['N2O_C', 'N2O_C_plus1', 'N2O_RF', "N2O_FS", "N2O_lifetime"]
        co2_f = {key: None for key in co2_key} 
        ch4_f = {key: None for key in ch4_key} 
        n2o_f = {key: None for key in n2o_key} 

        for gas in ["CO2",  "CH4", "N2O"]:
            concentration = self.f.concentration.loc[dict(scenario=self.scn)]
            gas_c = concentration.loc[dict(specie=gas)].values
            gas_c_plus1 =  concentration.loc[dict(specie=gas)].values + 1 
            gas_rf = np.ones_like(gas_c) * self.f.species_configs["baseline_concentration"].loc[dict(specie=gas)].values
            gas_fs = np.ones_like(gas_c_plus1) * self.f.species_configs["forcing_scale"].loc[dict(specie=gas)].values
            #n2o needed for CO2
            n2o_c =  concentration.loc[dict(specie="N2O")].values              
         
            # assign values to dict
            if gas == "CO2": 
                co2_f['CO2_C'], co2_f['CO2_C_plus1'], co2_f['CO2_RF'], co2_f['CO2_FS'], co2_f['N2O_C'] = gas_c , gas_c_plus1, gas_rf, gas_fs, n2o_c
                
            elif gas == "CH4":
                    """ 
                    # geting pert lifetime per scenario: 
                    # f.species_configs["unperturbed_lifetime"].loc[dict(specie = "CH4")]  (1001, 4)
                    # gasbox: 4 ?  
                    # unpertubed_lifetime for CH4 is 10.8537568 so this point value is used here 
                    """
                    alpha_s = self.f.alpha_lifetime.loc[dict(scenario=self.scn, specie = "CH4")][self.year_index].values 
                    ch4_unp =  10.8537568   
                    alpha_ch4 = ch4_unp  * alpha_s
                    ch4_f['CH4_lifetime'] = alpha_ch4
                    ch4_f['CH4_C'], ch4_f['CH4_C_plus1'], ch4_f['CH4_RF'], ch4_f['CH4_FS'], ch4_f['N2O_C'] =  gas_c , gas_c_plus1, gas_rf, gas_fs, n2o_c 

            elif gas == "N2O": 
                    alpha_s = self.f.alpha_lifetime.loc[dict(scenario=self.scn, specie = "N2O")][self.year_index].values
                    n2o_unp = 109     
                    alpha_n2o = n2o_unp  * alpha_s 
                    n2o_f['N2O_lifetime'] = alpha_n2o
                    n2o_f['N2O_C'], n2o_f['N2O_C_plus1'], n2o_f['N2O_RF'], n2o_f['N2O_FS'] = gas_c , gas_c_plus1, gas_rf, gas_fs 

        return co2_f, ch4_f, n2o_f 
    

    # step A.2: get co2 / ch4 / n2o ERF 
    # get_co2/ch4/n2o_meinshausen2020: using raw function from FAIR: src/fair/forcing/ghg.py
    # "https://github.com/OMS-NetZero/FAIR/blob/master/src/fair/forcing/ghg.py"
    def get_co2_meinshausen2020 (self, co2_c_orcplus, #n2o_c, 
                            a1=-2.4785e-07,b1=0.00075906,c1=-0.0021492,d1=5.2488,
                            a2=-0.00034197,b2=0.00025455,c2=-0.00024357,d2=0.12173,
                            a3=-8.9603e-05,b3=-0.00012462,d3=0.045194): 
        """
        Inputs: 
        -------
        run from above call_f_from_fair_gas() function
        gas_c_orcplus: f.output, gas_concentration or gas_concentrationplus1ppm  
        gas_rf: f.specie_config, reference_concentration array   
        n2o_c: f.output, n2o concentration  
        Returns:
        -------
        erf_co2: used in 
        """
                 
        co2 = co2_c_orcplus
        co2_base = self.co2_f["CO2_RF"] 
        n2o = self.n2o_c
        co2_fs = self.co2_f["CO2_FS"] 
        ca_max = co2_base - b1 / (2 * a1)
    
        #CO2
        where_central = np.asarray((co2_base < co2) & (co2 <= ca_max)).nonzero()  #get true 
        where_low = np.asarray((co2 <= co2_base)).nonzero()
        where_high = np.asarray((co2 > ca_max)).nonzero()

        alpha_p = np.ones_like(co2) * np.nan
        alpha_p[where_central] = (
                d1
                + a1 * (co2[where_central] - co2_base[where_central]) ** 2
                + b1 * (co2[where_central] - co2_base[where_central])
            )
        alpha_p[where_low] = d1
        alpha_p[where_high] = d1 - b1**2 / (4 * a1)

        alpha_n2o = c1 * np.sqrt(n2o)

        erf_co2  = (
                (alpha_p + alpha_n2o) * np.log(co2 / co2_base) * (co2_fs)
            )
        return erf_co2


    def get_ch4_meinshausen2020 (self, ch4_c_orcplus,   
                                 #n2o_c, 
                            a1=-2.4785e-07,b1=0.00075906,c1=-0.0021492,d1=5.2488,
                            a2=-0.00034197,b2=0.00025455,c2=-0.00024357,d2=0.12173,
                            a3=-8.9603e-05,b3=-0.00012462,d3=0.045194): 
        """
        Inputs: run from call_f_from_fair_gas() function
        """
    
        ch4 = ch4_c_orcplus
        ch4_base = self.ch4_f["CH4_RF"] 
        n2o = self.ch4_f["N2O_C"] # == self.n2o_f["N2O_C"] == self.n2o_f["N2O_C"] 
        ch4_fs = self.ch4_f["CH4_FS"] 

        # CH4
        erf_ch4  = (
            (a3 * np.sqrt(ch4) + b3 * np.sqrt(n2o) + d3)
            * (np.sqrt(ch4) - np.sqrt(ch4_base))
        ) *  ch4_fs
    
        return erf_ch4

        
    def get_n2o_meinshausen2020 (self, n2o_c_orcplus,   
                                 #n2o_c, 
                            a1=-2.4785e-07,b1=0.00075906,c1=-0.0021492,d1=5.2488,
                            a2=-0.00034197,b2=0.00025455,c2=-0.00024357,d2=0.12173,
                            a3=-8.9603e-05,b3=-0.00012462,d3=0.045194): 
        """
        Inputs: run from call_f_from_fair_gas() function
        """
    
        n2o = n2o_c_orcplus
        n2o_base = self.n2o_f["N2O_RF"] 
        ch4 = self.ch4_f["CH4_C"]
        co2 = self.co2_f["CO2_C"]
        n2o_fs = self.n2o_f["N2O_FS"] 

        # N2O
        erf_n2o = (
            (a2 * np.sqrt(co2) + b2 * np.sqrt(n2o) + c2 * np.sqrt(ch4) + d2)
            * (np.sqrt(n2o) - np.sqrt(n2o_base))
        ) * n2o_fs

        return erf_n2o
    

    # Step A.3 get erf per 1 unit pertubation [W m-2 ppm/b-1]
    def get_co2_1ppm_erf(self): 
        """
        erf_co2: from above functions get_co2_meinshausen2020(), 
        running get_co2_meinshausen2020() twice, for C then C_plus, to get one unit pluse emission 
        Returns:
        -------
        co2_erf_diff: used for Module_B, to replace the static radiative efficiency value
        """
        erf1 = self.get_co2_meinshausen2020(self.co2_f["CO2_C"] ) 
        erf2 = self.get_co2_meinshausen2020(self.co2_f["CO2_C_plus1"])  
        co2_erf_diff = erf2 - erf1
        return co2_erf_diff
 
        
    def get_ch4_1ppb_erf(self): 
        """
        """
        erf1 = self.get_ch4_meinshausen2020(self.ch4_f["CH4_C"] ) 
        erf2 = self.get_ch4_meinshausen2020(self.ch4_f["CH4_C_plus1"])  
        ch4_erf_diff = erf2 - erf1
        return ch4_erf_diff
    
    def get_n2o_1ppb_erf(self): 
        """
        """
        erf1 = self.get_n2o_meinshausen2020(self.n2o_f["N2O_C"] ) 
        erf2 = self.get_n2o_meinshausen2020(self.n2o_f["N2O_C_plus1"])  
        n2o_erf_diff = erf2 - erf1
        return n2o_erf_diff

        