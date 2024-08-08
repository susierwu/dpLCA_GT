import os 
import numpy as np
import pandas as pd
import re
from datetime import datetime


class minorghg_get_f_anal_metric:
    """ 
    for each minorghg (user_input e.g., 'CFC-11', first get fair output parameters, then, calculate metrics, combined in the same class for minorghg to re-use self.f output paramters
    For minorGHG, two versions of metric calculation: 
    1. fixed RE and alpha from FaIR (set metric_wh1001 = False), so final RF/AGWP/GWP only has point values, no 1001 ensemble as majorGHG
    2. running get_minorghg_1ppb_erf() below to get 1001 RE, then similar metric calculation as majorGHG with uncertainty
    all minorGHG RE / lifetime / mass.../ mapping GHG names in GHG spreadsheet and RCMIP: https://github.com/chrisroadmap/ar6/blob/main/src/ar6/constants/gases.py
https://github.com/chrisroadmap/ar6/blob/main/data_output/7sm/metrics_supplement_cleaned.csv 
    as of 2020-08-08, missing 13 minorGHG (/Users/susierwu/dpLCA_main/dpLCIA/LCIA/premise_gwp/lcia_gwp2021_100a_w_bio_map_minorGHG_AR6_fair.xlsx) can't be mapped to FaIR: CH2=CF2, n-C4H10, C2H6, CH2ClCH2Cl, CHClFCF4, CHCl2CF3, CH2ClF, CClF3, CHCl2F, CH3CH2Cl, C3H8, CCl2=CCl2, CHCl=CCl2
    """
    def __init__(
        self,
        f,
        minorghg,
        scn,
        H_max=100,
        ts_per_year = 1,
        metric_wh1001 = True,  # added for minorGHG, majorGHG only default with 1001
        fair_start_y = 1750,  # we used 2000 as the fair start year
        year_index = 269,  #1750 + 269 = year2019 for start year 1750
        # parameter needed for metric calc
        q = np.array([0.443767728883447, 0.313998206372015]), 
        d = np.array([3.424102092311, 285.003477841911]), 
        M_ATMOS = 5.1352E18,
        M_AIR = 28.97E-3
        
    ):
        """Initialise"""
        self.f = f
        self.minorghg = minorghg
        self.scn = scn
        self.H_max = H_max
        self.fair_start_y = fair_start_y
        self.year_index = year_index 
        self.ts_per_year = ts_per_year
        self.metric_wh1001 = metric_wh1001
        self.H = np.linspace(0, self.H_max, self.H_max * self.ts_per_year + 1)  
        self.gas_f = self.call_f_from_fair_minorgas()
        self.q = q
        self.d = d
        self.M_ATMOS = M_ATMOS
        self.M_AIR = M_AIR
        
    # step A.1: call from fair, get parameter needed for climate metric calculation 
    def call_f_from_fair_minorgas(self):
        """
        call external fair class, its output parameters are input parameters for metric calculation, gas_c vs. gas_c_plus1 is used to calculate pulse emission, all minorGHGs in units of ppt.
        Inputs:
        -------
        Returns: 
        -------
        gas_f: fair dict for one minorGHG: {gas_f}:
        """
        # for each gas, prepare empty dict, with keys, value to be added 
        gas_key = ['gas_C', 'gas_C_plus1', 'gas_RF', "gas_FS", "gas_lifetime"]
        gas_f = {key: None for key in gas_key}  
            
        concentration = self.f.concentration.loc[dict(scenario=self.scn)]
        gas_c = concentration.loc[dict(specie=self.minorghg)].values
        gas_c_plus1 =  concentration.loc[dict(specie=self.minorghg)].values + 1 
        gas_rf = np.ones_like(gas_c) * self.f.species_configs["baseline_concentration"].loc[dict(specie=self.minorghg)].values
        gas_fs = np.ones_like(gas_c_plus1) * self.f.species_configs["forcing_scale"].loc[dict(specie=self.minorghg)].values
        #unper_lifetime and scaling factor with shape of (1001,) for each gas, it's the same value but we keep the 1001   
        gas_unp = self.f.species_configs['unperturbed_lifetime'].loc[dict(specie = self.minorghg, gasbox=0)].data
        alpha_s = self.f.alpha_lifetime.loc[dict(scenario=self.scn, specie = self.minorghg)][self.year_index].values 
        alpha_gas = gas_unp  * alpha_s
        gas_f['gas_lifetime'] = alpha_gas
        gas_f['gas_C'], gas_f['gas_C_plus1'] =  gas_c , gas_c_plus1
        gas_f['gas_RF'], gas_f['gas_FS'] =  gas_rf, gas_fs
        return gas_f 
    
    
    # step A.2: get minorghg ERF: using raw function from FAIR: src/fair/forcing/ghg.py
    # "https://github.com/OMS-NetZero/FAIR/blob/master/src/fair/forcing/ghg.py"
    # unlike majorGHG, here RE is an input parameter, read from f:: gas_re_fixed 
    def get_minorghg_meinshausen2020 (self, minorghg_c_orcplus): 
        """
        to get get minorghg ERF: using raw function from FAIR:
        "https://github.com/OMS-NetZero/FAIR/blob/master/src/fair/forcing/ghg.py" 
        unlike majorGHG, here RE is an input parameter
        Inputs:
        -------
        minorghg_c_orcplus: run from call_f_from_fair_minorgas() function
        Returns: 
        -------
        1st returned value: ERF for either C or Cplus1
        2nd returned value: the fixed RE get from f.species_configs[]
        """
        gas = minorghg_c_orcplus
        gas_base = self.gas_f["gas_RF"] 
        minorghg_fs = self.gas_f["gas_FS"] 
        
        gas_re_fixed = self.f.species_configs["greenhouse_gas_radiative_efficiency"].loc[dict(specie=self.minorghg)].data
        
        erf_minorghg = (gas - gas_base ) * gas_re_fixed  * minorghg_fs # * 0.001
        
        # same eq. from the raw function below erf_out ... : linear for other gases
        '''
        erf_out[..., minor_greenhouse_gas_indices] = (
            (concentration[..., minor_greenhouse_gas_indices]
                - reference_concentration[..., minor_greenhouse_gas_indices]
            ) * radiative_efficiency[..., minor_greenhouse_gas_indices]
            * 0.001  # unit handling
        ) * (forcing_scaling[..., minor_greenhouse_gas_indices])
        '''
        return erf_minorghg, gas_re_fixed
    
    
    # Step A.3 get erf per 1 unit pertubation [W m-2 ppm/b-1] 
    def get_minorghg_1ppb_erf(self): 
        """
        If point-value metric calculation (minorghg_analytical_fixedRE), using fixed RE from FaIR directly, this function won't be used, it serves two purposes:  
        1. cross-validation purposes (gas_erf_diff calculated here vs. fixed RE)
        2. check (in the notebook) if RE re-calculated here is same for all 1001 ensemble 
        for minorGHG, RE is obtained as a fixed value from FaIR, e.g., gas_re_fixed from get_minorghg_meinshausen2020. this function re-calc ERF (same method as for major_GHG), running this equation get very slightly diff RE from FaIR fixed RE, it doens't change perSSP/MY, e.g., for 'HFC-236fa', RE from FaIR is 0.25069, here reversed calculation, median is 0.248691, same for all SSP[x] / Y[t]
        NOTE: did not * 0.001 in get_minorghg_meinshausen2020(), so output unit is ppb not ppt
        Inputs:
        -------        
        Returns:
        -------
        gas_erf_diff: shape of (551, 1001) e.g., 551 as starting 1750, ending 2300

        """
        erf1 = self.get_minorghg_meinshausen2020(self.gas_f["gas_C"])[0]  #first returned value erf_minorghg
        erf2 = self.get_minorghg_meinshausen2020(self.gas_f["gas_C_plus1"])[0]
        gas_erf_diff = erf2 - erf1
        return gas_erf_diff 
    
    
    # Step B (ModuleB): calculate metric, combined in the same class for minorghg to re-use self.f output paramters

    def minorghg_analytical (self, gas_erf_diff, halogen_ra=0 ):
        """
        In metric calculation, using fixed RE from FaIR directly (self.f.species_configs["greenhouse_gas_radiative_efficiency"])
        Impluse Response Function same as: https://github.com/chrisroadmap/ar6/blob/main/src/ar6/metrics/halogen_generic.py
        NOTE: for minorGHG, as RE and alpha are fixed, so final RF/AGWP only has point values, no 1001 ensemble as majorGHG
        Input:
        --------
        gas_erf_diff: RE with 1001 ensemble run from get_minorghg_1ppb_erf(), only used if metric_wh1001 = True
        metric_wh1001: if users wanna run fixed metrics calculation (by fixed RE) or with 1001 (using gas_erf_diff)
        Returns:
        --------
        rf, agwp, agtp, iagtp as in C.Smith (halogen_generic.py)
        """
        
        y = self.year_index 
        q = self.q
        d = self.d    
        H = self.H 
        # molecular_weight / 1000  to convert to kg/mol needed in metric calc.  
        mass = self.f.species_configs["molecular_weight"].loc[dict(specie=self.minorghg)].values / 1000 
        alpha = self.gas_f['gas_lifetime'][0] 
        
        if self.metric_wh1001 == False:
            re = self.f.species_configs["greenhouse_gas_radiative_efficiency"].loc[dict(specie=self.minorghg)].data[0]
        else:
            re = gas_erf_diff[y]   #RE (1001 ensemble) doesn't change per SSP
            # expanding H.shape, so that each agwp[t] it incl. fair model uncertainty, not a point agwp[t]
            H = np.repeat(H, len(re)).reshape(-1, len(re)) 
            
        ppb2kg = 1e-9 * (mass/ self.M_AIR) * self.M_ATMOS
        A = re/ ppb2kg * ( 1 + halogen_ra) # W/m2/ppm -> W/m2/k
        
        agtp = H*0.
        iagtp = H*0.
        rf = H*0.
        agwp = H*0.
        rf = rf+A*np.exp(-H / alpha)
        agwp = agwp + A * alpha * (1-np.exp(-H / alpha))
        for j in np.arange(2):
            agtp = agtp + A *  alpha * q[j] * (np.exp(-H/ alpha) -
                                            np.exp(-H/d[j])) /\
                                       ( alpha-d[j])
            iagtp = iagtp + A * alpha * q[j] * \
                (alpha * ( 1 - np.exp(-H/( alpha)))-d[j] *
                 (1-np.exp(-H/d[j]))) / \
                (alpha-d[j])
        return rf, agwp, agtp, iagtp
    

    def get_minorghg_GWP(self, self_agwp, agwp_co2_folder = "output/metrics/"):
        """
        to calculate dGWP for minorGHG, need to fetch agwp_co2 from majorGHG modules (saved to output/metrics folder) 
        Input:
        --------
        self_agwp: agwp calculated from minorghg_anlytical()
        agwp_co2_folder: CO2 output folder name, reading SSP/MY for the CO2, AGWP_point value only, as minorGHG has only point AGWP
        Returns:
        --------
        final_gwp: pdGWP (only point value) for each minorGHG
        """
        sheet_y =  str(self.fair_start_y + self.year_index)
        co2_excel_file_name = "agwp_dcf_gwp" + str(self.H_max) + "_tstep" + str(self.ts_per_year) +  "CO2" + "_" + self.scn + "_fair_start" + str(self.fair_start_y) + "MY" + sheet_y + '.xlsx' 
        co2_excel_file_path = os.path.join(agwp_co2_folder, co2_excel_file_name)
      
        get_sheet = "agwp_pointvalueCO2_" if not self.metric_wh1001 else "agwp_wh_ensmb_CO2"
        agwp_co2_df = pd.read_excel(co2_excel_file_path, sheet_name = get_sheet + str(sheet_y), usecols=lambda column: column != "Unnamed: 0")  #ignore first index_col  
        if self.metric_wh1001 == False:
            agwp_co2 = agwp_co2_df["CO2"].values # w/o 1001, get single_col value to ensure same shape with self.agwp should be (101,) not (101, 1)
        else:
            agwp_co2 = agwp_co2_df.values  # if wh 1001, shape like (101, 1001) 
        print(f"final dims of AGWP for assessed minorGHG and CO2 is: ", self_agwp.shape, agwp_co2.shape )
        assert self_agwp.shape == agwp_co2.shape, "final dims of AGWP-minorGHG and CO2 not same, check if same FaIR startyear, ModelYear, H_max defined"
        final_gwp = np.ones_like(self_agwp) * 0
        final_gwp[1:, ] = self_agwp[1:,]  / (agwp_co2[1:,] )   #YR[0] has 0 AGWP, starting from YR[1]
        return final_gwp
    
    
    def get_dcf_finaloutput_minorghg (self, rf, agwp, gwp):
        """
        once final metrics calculated, save all output as dynamic CFs to excel
        Input:
        --------
        rf/agwp/gwp for each minorgas 
        Returns:
        --------
        Two versions for minorghg: 
        final_rf, final_agwp, final_gwp: same as input, wo 1001 ensemble, point value, saving to 3 excel sheets
        final_agwp, final_agwp_single, final_gwp, final_gwp_single, final_rf, final_rf_single: wh 1001 ensemble
        """
        tstep = self.ts_per_year
        HT = self.H_max         
        whichgas = str(self.minorghg) # to be used as excel names 
        #whichgas = re.sub(r'[^a-zA-Z0-9\s]', '', whichgas)  # remove special characters for excel naming
        
        if (tstep == 1):
            final_rf = rf        
            final_agwp = agwp  
            final_gwp = gwp
        else:  
            pass
            '''
            TBD
            '''
        if self.metric_wh1001 == True:  # calc point from 1001 ensemble, and more excel sheets to write 
            # get point value for each year
            final_rf_single = np.median(final_rf, axis=1)
            final_agwp_single = np.median(final_agwp, axis=1)
            final_gwp_single = np.median(final_gwp, axis=1)
        else: 
            pass   #final_agwp / final_gwp / final_rf are point values, as input
            
        # write metrics to excel
        dat = datetime.now().strftime('%Y-%m-%d')
        folder_name = f"{'output/metrics/minorGHG'}_{dat}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        
        sheet_y =  str(self.fair_start_y + self.year_index)
        excel_file_name = "agwp_dcf_gwp" + str(self.H_max) + "_tstep" + str(self.ts_per_year) +  whichgas + "_" + self.scn + "_fair_start" + str(self.fair_start_y) + "MY" + sheet_y + '.xlsx' 
        excel_file_path = os.path.join(folder_name, excel_file_name)

        with pd.ExcelWriter(excel_file_path, engine='xlsxwriter') as writer:
            if self.metric_wh1001 == True: 
                pd.DataFrame(final_agwp).to_excel(writer, sheet_name='agwp_wh_ensmb_' + whichgas  + sheet_y, index=True)
                pd.DataFrame(final_gwp).to_excel(writer, sheet_name='gwp_wh_ensmb_' + whichgas + sheet_y, index=True)
                pd.DataFrame(final_rf).to_excel(writer, sheet_name='rf_wh_ensmb_' + whichgas + sheet_y, index=True)
                #adding a col name for point-value 
                pd.DataFrame(final_agwp_single).rename(columns={pd.DataFrame(final_agwp_single).columns[0]: whichgas}).to_excel(writer, sheet_name='agwp_pointvalue' + whichgas+ '_'+ sheet_y, index=True)
                pd.DataFrame(final_gwp_single).rename(columns={pd.DataFrame(final_gwp_single).columns[0]: whichgas}).to_excel(writer, sheet_name='gwp_pointvalue' + whichgas+ '_' + sheet_y, index=True)
                pd.DataFrame(final_rf_single).rename(columns={pd.DataFrame(final_rf_single).columns[0]: whichgas}).to_excel(writer, sheet_name='rf_pointvalue' + whichgas+ '_' + sheet_y, index=True)

            else: #only point-value
                pd.DataFrame(final_agwp).rename(columns={pd.DataFrame(final_agwp).columns[0]: whichgas}).to_excel(writer, sheet_name='agwp_pointvalue' + whichgas+ '_'+ sheet_y, index=True)
                pd.DataFrame(final_gwp).rename(columns={pd.DataFrame(final_gwp).columns[0]: whichgas}).to_excel(writer, sheet_name='gwp_pointvalue' + whichgas+ '_' + sheet_y, index=True)
                pd.DataFrame(final_rf).rename(columns={pd.DataFrame(final_rf).columns[0]: whichgas}).to_excel(writer, sheet_name='rf_pointvalue' + whichgas+ '_' + sheet_y, index=True)

        print(f'calculated metric saved to {excel_file_path}')
        
        if self.metric_wh1001 == False: 
            return final_agwp, final_gwp, final_rf
        else: 
            return final_agwp, final_agwp_single, final_gwp, final_gwp_single, final_rf, final_rf_single
