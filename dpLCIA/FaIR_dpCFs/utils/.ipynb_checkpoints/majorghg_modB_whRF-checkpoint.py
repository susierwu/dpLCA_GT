import os 
import numpy as np
import pandas as pd 

class majorghg_analysis:
    """ once fair is run and output parameters ready through module A, this module B runs following functions:
    1. co2_analytical2023 <- modifying co2_analytical from IPCC_AR6, adding new input parameters ERF_per_unit_emission 
    2. ch4_analytical2023 
    3. n2o_analytical2023
    4. carbon cycle (CC) adjustment for CH4 and N2O  
    5. prepare final RF (W m-2 kg-1), AGWP and GWP for each gas 
    6. calculate dynamic characterization factor (DCF) as in A. Levasseur paper for dynmic life cycle assessment\
        and saved all output to folder "output/metric", this DCF won't be used in our study though
    """

    def __init__(
        self,
        f,
        #majorghg = ["CO2",  "CH4", "N2O"],
        scn,
        H_max=100,
        ts_per_year = 1,
        scenarios = ["ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp460", "ssp534-over", "ssp585"],
        fair_start_y = 1750,  # we used 2000 as the fair start year 
        year_index = 269,     # 1750 + 269 = year2019;  this year_index ->>  ModelYear
        #y_range=[2020,2030,2040,2050,2060],
        M_ATMOS = 5.1352E18,
        M_AIR = 28.97E-3,
        M_CO2 = 44.01E-3,
        M_C = 12.0E-3,
        M_CH4 = 16.043E-3,  
        M_N2O = 44.0E-3 # 44.01?
    ):
        """Initialise"""
        #self.gas = majorghg
        self.f = f
        self.scn = scn
        self.H_max = H_max
        self.year_index = year_index
        self.fair_start_y = fair_start_y
        self.ts_per_year = ts_per_year
        self.H = np.linspace(0, self.H_max, self.H_max * self.ts_per_year + 1) #parameter for co2_analytical()
        self.M_ATMOS = M_ATMOS
        self.M_AIR = M_AIR
        self.M_CO2 = M_CO2 
        self.M_C = M_C 
        self.M_CH4 = M_CH4 
        self.M_N2O = M_N2O
        
        # parameters for CC adjustment if combining Module A / B, can call self. directly, w/o specifying input parameters        
        # self.co2_erf_diff_t = self.get_co2_1ppm_erf()  #from module A
        # self.rf_co2 = self.co2_analytical2023()[0]
        # self.agwp_co2 = self.co2_analytical2023()[1]
        # self.agtp_co2 = self.co2_analytical2023()[2]
        # self.iagtp_co2 = self.co2_analytical2023()[3]
    


    def co2_analytical2023(self,  co2_erf_diff_t, 
                       co2=409.85, n2o=332.091, co2_ra=0.05,  #2019 base_concentration as in IPCC_ar6 
                       d=np.array([3.424102092311, 285.003477841911]), 
                       q=np.array([0.443767728883447, 0.313998206372015]), 
                       a=np.array([0.2173, 0.2240, 0.2824, 0.2763]),   #parameters for AGWP, same a0/a1/a2/a3 as in AR5  
                       alpha_co2=np.array([0, 394.4, 36.54, 4.304], 
                       )
                    ):
        
        """
        modifying IPCC_ar6::co2_analytical(), RE is no longer calculated from old fair meinshausen(), as: re = meinshausen(np.array([co2+1, 1866.3, n2o]), np.array([co2, 1866.3, n2o]), scale_F2x=False)[0] * (1+co2_ra) 
        now supplied through get_co2_meinshausen2020() function as co2_erf_diff_t[y]: 
        Inputs:
        -------
        co2_erf_diff_t: from Module A
        co2=409.85, n2o=332.091 the 2019 baseline concentration is not used, changing per SSP[x] MY[t], by the new RE = co2_erf_diff_t[y]
        # to to: d / q need to be scenario-depending for AGTP,  how?
        Returns:
        --------
        (rf, agwp, agtp, iagtp) : tuple of float or `np.ndarray`
            rf : Effective radiative forcing from a 1 ppmv increase in CO2
            agwp : Absolute global warming potential of CO2, W m-2 yr kg-1
            agtp : Absolute global temperature change potential of CO2, K kg-1
            iagtp : Integrated absolute global temperature change potential, K kg-1
        """

        # the CH4 concentration does not affect CO2 forcing, using hardcoded 2019 value from co2_analytical()
        y = self.year_index 
        print("analyzing for year:", (self.fair_start_y + y))
        
        re = co2_erf_diff_t[y]
        H = self.H
        # expanding H.shape, so that each agwp[t] it incl. fair model uncertainty, not a point agwp[t]
        H = np.repeat(H, len(re)).reshape(-1, len(re))  

        ppm2kg = 1E-6*(self.M_CO2/self.M_AIR)*self.M_ATMOS
        A = re/ppm2kg  # W/m2/ppm -> W/m2/kg

        agtp = H*0.
        iagtp = H*0.
        rf = H*0.
        agwp = H*0.
        #print(rf.shape, agwp.shape)
        for j in np.arange(2):
            if (j == 0):
                rf = rf+A*a[0]
                agwp = agwp+A*a[0]*H
            agtp = agtp+A*a[0]*q[j]*(1-np.exp(-H/d[j]))
            iagtp = iagtp+A*a[0]*q[j]*(H-d[j] * (1-np.exp(-H/d[j])))

            for i in np.arange(1, 4):
                if (j == 0):
                    rf = rf+A*a[i]*np.exp(-H/alpha_co2[i])
                    agwp = agwp+A*a[i]*alpha_co2[i] *\
                        (1-np.exp(-H/alpha_co2[i]))
                agtp = agtp+A*a[i]*alpha_co2[i]*q[j] *\
                    (np.exp(-H/alpha_co2[i]) -
                     np.exp(-H/d[j]))/(alpha_co2[i]-d[j])
                iagtp = iagtp+A*a[i]*alpha_co2[i]*q[j] * \
                    (alpha_co2[i]*(1-np.exp(-H/alpha_co2[i])) -
                     d[j]*(1-np.exp(-H/d[j]))) / \
                    (alpha_co2[i]-d[j])

        return rf, agwp, agtp, iagtp
    
    
    def ch4_analytical2023(self, alpha_ch4, ch4_erf_diff_t, 
                   co2=409.85, ch4=1866.3275, n2o=332.091, ch4_ra=-0.14, ch4_o3=1.4e-4, ch4_h2o=0.00004, 
                   d= np.array([3.424102092311, 285.003477841911]),
                   q= np.array([0.443767728883447, 0.313998206372015])): 
        """Calculates metrics for a 1 ppb CH4 perturbation.
        Inputs:
        -------
        ch4_erf_diff_t : replacing re from meinshausen(), calculated from majorghg_get_f.get_ch4_1ppb_erf as W/m2/ppb CH4
        alpha_ch4: output from majorghg_get_f.call_f_from_fair_gas(), changing (not a single float) perturbation lifetime of CH4, for SSPx,year[x]
        d / q array from "335_chapter7_generate_metrics.ipynb", it should be also SSP-depending for AGTP?
        co2=409.85, ch4=1866.3275 ...not used, kept here as in original ch4_analytical() function, can delete
        """
        ppb2kg = 1e-9*(self.M_CH4/self.M_AIR)*self.M_ATMOS
        y = self.year_index 
        re = ch4_erf_diff_t[y]
        H =  self.H 
        H = np.repeat(H, len(re)).reshape(-1, len(re))  

        A = (re + ch4_o3 + ch4_h2o)/ppb2kg  #(1001,)
        #print(alpha_ch4.shape, re.shape, H.shape )    
        #print("A shape:", A.shape)

        agtp = H*0.
        iagtp = H*0.
        rf = H*0.
        agwp = H*0.

        rf = rf+A*np.exp(-H/(alpha_ch4))

        agwp = agwp+A*alpha_ch4*(1-np.exp(-H/alpha_ch4))
        for j in np.arange(2):
            agtp = agtp+A*alpha_ch4*q[j] *\
                (np.exp(-H/(alpha_ch4)) -
                 np.exp(-H/d[j]))/(alpha_ch4-d[j])
            iagtp = iagtp+A*alpha_ch4*q[j] * \
                (alpha_ch4*(1-np.exp(-H/(alpha_ch4))) -
                 d[j]*(1-np.exp(-H/d[j]))) / \
                (alpha_ch4-d[j])
        return rf, agwp, agtp, iagtp 
    
    
    def n2o_analytical2023(self, alpha_n2o, ch4_erf_diff_t, n2o_erf_diff_t,
                    co2=409.85, ch4=1866.3275, n2o=332.091, n2o_ra=0.07, 
                   n2o_o3=5.5e-4, f_n2o_ch4=-1.7, ch4_ra=-0.14, ch4_o3=1.4e-4, ch4_h2o=0.00004, 
                   d=np.array([3.424102092311, 285.003477841911]),
                   q = np.array([0.443767728883447, 0.313998206372015]),   # this is what Bill used 
                  ):
        """Calculates metrics for a 1 ppb N2O perturbation.
        Inputs:
        -------
        alpha_n2o : output from ModuleA, changing lifetime of N2O, for SSPx,year[x], we keep it as 109 here
        ch4_erf_diff_t, replacing re from meinshausen(), calculated from majorghg_get_f.get_ch4_1ppb_erf as W/m2/ppb CH4
        n2o__erf_diff_t,replacing re from meinshausen(), calculated from majorghg_get_f.get_n2o_1ppb_erf as W/m2/ppb N2O
        d / q array from "335_chapter7_generate_metrics.ipynb", it should be also ssp-depending for AGTP?
        Returns:
        --------
        (rf, agwp, agtp, iagtp) : tuple of float or `np.ndarray`
            rf : Effective radiative forcing from a 1 ppbv increase in CH4
            agwp : Absolute global warming potential of CH4, W m-2 yr kg-1
            agtp : Absolute global temperature change potential of CH4, K kg-1
            iagtp : Integrated absolute global temperature change potential, K kg-1
        """
        y = self.year_index 
        re_n2o = n2o_erf_diff_t[y]
        re_ch4 = ch4_erf_diff_t[y]
        H =  self.H 
        H = np.repeat(H, len(re_n2o)).reshape(-1, len(re_n2o))  
        ppb2kg = 1e-9*(self.M_N2O/self.M_AIR)*self.M_ATMOS
        # Add in a component for the destruction of methane from AR5 8.SM.11.3.3
        A = (re_n2o+f_n2o_ch4*re_ch4)/ppb2kg
        #print(alpha_n2o.shape, re_n2o.shape, re_ch4.shape, H.shape )    
        #print("A shape:", A.shape)

        agtp = H*0.
        iagtp = H*0.
        rf = H*0.
        agwp = H*0.
        rf = rf+A*np.exp(-H/alpha_n2o)
        agwp = agwp+A*alpha_n2o*(1-np.exp(-H/alpha_n2o))
        for j in np.arange(2):
            agtp = agtp+A*alpha_n2o*q[j]*(np.exp(-H/alpha_n2o) -
                                                      np.exp(-H/d[j])) /\
                                                 (alpha_n2o-d[j])
            iagtp = iagtp+A*alpha_n2o*q[j] * \
                (alpha_n2o*(1-np.exp(-H/(alpha_n2o))) -
                 d[j]*(1-np.exp(-H/d[j]))) / \
                (alpha_n2o-d[j])
        return rf, agwp, agtp, iagtp
    
    
    
    # adding carbon cycle adjustment 
    def carbon_cycle_adjustment(self, agtp, 
                                rf_co2, agwp_co2, agtp_co2, iagtp_co2 , 
                                co2=409.85, n2o=332.091, co2_ra=0.05, 
                               d= np.array([3.424102092311, 285.003477841911]),
                               q= np.array([0.443767728883447, 0.313998206372015])):
        """Calculates adjustment to metrics based on carbon cycle feedback
        Inputs:
        -------
        agtp: unadjusted AGTP for CH4/N2O, 3rd output from ch4_analytical2023()[2] & n2o_analytical2023()[2], with shape(H,1001) 
        all a / d /q / alpha/ gamma parameters same as "https://github.com/chrisroadmap/ar6/blob/main/src/ar6/metrics/gasser.py#L13"
        alpha [2.376, 30.14, 490.1] hardcoded below need to adjust these parameters per SSP/year? 
        Returns:
        --------
        (rf_cc, agwp_cc, agtp_cc) : tuple of `np.ndarray`
            rf_cc : Increase in effective radiative forcing due to carbon cycle adjustment
            agwp : Increase in absolute global warming potential due to carbon cycle adjustment, W m-2 yr kg-1
            agtp : Increase in absolute global temperature change potential due to carbon cycle adjustment, K kg-1
        """
        
        H = self.H                        
        H = np.repeat(H, agtp.shape[1]).reshape(-1, agtp.shape[1])  #each year with 1001 runs (agtp.shape[1])
        dts = H[1]
                                                
        #rf_co2, agwp_co2, agtp_co2, iagtp_co2  = self.rf_co2,  self.agwp_co2, self.agtp_co2, self.iagtp_co2
        #print(H.shape, dts.shape, agwp_co2.shape, agtp_co2.shape)
                                
        agtp_cc = H*0.
        agwp_cc = H*0.
        rf_cc = H*0.
        F_CO2 = H*0.
                                
        a = np.array([0.6368, 0.3322, 0.0310])  # Gasser et al. 2017
        alpha = np.array([2.376, 30.14, 490.1])
        gamma = 3.015*1E12  # kgCO2/yr/K  Gasser et al. 2017
        r_f = H*0.
        r_f[0] = np.sum(a)/dts
        for i in np.arange(0, 3):
            r_f = r_f-(a[i]/alpha[i])*np.exp(-H/alpha[i])

        for j in np.arange(H.shape[0]):  # "H.size" if point value for each year
            for i in np.arange(j+1):
                F_CO2[j] = F_CO2[j]+agtp[i]*gamma*r_f[j-i]*dts
        for j in np.arange(H.shape[0]):
            for i in np.arange(j+1):
                rf_cc[j] = rf_cc[j]+F_CO2[i]*rf_co2[j-i]*dts * \
                    (self.M_CO2/self.M_C)
                agwp_cc[j] = agwp_cc[j]+F_CO2[i]*agwp_co2[j-i]*dts * \
                    (self.M_CO2/self.M_C)
        return rf_cc, agwp_cc, agtp_cc
    
    
    
    ### get final RF, AGWP and GWP wh CC_adjustment for N2O and CH4
    def get_final_rf_agwp_gwp_ch4_n2o(self, rf_gas, rf_cc_gas, agwp_gas, agwp_cc_gas, agwp_co2):  
        """get final AGWP and GWP with CC for CH4, N2O
        Input:  
        -------
            rf_gas ::                               <- ch4_/n2o_analytical2023()[0]
            rf_cc_gas ::                            <- carbon_cycle_adjustment()[0]
            agwp_gas::agwp_ch4, agwp_n2o,           <- ch4_/n2o_analytical2023()[1]
            agwp_cc_gas: agwp_cc_ch4, agwp_cc_n2o,  <- carbon_cycle_adjustment()[1]
        Returns:  
        -------
        """
        final_rf = rf_gas + rf_cc_gas
        final_agwp = agwp_gas  +  agwp_cc_gas
        final_gwp = np.ones_like(final_agwp) * 0         #YR[0] has 0 AGWP_CO2 as denominator, using same [0,...] as AGWP
        final_gwp[1:, ] =  final_agwp[1:,]  / (agwp_co2[1:,] )   #starting from YR[1]
        return final_agwp,  final_gwp, final_rf
     
    def final_rf_agwp_gwp_co2(self, rf_co2, agwp_co2): 
        final_rf_moi = rf_co2
        final_agwp_moi = agwp_co2
        final_gwp_moi = np.ones_like(agwp_co2) * 1
        return final_agwp_moi, final_gwp_moi, final_rf_moi
        
    
    ### get instn. dynamic characterization factor (DCF[t]) for each gas as in Annie' paper: 
    """ https://pubs.acs.org/doi/full/10.1021/es9030003# 
    To obtain  the  instantaneous  value  of  the  characterization factor for any year following the emission, shown in Figure3b, 
    the time scale was divided into one-year time  periods and the integration boundaries were set for every time step,as shown in eq 5 
    """
    
    def get_dcf_finaloutput (self, whichgas, rf, agwp, gwp):

        """
        once final metrics calculated, save all output to excel
        When tstep == 1, final_agwp == agwp
        When tstep != 1, final_agwp should be annual data, not same agwp.dims, as each year has tstep
        Input:
        --------
         whichgas = ["CO2", "CH4", "N2O"]
         rf/agwp/gwp for the gas directly prepared in get_finalagwp_gwp_ch4_n2o() | final_agwp_gwp_co2()
        Returns:
        --------
        (final_rf, final_rf_single, final_agwp, final_agwp_single, final_dcf, final_dcf_single, final_gwp, final_gwp_single): `np.ndarray`
            final_rf: final adjusted inst. radiative forcing (IRF)
            final_agwp : final AGWP[t] with uncertainty from a unit increase in gas, W m-2 yr kg-1 
            final_agwp_single : median of (final_agwp) for AGWP[t]
            final_dcf: annaul dynamic CF with uncertainty for dGWP calculation (Annie's metric, we don't use it here)
            final_dcf_single:  median of (final_dcf)
        the final_single value can be used to cross-check against Annie's metrics:
            1 annual DCF data at "A.Levasseur (2010) https://pubs.acs.org/doi/epdf/10.1021/es9030003"  
                1.1 validation excel: https://ciraig.org/index.php/project/dynco2-dynamic-carbon-footprinter/ 
            2.AGWP20/100 value in IPCC AR6 Table 7.SM.7
        -------- 
        """
        tstep = self.ts_per_year
        HT = self.H_max         
        final_rf_single = np.ones_like(range(rf.shape[0])) * np.nan
        final_agwp_single = np.ones_like(range(agwp.shape[0])) * np.nan   #one point value for each single year
        final_gwp_single = np.ones_like(range(gwp.shape[0])) * np.nan 
        final_dcf_single = np.ones_like(range(agwp.shape[0])) * np.nan 
        
        if (tstep == 1):
            final_rf = rf        
            final_agwp = agwp   #from final_agwp_gwp_ch4_n2o(), adding AGWP_self + AGWP_CC
            final_gwp = gwp
        else:  
            pass
            '''
            # to improve, below code is not correct 
            shape = (self.H_max+1, 1001)   #1001 is the fair ensemble number can be fixed
            final_rf, final_agwp, final_gwp = np.full(shape, np.nan), np.full(shape, np.nan), np.full(shape, np.nan)
            #final_agwp = np.ones_like(agwp)/tstep * np.nan
            #final_gwp = np.ones_like(agwp)/tstep * np.nan 
            #print(final_agwp.shape)  
            for i in range(0,HT+1): 
                rf_t[i, ] = rf[i * tstep]
                final_rf[i, ] = rf_t
                agwp_t =  agwp[i * tstep]
                final_agwp[i, ] = agwp_t
                gwp_t = gwp[i * tstep]
                final_gwp[i, ] = gwp_t   #final_gwp.shape = (H_max*tstep, 1001)
            '''     
        # each year w/h N element: len = agwp36.shape[1], repeat N times 0 and adding to YEAR_0 
        final_dcf = np.diff(final_agwp, axis = 0) # axis = 0 diff by row : AGWP[t,j] - AGWP[t-1, j] 
        to_insert = np.array(np.repeat(0, final_agwp.shape[1]))
        final_dcf = np.insert(final_dcf, 0, to_insert, axis=0)  
        print(f"original AGWP and GWP dims is: ", agwp.shape, gwp.shape, "Hmax and ts per year is:",  HT, tstep )
        print(f"final dims is: ", final_dcf.shape, final_agwp.shape, final_gwp.shape )
        print("check dimentions:", final_dcf.shape == final_agwp.shape == final_gwp.shape ) 
        # get point value for each year
        final_rf_single = np.median(final_rf, axis=1)
        final_agwp_single = np.median(final_agwp, axis=1)
        final_dcf_single = np.median(final_dcf, axis=1)
        final_gwp_single = np.median(final_gwp, axis=1)
        
        # write metrics to excel
        if not os.path.exists('output/metrics'):
            os.makedirs('output/metrics')
            
        sheet_y =  str(self.fair_start_y + self.year_index)
        excel_file_name = "agwp_dcf_gwp" + str(self.H_max) + "_tstep" + str(self.ts_per_year) +  whichgas + "_" + self.scn + "_fair_start" + str(self.fair_start_y) + "MY" + sheet_y + '.xlsx' 
        output_folder = 'output/metrics'   
        excel_file_path = os.path.join(output_folder, excel_file_name)

        with pd.ExcelWriter(excel_file_path, engine='xlsxwriter') as writer:
            pd.DataFrame(final_agwp).to_excel(writer, sheet_name='agwp_wh_ensmb_' + whichgas  + sheet_y, index=True)
            pd.DataFrame(final_dcf).to_excel(writer, sheet_name='dcf_wh_ensmb_' + whichgas + sheet_y, index=True)
            pd.DataFrame(final_gwp).to_excel(writer, sheet_name='gwp_wh_ensmb_' + whichgas + sheet_y, index=True)
            pd.DataFrame(final_rf).to_excel(writer, sheet_name='rf_wh_ensmb_' + whichgas + sheet_y, index=True)
            
            #adding a col name for point-value AGWP and DCF
            pd.DataFrame(final_agwp_single).rename(columns={pd.DataFrame(final_agwp_single).columns[0]: whichgas}).to_excel(writer, sheet_name='agwp_pointvalue' + whichgas+ '_'+ sheet_y, index=True)
            pd.DataFrame(final_dcf_single).rename(columns={pd.DataFrame(final_dcf_single).columns[0]: whichgas}).to_excel(writer, sheet_name='dcf_pointvalue' + whichgas+ '_' + sheet_y, index=True)
            pd.DataFrame(final_gwp_single).rename(columns={pd.DataFrame(final_gwp_single).columns[0]: whichgas}).to_excel(writer, sheet_name='gwp_pointvalue' + whichgas+ '_' + sheet_y, index=True)
            pd.DataFrame(final_rf_single).rename(columns={pd.DataFrame(final_rf_single).columns[0]: whichgas}).to_excel(writer, sheet_name='rf_pointvalue' + whichgas+ '_' + sheet_y, index=True)

        print(f'calculated metric saved to {excel_file_path}')
        return final_agwp, final_agwp_single, final_dcf, final_dcf_single, final_gwp, final_gwp_single, final_rf, final_rf_single
    
