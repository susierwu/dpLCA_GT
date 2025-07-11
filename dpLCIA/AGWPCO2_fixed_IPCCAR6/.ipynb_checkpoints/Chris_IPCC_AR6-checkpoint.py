"""
All functions from : Chris IPCC roadmap
Module for calculating metrics from CO2, usually as a baseline to compare other gases.
Author: Bill Collins (UK)
Adapted by Chris Smith
"""

import numpy as np
from fair.forcing.ghg import meinshausen
from fair.defaults.thermal import q, d         #hard-coded below
from fair.defaults.carbon import a, tau as alpha_co2
d = np.array([3.424102092311, 285.003477841911]),        #parameters used to calculate AGTP
q = np.array([0.443767728883447, 0.313998206372015])  

a=np.array([0.2173, 0.2240, 0.2824, 0.2763]) 
alpha_co2=np.array([0, 394.4, 36.54, 4.304])


M_ATMOS = 5.1352E18
M_AIR = 28.97E-3
M_CO2 = 44.01E-3
M_C = 12.0E-3  # 12.01?
M_CH4 = 16.043E-3
M_N2O = 44.0E-3 # 44.01?


# TODO: make better variable names
def co2_analytical(H, co2=409.85, n2o=332.091, co2_ra=0.05, d=d, q=q, a=a, alpha_co2=alpha_co2):
    """Calculates baseline metrics for a 1 ppm CO2 perturbation.
    
    Inputs:
    -------
    H : float or `np.ndarray`
        time horizon(s) of interest
    co2 : float, optional
        baseline concentrations of CO2, ppmv
    n2o : float, optional
        baseline concentrations of N2O, ppbv
    co2_ra : float, optional
        tropospheric rapid adjustment enhancement of CO2 forcing, expressed as a decimal
    d : `np.ndarray`, optional
        2-element array of fast and slow timescales to climate warming impulse response function
    q : `np.ndarray`, optional
        2-element array of fast and slow contributions to climate warming impulse response function
    a : `np.ndarray`, optional
        4-element array of partition fractions of CO2 atmospheric boxes, slow to fast
    alpha_co2 : `np.ndarray`, optional
        4-element array of time constants of CO2 atmospheric boxes, slow to fast
    
    Returns:
    --------
    (rf, agwp, agtp, iagtp) : tuple of float or `np.ndarray`
        rf : Effective radiative forcing from a 1 ppmv increase in CO2
        agwp : Absolute global warming potential of CO2, W m-2 yr kg-1
        agtp : Absolute global temperature change potential of CO2, K kg-1
        iagtp : Integrated absolute global temperature change potential, K kg-1
    """
    # the CH4 concentration does not affect CO2 forcing, so we hardcode an approximate 2019 value
    re = meinshausen(np.array([co2+1, 1866.3, n2o]), np.array([co2, 1866.3, n2o]), scale_F2x=False)[0] * (1+co2_ra)

    ppm2kg = 1E-6*(M_CO2/M_AIR)*M_ATMOS
    A = re/ppm2kg  # W/m2/kg

    agtp = H*0.
    iagtp = H*0.
    rf = H*0.
    agwp = H*0.
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






def ch4_analytical(H, co2=409.85, ch4=1866.3275, n2o=332.091, ch4_ra=-0.14, ch4_o3=1.4e-4, ch4_h2o=0.00004, d=d, q=q, alpha_ch4=11.8):
    """Calculates metrics for a 1 ppb CH4 perturbation.
    
    Inputs:
    -------
    H : float or `np.ndarray`
        time horizon(s) of interest
    co2 : float, optional
        baseline concentrations of CO2, ppmv
    ch4: float, optional
        baseline concentrations of CH4, ppbv
    n2o : float, optional
        baseline concentrations of N2O, ppbv
    ch4_ra : float, optional
        tropospheric rapid adjustment enhancement of CH4 forcing
    ch4_o3 : float, optional
        radiative efficiency increase of CH4 emissions due to O3 formation, W m-2 (ppb CH4)-1
    ch4_h2o : float, optional
        radiative efficiency increase of CH4 emissions due to stratospheric H2O formation, W m-2 (ppb CH4)-1 
    d : `np.ndarray`, optional
        2-element array of fast and slow timescales to climate warming impulse response function
    q : `np.ndarray`, optional
        2-element array of fast and slow contributions to climate warming impulse response function
    alpha_ch4 : float
        perturbation lifetime of CH4, years
    
    Returns:
    --------
    (rf, agwp, agtp, iagtp) : tuple of float or `np.ndarray`
        rf : Effective radiative forcing from a 1 ppbv increase in CH4
        agwp : Absolute global warming potential of CH4, W m-2 yr kg-1
        agtp : Absolute global temperature change potential of CH4, K kg-1
        iagtp : Integrated absolute global temperature change potential, K kg-1
    """
    re = meinshausen(np.array([co2, ch4+1, n2o]), np.array([co2, ch4, n2o]), scale_F2x=False)[1] * (1+ch4_ra)
    ppb2kg = 1e-9*(M_CH4/M_AIR)*M_ATMOS
    A = (re + ch4_o3 + ch4_h2o)/ppb2kg

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







# TODO: make better variable names
def n2o_analytical(H, co2=409.85, ch4=1866.3275, n2o=332.091, n2o_ra=0.07, n2o_o3=5.5e-4, f_n2o_ch4=-1.7, ch4_ra=-0.14, ch4_o3=1.4e-4, ch4_h2o=0.00004, d=d, q=q, alpha_n2o=109):
    """Calculates metrics for a 1 ppb N2O perturbation.
    
    Inputs:
    -------
    H : float or `np.ndarray`
        time horizon(s) of interest
    co2 : float, optional
        baseline concentrations of CO2, ppmv
    ch4: float, optional
        baseline concentrations of CH4, ppbv
    n2o : float, optional
        baseline concentrations of N2O, ppbv
    n2o_ra : float, optional
        tropospheric rapid adjustment enhancement of N2O forcing
    n2o_o3 : float, optional
        radiative efficiency increase of N2O emissions due to O3 formation, W m-2 (ppb N2O)-1
    f_n2o_ch4 : float, optional
        feedback change in methane lifetime due to N2O emissions, (ppb CH4) (ppb N2O)-1
    ch4_ra : float, optional
        tropospheric rapid adjustment enhancement of CH4 forcing
    ch4_o3 : float, optional
        radiative efficiency increase of CH4 emissions due to O3 formation, W m-2 (ppb CH4)-1
    ch4_h2o : float, optional
        radiative efficiency increase of CH4 emissions due to stratospheric H2O formation, W m-2 (ppb CH4)-1 
    d : `np.ndarray`, optional
        2-element array of fast and slow timescales to climate warming impulse response function
    q : `np.ndarray`, optional
        2-element array of fast and slow contributions to climate warming impulse response function
    alpha_n2o : float
        perturbation lifetime of N2O, years
    
    Returns:
    --------
    (rf, agwp, agtp, iagtp) : tuple of float or `np.ndarray`
        rf : Effective radiative forcing from a 1 ppbv increase in CH4
        agwp : Absolute global warming potential of CH4, W m-2 yr kg-1
        agtp : Absolute global temperature change potential of CH4, K kg-1
        iagtp : Integrated absolute global temperature change potential, K kg-1
    """
    re_n2o = meinshausen(np.array([co2, ch4, n2o+1]), np.array([co2, ch4, n2o]), scale_F2x=False)[2] * (1+n2o_ra) + n2o_o3
    re_ch4 = meinshausen(np.array([co2, ch4+1, n2o]), np.array([co2, ch4, n2o]), scale_F2x=False)[1] * (1+ch4_ra) + ch4_o3+ch4_h2o
    ppb2kg = 1e-9*(M_N2O/M_AIR)*M_ATMOS
# Add in a component for the destruction of methane from AR5 8.SM.11.3.3
    A = (re_n2o+f_n2o_ch4*re_ch4)/ppb2kg

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