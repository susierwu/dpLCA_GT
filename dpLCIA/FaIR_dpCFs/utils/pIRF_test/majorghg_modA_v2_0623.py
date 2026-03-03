import os
import numpy as np
import pandas as pd
# ======================================================================
# CHANGE (fix): prospective decay side (`pIRF`) via perturbed FaIR rerun
#   - FIX 1: robust matching for FaIR float timepoints (e.g., 2020.5)
#   - FIX 2: normalize IRF using the first post-pulse ΔC (discrete grid),
#            and enforce IRF[0] = 1 to avoid early-time IRF > 1 artifacts.
# ======================================================================
from typing import Optional, Union
import numpy as np
import copy

"""
def _as_float_array(x) -> np.ndarray:
    return np.asarray(x, dtype=float)

def _find_index_near(arr: Union[np.ndarray, list], target: float, tol: float = 1e-6) -> int:
    #Find index i where arr[i] is closest to target, within tolerance.
    #Works for float arrays like FaIR timepoints (1750.5, 1751.5, ...).
    a = _as_float_array(arr)
    i = int(np.argmin(np.abs(a - float(target))))
    if abs(a[i] - float(target)) > tol:
        raise ValueError(f"Target={target} not found within tol={tol}. Closest={a[i]} at i={i}.")
    return i

def _find_index_exact_int(arr: Union[np.ndarray, list], year_int: int) -> int:
    
    #Exact match on integer-like year arrays (FaIR timebounds typically 1750., 1751., ...).
    a = _as_float_array(arr)
    hits = np.where(a == float(int(year_int)))[0]
    if len(hits) == 0:
        raise ValueError(f"Year {year_int} not found in year axis.")
    return int(hits[0])


def _get_dt_years(ts_per_year: int) -> float:
    return 1.0 / float(ts_per_year)

"""

class majorghg_get_f:
    """In order to calculate climate metrics, this module A runs:
    1) extract major GHG FaIR output parameters (C, baseline, forcing_scale, lifetime scale)
    2) compute ERF using Meinshausen2020 parameterization (FAIR: src/fair/forcing/ghg.py)
    3) compute ERF per unit concentration perturbation ([W m-2 ppm-1] CO2, [W m-2 ppb-1] CH4/N2O)
    4) CHANGE/ADDED in v2: compute prospective impulse response functions pIRF via a perturbed FaIR re-run:
       DeltaC(t) = C_pert(t) - C_base(t), IRF(t) = DeltaC(t) / DeltaC(0)
    """

    def __init__(
        self,
        f,
        scn,
        H_max=100,
        ts_per_year=1,
        scenarios=[
            "ssp119", "ssp126", "ssp245", "ssp370",
            "ssp434", "ssp460", "ssp534-over", "ssp585"
        ],
        year_index=269,  # 1750 + 269 = year2019 for annual time axis
    ):
        """Initialise"""
        self.f = f
        self.scn = scn
        self.H_max = H_max
        self.year_index = year_index
        self.ts_per_year = ts_per_year
        self.H = np.linspace(0, self.H_max, self.H_max * self.ts_per_year + 1)

        # core FaIR-derived parameter dicts
        self.co2_f, self.ch4_f, self.n2o_f = self.call_f_from_fair_gas()
        self.n2o_c = self.n2o_f["N2O_C"]  # retained for CO2 Meinshausen coupling

    # =========================
    # step A.1: call from fair, get parameter needed for climate metric calculation
    # =========================
    def call_f_from_fair_gas(self):
        """
        Call FaIR outputs required for metric calculation.

        Returns
        -------
        co2_f, ch4_f, n2o_f : dict
            Each dict contains:
              *_C: concentration (timebounds x ensemble)
              *_C_plus1: concentration + 1 ppm/ppb (local derivative trick for RE)
              *_RF: baseline_concentration broadcast to output shape
              *_FS: forcing_scale broadcast to output shape
              *_lifetime: (CH4/N2O) scenario-conditioned lifetime at year_index (ensemble)
            CHANGE: also reserves keys for prospective DeltaC/IRF (pIRF) outputs.
        """
        # =========================
        # CHANGE: add DeltaC/IRF keys for all gases
        # =========================
        co2_key = [
            "CO2_C", "CO2_C_plus1", "CO2_RF", "CO2_FS", "N2O_C",
            "CO2_dC_ppm_t",       # added
            "CO2_IRF_fair_t",     # added
        ]
        ch4_key = [
            "CH4_C", "CH4_C_plus1", "CH4_RF", "CH4_FS", "N2O_C", "CH4_lifetime",
            "CH4_dC_ppb_t",       # added
            "CH4_IRF_fair_t",     # added
        ]
        n2o_key = [
            "N2O_C", "N2O_C_plus1", "N2O_RF", "N2O_FS", "N2O_lifetime",
            "N2O_dC_ppb_t",        # added
            "N2O_IRF_fair_t",      # added
        ]

        co2_f = {key: None for key in co2_key}
        ch4_f = {key: None for key in ch4_key}
        n2o_f = {key: None for key in n2o_key}

        concentration = self.f.concentration.loc[dict(scenario=self.scn)]

        for gas in ["CO2", "CH4", "N2O"]:
            gas_c = concentration.loc[dict(specie=gas)].values
            gas_c_plus1 = gas_c + 1

            gas_rf = (
                np.ones_like(gas_c)
                * self.f.species_configs["baseline_concentration"].loc[dict(specie=gas)].values
            )
            gas_fs = (
                np.ones_like(gas_c_plus1)
                * self.f.species_configs["forcing_scale"].loc[dict(specie=gas)].values
            )

            # N2O needed for CO2
            n2o_c = concentration.loc[dict(specie="N2O")].values

            if gas == "CO2":
                co2_f["CO2_C"] = gas_c
                co2_f["CO2_C_plus1"] = gas_c_plus1
                co2_f["CO2_RF"] = gas_rf
                co2_f["CO2_FS"] = gas_fs
                co2_f["N2O_C"] = n2o_c

            elif gas == "CH4":
                # lifetime scaling at pulse year
                alpha_s = self.f.alpha_lifetime.loc[dict(scenario=self.scn, specie="CH4")][self.year_index].values
                ch4_unp = 10.8537568
                alpha_ch4 = ch4_unp * alpha_s
                ch4_f["CH4_lifetime"] = alpha_ch4

                ch4_f["CH4_C"] = gas_c
                ch4_f["CH4_C_plus1"] = gas_c_plus1
                ch4_f["CH4_RF"] = gas_rf
                ch4_f["CH4_FS"] = gas_fs
                ch4_f["N2O_C"] = n2o_c

            elif gas == "N2O":
                alpha_s = self.f.alpha_lifetime.loc[dict(scenario=self.scn, specie="N2O")][self.year_index].values
                n2o_unp = 109.0
                alpha_n2o = n2o_unp * alpha_s
                n2o_f["N2O_lifetime"] = alpha_n2o

                n2o_f["N2O_C"] = gas_c
                n2o_f["N2O_C_plus1"] = gas_c_plus1
                n2o_f["N2O_RF"] = gas_rf
                n2o_f["N2O_FS"] = gas_fs

        return co2_f, ch4_f, n2o_f

    # =========================
    # step A.2: get co2 / ch4 / n2o ERF (Meinshausen2020)
    # =========================
    def get_co2_meinshausen2020(
        self,
        co2_c_orcplus,
        a1=-2.4785e-07, b1=0.00075906, c1=-0.0021492, d1=5.2488,
        a2=-0.00034197, b2=0.00025455, c2=-0.00024357, d2=0.12173,
        a3=-8.9603e-05, b3=-0.00012462, d3=0.045194,
    ):
        """
        Inputs
        ------
        co2_c_orcplus: concentration array (timebounds x ensemble), baseline or +1 ppm

        Returns
        -------
        erf_co2: ERF array (timebounds x ensemble)
        """
        co2 = co2_c_orcplus
        co2_base = self.co2_f["CO2_RF"]
        n2o = self.n2o_c
        co2_fs = self.co2_f["CO2_FS"]

        ca_max = co2_base - b1 / (2 * a1)

        where_central = np.asarray((co2_base < co2) & (co2 <= ca_max)).nonzero()
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

        erf_co2 = (alpha_p + alpha_n2o) * np.log(co2 / co2_base) * co2_fs
        return erf_co2

    def get_ch4_meinshausen2020(
        self,
        ch4_c_orcplus,
        a1=-2.4785e-07, b1=0.00075906, c1=-0.0021492, d1=5.2488,
        a2=-0.00034197, b2=0.00025455, c2=-0.00024357, d2=0.12173,
        a3=-8.9603e-05, b3=-0.00012462, d3=0.045194,
    ):
        """
        Inputs
        ------
        ch4_c_orcplus: CH4 concentration array (timebounds x ensemble), baseline or +1 ppb
        """
        ch4 = ch4_c_orcplus
        ch4_base = self.ch4_f["CH4_RF"]
        n2o = self.ch4_f["N2O_C"]
        ch4_fs = self.ch4_f["CH4_FS"]

        erf_ch4 = (
            (a3 * np.sqrt(ch4) + b3 * np.sqrt(n2o) + d3)
            * (np.sqrt(ch4) - np.sqrt(ch4_base))
        ) * ch4_fs

        return erf_ch4

    def get_n2o_meinshausen2020(
        self,
        n2o_c_orcplus,
        a1=-2.4785e-07, b1=0.00075906, c1=-0.0021492, d1=5.2488,
        a2=-0.00034197, b2=0.00025455, c2=-0.00024357, d2=0.12173,
        a3=-8.9603e-05, b3=-0.00012462, d3=0.045194,
    ):
        """
        Inputs
        ------
        n2o_c_orcplus: N2O concentration array (timebounds x ensemble), baseline or +1 ppb
        """
        n2o = n2o_c_orcplus
        n2o_base = self.n2o_f["N2O_RF"]
        ch4 = self.ch4_f["CH4_C"]
        co2 = self.co2_f["CO2_C"]
        n2o_fs = self.n2o_f["N2O_FS"]

        erf_n2o = (
            (a2 * np.sqrt(co2) + b2 * np.sqrt(n2o) + c2 * np.sqrt(ch4) + d2)
            * (np.sqrt(n2o) - np.sqrt(n2o_base))
        ) * n2o_fs

        return erf_n2o

    # =========================
    # step A.3 get erf per 1 unit concentration perturbation [W m-2 ppm/ppb-1]
    # =========================
    def get_co2_1ppm_erf(self):
        """
        Returns
        -------
        co2_erf_diff: ERF difference for +1 ppm (timebounds x ensemble)
        """
        erf1 = self.get_co2_meinshausen2020(self.co2_f["CO2_C"])
        erf2 = self.get_co2_meinshausen2020(self.co2_f["CO2_C_plus1"])
        return erf2 - erf1

    def get_ch4_1ppb_erf(self):
        erf1 = self.get_ch4_meinshausen2020(self.ch4_f["CH4_C"])
        erf2 = self.get_ch4_meinshausen2020(self.ch4_f["CH4_C_plus1"])
        return erf2 - erf1

    def get_n2o_1ppb_erf(self):
        erf1 = self.get_n2o_meinshausen2020(self.n2o_f["N2O_C"])
        erf2 = self.get_n2o_meinshausen2020(self.n2o_f["N2O_C_plus1"])
        return erf2 - erf1

    # =========================
    # step A.4-A: alpha_lifetime for CH4 and N2O
    # =========================
    def get_ch4_n2o_alpha(self):
        """Return lifetime arrays (ensemble) for CH4 and N2O at year_index."""
        return self.ch4_f["CH4_lifetime"], self.n2o_f["N2O_lifetime"]

    # =========================
    # step A.4-B: concentration arrays (timebounds x ensemble)
    # =========================
    def get_majorghg_concentration(self):
        """Return concentrations (timebounds x ensemble) for CO2, CH4, N2O."""
        return self.co2_f["CO2_C"], self.ch4_f["CH4_C"], self.n2o_f["N2O_C"]

    
    # -----------------------------------------------------------------------------------------------------------------------------
    # FIXED pulse-year handling
    # -----------------------------------------------------------------------------------------------------------------------------
# ======================================================================
# CHANGE: prospective decay side (`pIRF`) via perturbed FaIR rerun (FAST)
# ======================================================================
    def _get_pulse_year(self) -> int:
        """
        CHANGE (fix): pulse year is the *calendar year* on timebounds at year_index.
        With annual FaIR grids: timebounds[k] = year, timepoints[k] = year + 0.5.
        """
        if hasattr(self.f, "timebounds"):
            return int(np.asarray(self.f.timebounds)[self.year_index])
        return int(self.year_index)


    @staticmethod
    def _find_index_by_value(arr, target, tol: float = 1e-9) -> int:
        """
        CHANGE (fix): robust index lookup for float year axes (supports 2020.0 and 2020.5).
        """
        a = np.asarray(arr, dtype=float)
        hits = np.where(np.isclose(a, float(target), atol=tol, rtol=0.0))[0]
        if hits.size == 0:
            raise ValueError(f"Target {target} not found in axis (min={a.min()}, max={a.max()}).")
        return int(hits[0])


    def _init_pIRF_cache(self):
        """
        CHANGE (new): one-time cache for indices, baseline slices, and a reusable perturbed FaIR object.
        This runs once per majorghg_get_f instance (i.e., per SSP + MY).
        """
        if getattr(self, "_pirf_ready", False):
            return

        # ---- pulse year and indices ----
        pulse_year = self._get_pulse_year()

        # concentration axis (timebounds): integer years
        if not hasattr(self.f, "timebounds"):
            raise AttributeError("FaIR object missing timebounds; cannot compute pIRF safely.")
        tb = np.asarray(self.f.timebounds, dtype=float)
        tb0 = self._find_index_by_value(tb, float(pulse_year))

        # emissions axis (timepoints): year + 0.5
        if not hasattr(self.f, "timepoints"):
            raise AttributeError("FaIR object missing timepoints; cannot apply pulse on emissions axis.")
        tp = np.asarray(self.f.timepoints, dtype=float)
        tp0 = self._find_index_by_value(tp, float(pulse_year) + 0.5)

        # ---- horizon slice on timebounds ----
        nH = self.H_max * self.ts_per_year + 1
        tb1 = tb0 + nH
        if tb1 > len(tb):
            raise ValueError(
                f"Horizon exceeds available timebounds: tb0={tb0}, nH={nH}, tb1={tb1}, available={len(tb)}"
            )
        tb_slice = slice(tb0, tb1)

        # ---- cache baseline concentration slices (per gas) ----
        base_conc_h = {}
        for specie in ("CO2", "CH4", "N2O"):
            c = np.asarray(self.f.concentration.loc[dict(scenario=self.scn, specie=specie)].values)
            base_conc_h[specie] = c[tb_slice, ...]  # (H, ensemble/config)

        # ---- cache baseline emissions arrays (per gas) for fast restore ----
        base_emis = {}
        for specie in ("CO2", "CH4", "N2O"):
            e = np.asarray(self.f.emissions.loc[dict(scenario=self.scn, specie=specie)].values)
            base_emis[specie] = e.copy()

        # ---- reusable perturbed FaIR object (deepcopy ONCE per SSP+MY) ----
        f_pert = copy.deepcopy(self.f)

        # store
        self._pulse_year = pulse_year
        self._tb0 = tb0
        self._tp0 = tp0
        self._tb_slice = tb_slice
        self._nH = nH
        self._base_conc_h = base_conc_h
        self._base_emis = base_emis
        self._f_pert = f_pert
        self._pirf_ready = True


    def _restore_baseline_emissions(self, specie: str):
        """
        CHANGE (new): restore baseline emissions for one specie on the reusable f_pert.
        """
        base = self._base_emis[specie]
        # assign back into the xarray slice
        em_slice = self._f_pert.emissions.loc[dict(scenario=self.scn, specie=specie)].copy()
        em_slice.values = base
        self._f_pert.emissions.loc[dict(scenario=self.scn, specie=specie)] = em_slice


    def get_gas_deltaC_t(self, specie: str, pulse_size: float = 1.0, store: bool = True):
        """
        CHANGE (fast): rerun FaIR with an emissions pulse at pulse_year, return ΔC(t)=C_pert-C_base
        over horizon (H_max*ts_per_year+1), evaluated on timebounds.

        Returns
        -------
        dC_h : np.ndarray  shape (H, ensemble/config)
        """
        self._init_pIRF_cache()

        if specie not in ("CO2", "CH4", "N2O"):
            raise ValueError("specie must be one of: 'CO2', 'CH4', 'N2O'")

        # 1) restore baseline emissions for this specie
        self._restore_baseline_emissions(specie)

        # 2) apply pulse on timepoints index
        tp0 = self._tp0
        em_slice = self._f_pert.emissions.loc[dict(scenario=self.scn, specie=specie)].copy()
        em_vals = np.asarray(em_slice.values)

        if em_vals.ndim == 1:
            em_vals[tp0] += pulse_size
        else:
            em_vals[tp0, ...] += pulse_size

        em_slice.values = em_vals
        self._f_pert.emissions.loc[dict(scenario=self.scn, specie=specie)] = em_slice

        # 3) run perturbed FaIR
        self._f_pert.run()

        # 4) read concentration response on timebounds and compute ΔC only on horizon slice
        conc_pert = np.asarray(
            self._f_pert.concentration.loc[dict(scenario=self.scn, specie=specie)].values
        )
        conc_pert_h = conc_pert[self._tb_slice, ...]
        dC_h = conc_pert_h - self._base_conc_h[specie]

        # optional store for later inspection/export
        if store:
            if specie == "CO2":
                self.co2_f["CO2_dC_ppm_t"] = dC_h
            elif specie == "CH4":
                self.ch4_f["CH4_dC_ppb_t"] = dC_h
            else:  # N2O
                self.n2o_f["N2O_dC_ppb_t"] = dC_h

        return dC_h


    def get_gas_IRF_fair_t(
        self,
        specie: str,
        pulse_size: float = 1.0,
        store: bool = True,
        eps: float = 1e-30,
    ):
        """
        IRF(H) = ΔC(H) / ΔC_init, evaluated on timebounds.
        ΔC_init = first fully post-pulse concentration response.
        """
        # 1) compute ΔC
        dC_h = self.get_gas_deltaC_t(
            specie=specie,
            pulse_size=pulse_size,
            store=store,
        )
        # 2) use first post-pulse ΔC for normalization
        # pulse applied at timepoints = pulse_year + 0.5
        # first full concentration response on timebounds is index 1
        init_idx = 1 if dC_h.shape[0] > 1 else 0
        
        dC_init = np.asarray(dC_h[init_idx, ...]) # this should be run by the model 
        
        # then detect problematic normalization values, eps = e-30
        mask_small = np.abs(dC_init) < eps

        if np.any(mask_small):
            n_bad = np.sum(mask_small)
            print(
                f"[WARNING] {specie} pIRF normalization: "
                f"{n_bad} ensemble member(s) have |ΔC_init| < {eps:.1e}. "
                f"Replacing with eps to avoid division instability."
            )
        # safe denominator
        dC_init_safe = np.where(mask_small, eps, dC_init)

        irf_h = dC_h / dC_init_safe

        # enforce IRF start = 1 (convention)
        irf_h[0, ...] = 1.0

        # 3) store
        if store:
            if specie == "CO2":
                self.co2_f["CO2_IRF_fair_t"] = irf_h
            elif specie == "CH4":
                self.ch4_f["CH4_IRF_fair_t"] = irf_h
            else:
                self.n2o_f["N2O_IRF_fair_t"] = irf_h

        return irf_h

    # convenience wrappers (unchanged API)
    def get_co2_pIRF(self, pulse_size: float = 1.0):
        return self.get_gas_IRF_fair_t("CO2", pulse_size=pulse_size, store=True)

    def get_ch4_pIRF(self, pulse_size: float = 1.0):
        return self.get_gas_IRF_fair_t("CH4", pulse_size=pulse_size, store=True)

    def get_n2o_pIRF(self, pulse_size: float = 1.0):
        return self.get_gas_IRF_fair_t("N2O", pulse_size=pulse_size, store=True)