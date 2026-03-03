import os
import copy
from typing import Optional, Union

import numpy as np
import pandas as pd


class majorghg_get_f:
    """Module A: extract FaIR parameters and compute climate-metric inputs.

    Key outputs
    -----------
    - Concentrations (C), baseline concentration (RF), forcing_scale (FS)
    - Scenario-conditioned lifetimes (alpha) for CH4/N2O (CO2 handled elsewhere)
    - Prospective impulse response functions (pIRF) computed via a *perturbed FaIR rerun*

    Important FaIR-species handling (CO2)
    ------------------------------------
    Many FaIR setups expose *three* CO2-related species in the coords:
      - emissions:  'CO2 FFI', 'CO2 AFOLU' (and sometimes 'CO2')
      - concentration: 'CO2'

    For pIRF, the concentration response must be read from 'CO2', but the pulse must be
    applied to the *emissions* species that FaIR actually uses to drive the carbon cycle.
    By default we pulse 'CO2 FFI' if present; otherwise we fall back to 'CO2'.
    """

    def __init__(
        self,
        f,
        scn,
        H_max: int = 100,
        ts_per_year: int = 1,
        scenarios=None,
        year_index: int = 269,
        co2_emission_specie: Optional[str] = None,
        debug: bool = True,
    ):
        self.f = f
        self.scn = scn
        self.H_max = int(H_max)
        self.year_index = int(year_index)
        self.ts_per_year = int(ts_per_year)
        self.H = np.linspace(0, self.H_max, self.H_max * self.ts_per_year + 1)

        self.debug = bool(debug)

        # Optional override: let caller force which emissions specie to pulse for CO2.
        self.co2_emission_specie = co2_emission_specie

        # core FaIR-derived parameter dicts
        self.co2_f, self.ch4_f, self.n2o_f = self.call_f_from_fair_gas()
        self.n2o_c = self.co2_f["N2O_C"]

    # =========================
    # step A.1: call from fair, get parameter needed for climate metric calculation
    # =========================
    def call_f_from_fair_gas(self):
        # Reserve keys for prospective DeltaC/IRF outputs.
        co2_key = [
            "CO2_C", "CO2_C_plus1", "CO2_RF", "CO2_FS", "N2O_C",
            "CO2_dC_ppm_t", "CO2_IRF_fair_t",
        ]
        ch4_key = [
            "CH4_C", "CH4_C_plus1", "CH4_RF", "CH4_FS", "N2O_C", "CH4_lifetime",
            "CH4_dC_ppb_t", "CH4_IRF_fair_t",
        ]
        n2o_key = [
            "N2O_C", "N2O_C_plus1", "N2O_RF", "N2O_FS", "N2O_lifetime",
            "N2O_dC_ppb_t", "N2O_IRF_fair_t",
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

            n2o_c = concentration.loc[dict(specie="N2O")].values

            if gas == "CO2":
                co2_f["CO2_C"] = gas_c
                co2_f["CO2_C_plus1"] = gas_c_plus1
                co2_f["CO2_RF"] = gas_rf
                co2_f["CO2_FS"] = gas_fs
                co2_f["N2O_C"] = n2o_c

            elif gas == "CH4":
                alpha_s = self.f.alpha_lifetime.loc[dict(scenario=self.scn, specie="CH4")][self.year_index].values
                ch4_unp = 10.8537568
                ch4_f["CH4_lifetime"] = ch4_unp * alpha_s

                ch4_f["CH4_C"] = gas_c
                ch4_f["CH4_C_plus1"] = gas_c_plus1
                ch4_f["CH4_RF"] = gas_rf
                ch4_f["CH4_FS"] = gas_fs
                ch4_f["N2O_C"] = n2o_c

            elif gas == "N2O":
                alpha_s = self.f.alpha_lifetime.loc[dict(scenario=self.scn, specie="N2O")][self.year_index].values
                n2o_unp = 109.0
                n2o_f["N2O_lifetime"] = n2o_unp * alpha_s

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
        """Meinshausen (2020)-style ERF for N2O (FaIR implementation).

        Inputs
        ------
        n2o_c_orcplus : array
            N2O concentration array (timebounds x ensemble), baseline or +1 ppb.
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
    # step A.3: ERF per unit concentration perturbation (RE)
    # =========================
    def get_majorghg_RE(self):
        # CO2
        erf_co2 = self.get_co2_meinshausen2020(self.co2_f["CO2_C"])
        erf_co2_plus1 = self.get_co2_meinshausen2020(self.co2_f["CO2_C_plus1"])
        co2_re = erf_co2_plus1 - erf_co2

        # CH4
        erf_ch4 = self.get_ch4_meinshausen2020(self.ch4_f["CH4_C"])
        erf_ch4_plus1 = self.get_ch4_meinshausen2020(self.ch4_f["CH4_C_plus1"])
        ch4_re = erf_ch4_plus1 - erf_ch4

        # N2O
        erf_n2o = self.get_n2o_meinshausen2020(self.n2o_f["N2O_C"])
        erf_n2o_plus1 = self.get_n2o_meinshausen2020(self.n2o_f["N2O_C_plus1"])
        n2o_re = erf_n2o_plus1 - erf_n2o

        return co2_re, ch4_re, n2o_re

    # =========================
    # step A.4: alpha_lifetime for CH4 and N2O
    # =========================
    def get_ch4_n2o_alpha(self):
        return self.ch4_f["CH4_lifetime"], self.n2o_f["N2O_lifetime"]

    def get_majorghg_concentration(self):
        return self.co2_f["CO2_C"], self.ch4_f["CH4_C"], self.n2o_f["N2O_C"]

    # ======================================================================
    # CHANGE: prospective decay side (`pIRF`) via perturbed FaIR rerun
    # ======================================================================
    def _get_pulse_year(self) -> int:
        if hasattr(self.f, "timebounds"):
            tb0 = float(np.asarray(self.f.timebounds, dtype=float)[0])
            Y0 = int(round(tb0))
            return int(Y0 + int(self.year_index))
        return int(self.year_index)

    @staticmethod
    def _find_index_by_value(arr, target, tol: float = 1e-9) -> int:
        a = np.asarray(arr, dtype=float)
        hits = np.where(np.isclose(a, float(target), atol=tol, rtol=0.0))[0]
        if hits.size == 0:
            raise ValueError(f"Target {target} not found in axis (min={a.min()}, max={a.max()}).")
        return int(hits[0])

    def _resolve_emission_specie(self, specie: str) -> str:
        """Return the *emissions* specie name to pulse for a given logical gas."""
        if specie != "CO2":
            return specie

        # User override wins if it exists in the FaIR object.
        if self.co2_emission_specie is not None:
            if self.co2_emission_specie in set(self.f.emissions.coords["specie"].values.tolist()):
                return self.co2_emission_specie
            raise ValueError(
                f"co2_emission_specie='{self.co2_emission_specie}' not found in f.emissions specie coords."
            )

        emis_species = set(self.f.emissions.coords["specie"].values.tolist())
        if "CO2 FFI" in emis_species:
            return "CO2 FFI"
        if "CO2" in emis_species:
            return "CO2"
        # Last resort: try AFOLU if that's all there is.
        if "CO2 AFOLU" in emis_species:
            return "CO2 AFOLU"
        raise ValueError("No CO2 emissions specie found among {'CO2 FFI','CO2','CO2 AFOLU'}.")

    def get_gas_deltaC_t(self, specie: str, pulse_size: float = 1.0, store: bool = True):
        """Compute prospective concentration response ΔC(t) over the moving horizon.

        Notes
        -----
        - Concentration is always read from the *concentration* species:
            'CO2', 'CH4', 'N2O'
        - The pulse is applied to the resolved *emissions* species:
            CO2 -> typically 'CO2 FFI' (or fallback)
            CH4/N2O -> same name
        - H=0 corresponds to the *first realized* response at (pulse_year + 1) on timebounds.
        """
        if specie not in ("CO2", "CH4", "N2O"):
            raise ValueError("specie must be one of: 'CO2', 'CH4', 'N2O'")

        pulse_year = self._get_pulse_year()

        if not hasattr(self.f, "timebounds") or not hasattr(self.f, "timepoints"):
            raise AttributeError("FaIR object must provide both timebounds (concentration) and timepoints (emissions).")

        tb = np.asarray(self.f.timebounds, dtype=float)
        tp = np.asarray(self.f.timepoints, dtype=float)
        tb_pulse = self._find_index_by_value(tb, float(pulse_year))
        tp_pulse = self._find_index_by_value(tp, float(pulse_year) + 0.5)
        tb_resp0 = tb_pulse + 1

        conc_base = np.asarray(self.f.concentration.loc[dict(scenario=self.scn, specie=specie)].values)

        f_pert = copy.deepcopy(self.f)

        em_specie = self._resolve_emission_specie(specie)

        em_slice = f_pert.emissions.loc[dict(scenario=self.scn, specie=em_specie)].copy()
        em_vals = np.asarray(em_slice.values)

        if self.debug:
            em_before = em_vals.copy()

        if em_vals.ndim == 1:
            base_em_med = float(np.median(em_vals[tp_pulse]))
            em_vals[tp_pulse] = em_vals[tp_pulse] + pulse_size
            em_after_med = float(np.median(em_vals[tp_pulse]))
        else:
            base_em_med = float(np.median(em_vals[tp_pulse, ...]))
            em_vals[tp_pulse, ...] = em_vals[tp_pulse, ...] + pulse_size
            em_after_med = float(np.median(em_vals[tp_pulse, ...]))

        em_slice.values = em_vals
        f_pert.emissions.loc[dict(scenario=self.scn, specie=em_specie)] = em_slice

        # Run perturbed model.
        f_pert.run()

        conc_pert = np.asarray(f_pert.concentration.loc[dict(scenario=self.scn, specie=specie)].values)

        if self.debug:
            conc_diff_max = float(np.nanmax(np.abs(conc_pert - conc_base)))
            print(f"[SANITY] {specie} pulse applied to emissions specie='{em_specie}'.")
            print(f"[SANITY] {specie} emis@tp_pulse median before/after: {base_em_med:.6g} -> {em_after_med:.6g}")
            print(f"[SANITY] {specie} max |conc_after - conc_before| = {conc_diff_max:.6g}")

        dC = conc_pert - conc_base

        nH = self.H_max * self.ts_per_year + 1
        tb1 = tb_resp0 + nH
        if tb1 > dC.shape[0]:
            raise ValueError(
                "Horizon exceeds available timebounds. "
                f"tb_resp0={tb_resp0}, nH={nH}, tb1={tb1}, available={dC.shape[0]}"
            )

        dC_h = dC[tb_resp0:tb1, ...]

        if self.debug:
            dC0_med = float(np.nanmedian(dC_h[0, ...]))
            dCend_med = float(np.nanmedian(dC_h[-1, ...]))
            print(
                f"[pIRF-debug] scn={self.scn} year_index={self.year_index} pulse_year={pulse_year} gas={specie} "
                f"tb_pulse={tb_pulse}(year={tb[tb_pulse]:.1f}) "
                f"tp_pulse={tp_pulse}(year={tp[tp_pulse]:.1f}) "
                f"tb_resp0={tb_resp0}(year={tb[tb_resp0]:.1f}) "
                f"ΔC0~{dC0_med:.3g} ΔC_end(H={nH-1})~{dCend_med:.3g}"
            )

        if store:
            if specie == "CO2":
                self.co2_f["CO2_dC_ppm_t"] = dC_h
            elif specie == "CH4":
                self.ch4_f["CH4_dC_ppb_t"] = dC_h
            else:
                self.n2o_f["N2O_dC_ppb_t"] = dC_h

        return dC_h

    def get_gas_IRF_fair_t(
        self,
        specie: str,
        pulse_size: float = 1.0,
        store: bool = True,
        eps: float = 1e-30,
    ):
        """Prospective impulse response function from FaIR reruns: IRF(t) = ΔC(t)/ΔC(t0)."""
        dC_h = self.get_gas_deltaC_t(specie=specie, pulse_size=pulse_size, store=store)

        dC0 = np.asarray(dC_h[0, ...])
        mask_small = np.abs(dC0) < eps
        if np.any(mask_small):
            n_bad = int(np.sum(mask_small))
            print(
                f"[WARNING] {specie} pIRF: {n_bad} ensemble member(s) have |ΔC0| < {eps:.1e} "
                f"(scn={self.scn}, MY={self._get_pulse_year()}). Replacing with eps."
            )
        dC0_safe = np.where(mask_small, eps, dC0)

        irf_h = dC_h / dC0_safe
        irf_h[0, ...] = 1.0

        if store:
            if specie == "CO2":
                self.co2_f["CO2_IRF_fair_t"] = irf_h
            elif specie == "CH4":
                self.ch4_f["CH4_IRF_fair_t"] = irf_h
            else:
                self.n2o_f["N2O_IRF_fair_t"] = irf_h

        return irf_h

    def get_co2_pIRF(self, pulse_size: float = 1.0):
        return self.get_gas_IRF_fair_t("CO2", pulse_size=pulse_size, store=True)

    def get_ch4_pIRF(self, pulse_size: float = 1.0):
        return self.get_gas_IRF_fair_t("CH4", pulse_size=pulse_size, store=True)

    def get_n2o_pIRF(self, pulse_size: float = 1.0):
        return self.get_gas_IRF_fair_t("N2O", pulse_size=pulse_size, store=True)
