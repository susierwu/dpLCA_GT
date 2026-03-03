# moduleA.py
# =============================================================================
# Module A: extract FaIR outputs for major GHGs and compute:
#   (i) concentration-dependent ERF deltas (Meinshausen2020; FaIR implementation)
#   (ii) prospective impulse response functions (pIRF) via a perturbed FaIR rerun
#
# Notes on pIRF:
#   - FaIR is emissions-driven in the core run; pIRF is obtained by an additional
#     *emissions pulse* experiment and differencing concentrations:
#         ΔC(t) = C_pert(t) - C_base(t)
#         IRF(t) = ΔC(t) / ΔC(t0)
#     where t0 is the first post-pulse concentration timebound (annual grid: Y+1).
#   - For CO2, FaIR commonly uses emissions species "CO2 FFI" and/or "CO2 AFOLU"
#     to drive atmospheric CO2, while concentration is reported as "CO2".
#     Therefore, pIRF must allow decoupled (specie_emis, specie_conc).
# =============================================================================

from __future__ import annotations

from typing import Dict, Tuple, Optional, Union, Sequence
import copy
import numpy as np


def _as_float_array(x) -> np.ndarray:
    return np.asarray(x, dtype=float)


def _find_index_near(arr: Union[np.ndarray, Sequence[float]], target: float, tol: float = 1e-6) -> int:
    """
    Find index i where arr[i] is closest to target within absolute tolerance.
    Works for float axes like FaIR timepoints (1750.5, 1751.5, ...).
    """
    a = _as_float_array(arr)
    i = int(np.argmin(np.abs(a - float(target))))
    if abs(a[i] - float(target)) > tol:
        raise ValueError(f"Target={target} not found within tol={tol}. Closest={a[i]} at i={i}.")
    return i


def _find_index_exact(arr: Union[np.ndarray, Sequence[float]], value: float, tol: float = 1e-9) -> int:
    """
    Exact-ish match for float axes (supports 2030.0, 2030.5).
    """
    a = _as_float_array(arr)
    hits = np.where(np.isclose(a, float(value), atol=tol, rtol=0.0))[0]
    if hits.size == 0:
        raise ValueError(f"Value={value} not found in axis (min={a.min()}, max={a.max()}).")
    return int(hits[0])


def _get_dt_years(ts_per_year: int) -> float:
    return 1.0 / float(ts_per_year)


class majorghg_get_f:
    """
    In order to calculate climate metrics, this module A runs:
      1) extract major GHG FaIR output parameters (C, baseline, forcing_scale, lifetime scale)
      2) compute ERF using Meinshausen2020 parameterization (FaIR forcing/ghg.py)
      3) compute ERF per unit concentration perturbation ([W m-2 ppm-1] CO2, [W m-2 ppb-1] CH4/N2O)
      4) compute prospective impulse response functions pIRF via a perturbed FaIR re-run:
         DeltaC(t) = C_pert(t) - C_base(t), IRF(t) = DeltaC(t) / DeltaC(t0)

    Parameters
    ----------
    f : FaIR object
        A configured FaIR instance with emissions populated and baseline run already executed (f.run()).
    scn : str
        Scenario label matching f's scenario coordinate.
    H_max : int
        Horizon (years) for pIRF extraction (H=0..H_max).
    ts_per_year : int
        Timesteps per year in FaIR run (annual=1). This implementation assumes annual (1) unless
        you are sure your timepoints/timebounds follow the same conventions.
    year_index : int
        Index on timebounds corresponding to model year (MY - FAIR_START_Y) for annual runs.
    debug : bool
        If True, prints detailed sanity/debug lines during pIRF runs.
    """

    def __init__(
        self,
        f,
        scn: str,
        H_max: int = 100,
        ts_per_year: int = 1,
        scenarios: Optional[Sequence[str]] = None,
        year_index: int = 269,
        debug: bool = False,
    ):
        self.f = f
        self.scn = scn
        self.H_max = int(H_max)
        self.year_index = int(year_index)
        self.ts_per_year = int(ts_per_year)
        self.debug = bool(debug)

        self.H = np.linspace(0, self.H_max, self.H_max * self.ts_per_year + 1)

        # core FaIR-derived parameter dicts
        self.co2_f, self.ch4_f, self.n2o_f = self.call_f_from_fair_gas()
        self.n2o_c = self.n2o_f["N2O_C"]  # retained for CO2 Meinshausen coupling

        # cache pulse year for debug printing (derived from timebounds)
        self._pulse_year_cache: Optional[int] = None

    # =========================
    # A.1: extract from FaIR
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
            Also reserves keys for prospective DeltaC/IRF (pIRF) outputs.
        """
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
            gas_c_plus1 = gas_c + 1.0

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
    # A.2: Meinshausen2020 ERF
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
        a2=-0.00034197, b2=0.00025455, c2=-0.00024357, d2=0.12173,
    ):
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
    # A.3: ERF per 1 unit concentration perturbation
    # =========================
    def get_co2_1ppm_erf(self):
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
    # A.4: lifetimes & concentrations
    # =========================
    def get_ch4_n2o_alpha(self):
        return self.ch4_f["CH4_lifetime"], self.n2o_f["N2O_lifetime"]

    def get_majorghg_concentration(self):
        return self.co2_f["CO2_C"], self.ch4_f["CH4_C"], self.n2o_f["N2O_C"]

    # =============================================================================
    # pIRF utilities
    # =============================================================================
    def _get_pulse_year(self) -> int:
        """
        Pulse year is derived from timebounds at year_index (calendar year for MY).
        For annual runs: timebounds[k] = year, timepoints[k] = year + 0.5.
        """
        if self._pulse_year_cache is not None:
            return self._pulse_year_cache
        if not hasattr(self.f, "timebounds"):
            # fallback: interpret year_index as a calendar year
            self._pulse_year_cache = int(self.year_index)
            return self._pulse_year_cache
        tb = _as_float_array(self.f.timebounds)
        self._pulse_year_cache = int(tb[int(self.year_index)])
        return self._pulse_year_cache

    def _get_timebound_index(self, pulse_year: int) -> int:
        if not hasattr(self.f, "timebounds"):
            return int(self.year_index)
        return _find_index_exact(self.f.timebounds, float(pulse_year), tol=1e-9)

    def _get_timepoint_index(self, pulse_year: int) -> int:
        """
        For annual TS_PER_YEAR=1: timepoints are 1750.5, 1751.5, ...
        A pulse in calendar year Y is applied at timepoint (Y + 0.5*dt).
        """
        if not hasattr(self.f, "timepoints"):
            # fallback: approximate from timebounds
            tb0 = self._get_timebound_index(pulse_year)
            return int(max(tb0 - 1, 0))
        dt = _get_dt_years(self.ts_per_year)
        target_tp = float(pulse_year) + 0.5 * dt
        return _find_index_near(self.f.timepoints, target_tp, tol=1e-4)

    @staticmethod
    def _normalize_pulse_split(pulse_size: float, n: int, pulse_split: str) -> np.ndarray:
        """
        Split a scalar pulse_size across n emissions channels.
        pulse_split:
          - "equal": each gets pulse_size/n so total added = pulse_size
          - "none": each gets pulse_size (total added = n*pulse_size)
        """
        if n <= 0:
            raise ValueError("n must be positive.")
        if pulse_split not in ("equal", "none"):
            raise ValueError("pulse_split must be one of: 'equal', 'none'")
        if pulse_split == "equal":
            return np.full(n, float(pulse_size) / float(n))
        return np.full(n, float(pulse_size))

    def get_gas_deltaC_t(
        self,
        specie_conc: str,
        specie_emis: Union[str, Tuple[str, ...]],
        pulse_size: float = 1.0,
        store: bool = True,
        pulse_split: str = "equal",
    ) -> np.ndarray:
        """
        Rerun FaIR with a one-timestep emissions pulse at pulse_year; return ΔC on timebounds
        for horizons H=0..H_max (length nH = H_max*ts_per_year+1).

        The concentration response is measured for specie_conc on f.concentration.
        The pulse is applied to specie_emis on f.emissions (supports CO2 split channels).

        H=0 is defined at the first realized post-pulse concentration timebound:
            tb_resp0 = tb_pulse + 1
        for annual grids where emissions are at mid-year timepoints.
        """
        # ---- indices / calendar mapping ----
        pulse_year = self._get_pulse_year()
        tb_pulse = self._get_timebound_index(pulse_year)
        tp_pulse = self._get_timepoint_index(pulse_year)
        tb_resp0 = tb_pulse + 1  # first post-pulse timebound on annual grids

        if not hasattr(self.f, "timebounds") or not hasattr(self.f, "timepoints"):
            raise AttributeError("FaIR object must expose both timebounds and timepoints for pIRF.")

        tb_axis = _as_float_array(self.f.timebounds)
        tp_axis = _as_float_array(self.f.timepoints)

        # ---- clone base state (critical for correctness) ----
        f_pert = copy.deepcopy(self.f)

        # ---- baseline concentration from the SAME object we will perturb ----
        conc_before = np.asarray(
            f_pert.concentration.loc[dict(scenario=self.scn, specie=specie_conc)].values
        ).copy()

        # ---- apply pulse(s) on emissions axis ----
        if isinstance(specie_emis, str):
            emis_species = (specie_emis,)
        else:
            emis_species = tuple(specie_emis)

        pulse_parts = self._normalize_pulse_split(pulse_size=float(pulse_size), n=len(emis_species), pulse_split=pulse_split)

        base_em_before_list = []
        base_em_after_list = []

        for sp_e, p_add in zip(emis_species, pulse_parts):
            em_slice = f_pert.emissions.loc[dict(scenario=self.scn, specie=sp_e)].copy()
            em_vals = np.asarray(em_slice.values)

            if em_vals.ndim == 1:
                base_before = float(em_vals[tp_pulse])
                em_vals[tp_pulse] += p_add
                base_after = float(em_vals[tp_pulse])
            else:
                base_before = float(np.median(em_vals[tp_pulse, ...]))
                em_vals[tp_pulse, ...] += p_add
                base_after = float(np.median(em_vals[tp_pulse, ...]))

            em_slice.values = em_vals
            f_pert.emissions.loc[dict(scenario=self.scn, specie=sp_e)] = em_slice

            base_em_before_list.append((sp_e, base_before))
            base_em_after_list.append((sp_e, base_after))

        # ---- run perturbed model ----
        f_pert.run()

        # ---- perturbed concentration from the SAME object ----
        conc_after = np.asarray(
            f_pert.concentration.loc[dict(scenario=self.scn, specie=specie_conc)].values
        ).copy()

        # ---- Delta concentration defined consistently ----
        dC = conc_after - conc_before

        # ---- horizon slice ----
        nH = self.H_max * self.ts_per_year + 1
        tb1 = tb_resp0 + nH
        if tb1 > dC.shape[0]:
            raise ValueError(
                f"Horizon exceeds available timebounds. tb_resp0={tb_resp0}, nH={nH}, tb1={tb1}, available={dC.shape[0]}"
            )
        dC_h = dC[tb_resp0:tb1, ...]

        # ---- debug prints (consistent with dC_h) ----
        if self.debug:
            # emissions sanity
            for (sp_e, b0), (_, b1) in zip(base_em_before_list, base_em_after_list):
                print(f"[SANITY] {sp_e} emis@tp_pulse median before/after: {b0:.6g} -> {b1:.6g}")

            # report calendar years
            print(
                f"[pIRF-debug] scn={self.scn} year_index={self.year_index} pulse_year={pulse_year} "
                f"gas={specie_conc} tb_pulse={tb_pulse}(year={tb_axis[tb_pulse]:.1f}) "
                f"tp_pulse={tp_pulse}(year={tp_axis[tp_pulse]:.1f}) "
                f"tb_resp0={tb_resp0}(year={tb_axis[tb_resp0]:.1f}) "
                f"base_em(tp)~{base_em_before_list[0][1]:.6g}"
            )

            dC0_med = float(np.nanmedian(dC_h[0, ...]))
            dCend_med = float(np.nanmedian(dC_h[-1, ...]))
            max_abs_dC_full = float(np.nanmax(np.abs(dC)))
            max_abs_dC_slice = float(np.nanmax(np.abs(dC_h)))
            print(f"[SANITY] {specie_conc} max |ΔC| full={max_abs_dC_full:.6g} slice={max_abs_dC_slice:.6g}")
            print(f"[SANITY] {specie_conc} ΔC_h median: ΔC0={dC0_med:.6g} ΔCend(H={self.H_max})={dCend_med:.6g}")

        # ---- optional store ----
        if store:
            if specie_conc == "CO2":
                self.co2_f["CO2_dC_ppm_t"] = dC_h
            elif specie_conc == "CH4":
                self.ch4_f["CH4_dC_ppb_t"] = dC_h
            elif specie_conc == "N2O":
                self.n2o_f["N2O_dC_ppb_t"] = dC_h

        return dC_h

    def get_gas_IRF_fair_t(
        self,
        specie_conc: str,
        specie_emis: Union[str, Tuple[str, ...]],
        pulse_size: float = 1.0,
        store: bool = True,
        eps: float = 1e-30,
        pulse_split: str = "equal",
    ) -> np.ndarray:
        """
        IRF(H) = ΔC(H)/ΔC(0), where H=0 is defined at the first post-pulse timebound (tb_pulse+1).
        """
        dC_h = self.get_gas_deltaC_t(
            specie_conc=specie_conc,
            specie_emis=specie_emis,
            pulse_size=pulse_size,
            store=store,
            pulse_split=pulse_split,
        )

        dC0 = np.asarray(dC_h[0, ...])
        mask_small = np.abs(dC0) < eps
        if np.any(mask_small) and self.debug:
            n_bad = int(np.sum(mask_small))
            print(
                f"[WARNING] {specie_conc} pIRF: {n_bad} ensemble member(s) have |ΔC0| < {eps:.1e} "
                f"(scn={self.scn}, pulse_year={self._get_pulse_year()}). Replacing with eps."
            )
        dC0_safe = np.where(mask_small, eps, dC0)

        irf_h = dC_h / dC0_safe

        if store:
            if specie_conc == "CO2":
                self.co2_f["CO2_IRF_fair_t"] = irf_h
            elif specie_conc == "CH4":
                self.ch4_f["CH4_IRF_fair_t"] = irf_h
            elif specie_conc == "N2O":
                self.n2o_f["N2O_IRF_fair_t"] = irf_h

        return irf_h

    # =============================================================================
    # Convenience wrappers
    # =============================================================================
    def get_co2_pIRF(
        self,
        pulse_size: float = 1.0,
        co2_emis: Tuple[str, ...] = ("CO2 FFI",),
        pulse_split: str = "equal",
        store: bool = True,
    ) -> np.ndarray:
        """
        Default: pulse emissions on ("CO2 FFI",) and measure concentration response on "CO2".
        If you want total anthropogenic CO2, pass co2_emis=("CO2 FFI","CO2 AFOLU") and pulse_split="equal".
        """
        return self.get_gas_IRF_fair_t(
            specie_conc="CO2",
            specie_emis=co2_emis,
            pulse_size=pulse_size,
            store=store,
            pulse_split=pulse_split,
        )

    def get_ch4_pIRF(self, pulse_size: float = 1.0, store: bool = True) -> np.ndarray:
        return self.get_gas_IRF_fair_t(
            specie_conc="CH4",
            specie_emis="CH4",
            pulse_size=pulse_size,
            store=store,
        )

    def get_n2o_pIRF(self, pulse_size: float = 1.0, store: bool = True) -> np.ndarray:
        return self.get_gas_IRF_fair_t(
            specie_conc="N2O",
            specie_emis="N2O",
            pulse_size=pulse_size,
            store=store,
        )
