from __future__ import annotations
import copy
from typing import Dict, Tuple, Union
import numpy as np


def _as_float_array(x) -> np.ndarray:
    return np.asarray(x, dtype=float)


def _find_index_near(arr: Union[np.ndarray, list], target: float, tol: float = 1e-6) -> int:
    a = _as_float_array(arr)
    i = int(np.argmin(np.abs(a - float(target))))
    if abs(a[i] - float(target)) > tol:
        raise ValueError(f"Target={target} not found within tol={tol}. Closest={a[i]} at i={i}.")
    return i


def _find_index_exact_year(arr: Union[np.ndarray, list], year_int: int) -> int:
    a = _as_float_array(arr)
    hits = np.where(a == float(int(year_int)))[0]
    if hits.size == 0:
        raise ValueError(f"Year {year_int} not found in axis.")
    return int(hits[0])


def _first_index_ge(arr: Union[np.ndarray, list], target: float) -> int:
    a = _as_float_array(arr)
    hits = np.where(a >= float(target))[0]
    if hits.size == 0:
        raise ValueError(f"Target {target} beyond axis max {a.max()}.")
    return int(hits[0])


class majorghg_get_f:
    """
    FaIR pIRF calculator
    Uses full deepcopy for both baseline and perturbed to prevent FaIR state bleed-through.
    """

    def __init__(
        self,
        f,
        scn: str,
        H_max: int = 100,
        ts_per_year: int = 1,
        year_index: int = 269,
        debug: bool = False,
    ):
        self.f = f
        self.scn = str(scn)
        self.H_max = int(H_max)
        self.ts_per_year = int(ts_per_year)
        self.year_index = int(year_index)
        self.debug = bool(debug)

        self.H = np.linspace(0, self.H_max, self.H_max * self.ts_per_year + 1)

        # Optional storage
        self.co2_f = {"CO2_dC_ppm_t": None, "CO2_IRF_fair_t": None}
        self.ch4_f = {"CH4_dC_ppb_t": None, "CH4_IRF_fair_t": None}
        self.n2o_f = {"N2O_dC_ppb_t": None, "N2O_IRF_fair_t": None}

    # -------------------------
    # time helpers
    # -------------------------
    def _get_pulse_year(self) -> int:
        tb = _as_float_array(self.f.timebounds)
        return int(tb[int(self.year_index)])

    def _get_timebound_index(self, pulse_year: int) -> int:
        return _find_index_exact_year(self.f.timebounds, pulse_year)

    def _get_timepoint_index(self, pulse_year: int) -> int:
        dt = 1.0 / float(self.ts_per_year)
        target_tp = float(pulse_year) + 0.5 * dt
        return _find_index_near(self.f.timepoints, target_tp, tol=1e-4)

    # -------------------------
    # core ΔC and IRF
    # -------------------------
    def get_gas_deltaC_t(
        self,
        specie_conc: str,
        specie_emis: Union[str, Tuple[str, ...]],
        pulse_size: float = 1.0,
        pulse_split: str = "equal",     # "equal" or "none"
        pre_window_years: int = 20,
        eps_prepulse: float = 1e-3,
        store: bool = True,
    ) -> np.ndarray:
        """
        Baseline/perturbed computed from two independent deepcopies of self.f:
          f_base = deepcopy(self.f); run
          f_pert = deepcopy(self.f); pulse; run
        ΔC = conc_pert - conc_base (same scenario/specie)
        Enforces pre-pulse ΔC ~ 0.
        """

        # normalize emissions species
        if isinstance(specie_emis, str):
            emis_species = (specie_emis,)
        else:
            emis_species = tuple(specie_emis)

        if pulse_split not in ("equal", "none"):
            raise ValueError("pulse_split must be 'equal' or 'none'.")

        if pulse_split == "equal":
            pulse_parts = [float(pulse_size) / float(len(emis_species))] * len(emis_species)
        else:
            pulse_parts = [float(pulse_size)] * len(emis_species)

        # indices based on self.f axes
        pulse_year = self._get_pulse_year()
        tb_pulse = self._get_timebound_index(pulse_year)
        tp_pulse = self._get_timepoint_index(pulse_year)

        dt = 1.0 / float(self.ts_per_year)
        tp_axis_self = _as_float_array(self.f.timepoints)
        tb_axis_self = _as_float_array(self.f.timebounds)
        tp_pulse_time = float(tp_axis_self[tp_pulse])

        # response starts at first timebound >= tp_pulse_time + dt
        # so “response starting one timestep after the instant physical response” 
        # but should be more precisely as:  tb_resp0 = _first_index_ge(tb_axis_self, tp_pulse_time - 1e-12)
        tb_resp0 = _first_index_ge(tb_axis_self, tp_pulse_time + dt - 1e-12)

        nH = self.H_max * self.ts_per_year + 1
        tb1 = tb_resp0 + nH
        if tb1 > tb_axis_self.size:
            raise ValueError(f"Horizon exceeds timebounds: tb_resp0={tb_resp0}, nH={nH}, tb1={tb1}, avail={tb_axis_self.size}")

        # ---- baseline run (independent) ----
        f_base = copy.deepcopy(self.f)
        f_base.run()
        conc_base = np.asarray(
            f_base.concentration.loc[dict(scenario=self.scn, specie=specie_conc)].values
        ).copy()

        # ---- perturbed run (independent) ----
        f_pert = copy.deepcopy(self.f)

        sanity_before_after = []
        for sp_e, p_add in zip(emis_species, pulse_parts):
            em_slice = f_pert.emissions.loc[dict(scenario=self.scn, specie=sp_e)].copy()
            em_vals = np.asarray(em_slice.values)

            if em_vals.ndim == 1:
                before = float(em_vals[tp_pulse])
                em_vals[tp_pulse] = before + p_add
                after = float(em_vals[tp_pulse])
            else:
                before = float(np.nanmedian(em_vals[tp_pulse, ...]))
                em_vals[tp_pulse, ...] = em_vals[tp_pulse, ...] + p_add
                after = float(np.nanmedian(em_vals[tp_pulse, ...]))

            em_slice.values = em_vals
            f_pert.emissions.loc[dict(scenario=self.scn, specie=sp_e)] = em_slice
            sanity_before_after.append((sp_e, before, after, p_add))

        if self.debug:
            tb_axis = _as_float_array(f_base.timebounds)
            tp_axis = _as_float_array(f_base.timepoints)
            print(
                f"[pIRF-base] scn={self.scn} year_index={self.year_index} pulse_year={pulse_year} "
                f"tb_pulse={tb_pulse}(tb={tb_axis[tb_pulse]:.1f}) tp_pulse={tp_pulse}(tp={tp_axis[tp_pulse]:.3f})"
            )
            print(
                f"[pIRF-debug] scn={self.scn} year_index={self.year_index} pulse_year={pulse_year} "
                f"gas={specie_conc} emis={','.join(emis_species)} "
                f"tb_pulse={tb_pulse}(tb={tb_axis[tb_pulse]:.1f}) tp_pulse={tp_pulse}(tp={tp_axis[tp_pulse]:.3f}) "
                f"dt={dt:g} tb_resp0={tb_resp0}(tb={tb_axis[tb_resp0]:.1f})"
            )
            for (sp_e, b0, b1, p_add) in sanity_before_after:
                print(f"[SANITY] {sp_e} emis@tp_pulse median before/after: {b0:.6g} -> {b1:.6g} (add {p_add:.6g})")

        f_pert.run()
        conc_pert = np.asarray(
            f_pert.concentration.loc[dict(scenario=self.scn, specie=specie_conc)].values
        ).copy()

        if conc_pert.shape != conc_base.shape:
            raise RuntimeError(f"Shape mismatch: base={conc_base.shape}, pert={conc_pert.shape}")

        dC = conc_pert - conc_base

        # pre-pulse mismatch check
        pre0 = max(tb_pulse - int(pre_window_years * self.ts_per_year), 0)
        pre_slice = dC[pre0:tb_pulse, ...]
        pre_med = float(np.nanmedian(pre_slice)) if pre_slice.size else 0.0
        pre_max = float(np.nanmax(np.abs(pre_slice))) if pre_slice.size else 0.0
        full_max = float(np.nanmax(np.abs(dC)))

        if self.debug:
            print(f"[SANITY-ΔC] pre-pulse median={pre_med:.6g}, pre-pulse maxabs={pre_max:.6g}, full maxabs={full_max:.6g}")

        if pre_max > float(eps_prepulse):
            raise RuntimeError(
                "Baseline mismatch detected: ΔC is not ~0 before pulse.\n"
                f"pre_window_years={pre_window_years}, pre_maxabs={pre_max:.6g} (thresh {eps_prepulse:.6g}), pre_median={pre_med:.6g}\n"
                "This indicates baseline/perturbed are not consistent (scenario selection/run state mismatch)."
            )

        dC_h = dC[tb_resp0:tb1, ...].copy()

        if self.debug:
            dC0_med = float(np.nanmedian(dC_h[0, ...]))
            dCend_med = float(np.nanmedian(dC_h[-1, ...]))
            print(f"[SANITY] {specie_conc} ΔC_h median: ΔC0={dC0_med:.6g} ΔCend(H={self.H_max})={dCend_med:.6g}")

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
        pulse_split: str = "equal",
        eps: float = 1e-18,
        store: bool = True,
    ) -> np.ndarray:
        dC_h = self.get_gas_deltaC_t(
            specie_conc=specie_conc,
            specie_emis=specie_emis,
            pulse_size=pulse_size,
            pulse_split=pulse_split,
            store=store,
        )

        dC0 = np.asarray(dC_h[0, ...])
        small = np.abs(dC0) < eps
        if np.any(small):
            n_bad = int(np.sum(small))
            print(f"[WARNING] {specie_conc} pIRF: {n_bad} member(s) have |ΔC0|<{eps:.1e}; using eps.")

        dC0_safe = np.where(small, eps, dC0)
        irf_h = dC_h / dC0_safe
        irf_h[0, ...] = 1.0

        if store:
            if specie_conc == "CO2":
                self.co2_f["CO2_IRF_fair_t"] = irf_h
            elif specie_conc == "CH4":
                self.ch4_f["CH4_IRF_fair_t"] = irf_h
            elif specie_conc == "N2O":
                self.n2o_f["N2O_IRF_fair_t"] = irf_h

        return irf_h
    
    
    # wrappers
    def get_co2_pIRF(
        self,
        pulse_size: float = 1.0,
        co2_emis: Tuple[str, ...] = ("CO2 FFI",),
        pulse_split: str = "equal",
    ):
        return self.get_gas_IRF_fair_t(
            specie_conc="CO2",
            specie_emis=co2_emis,
            pulse_size=pulse_size,
            pulse_split=pulse_split,
            store=True,
        )

    def get_ch4_pIRF(self, pulse_size: float = 1.0):
        return self.get_gas_IRF_fair_t("CH4", "CH4", pulse_size=pulse_size, pulse_split="equal", store=True)

    def get_n2o_pIRF(self, pulse_size: float = 1.0):
        return self.get_gas_IRF_fair_t("N2O", "N2O", pulse_size=pulse_size, pulse_split="equal", store=True)
    
    
    ############### above get_x_pIRF to run each gas seperately for testing, now we run all three gas together 
    def _run_baseline_once(self):
        """Run and cache a baseline FaIR run (deepcopy-safe)."""
        f_base = copy.deepcopy(self.f)
        f_base.run()
        conc_base = {
            "CO2": np.asarray(f_base.concentration.loc[dict(scenario=self.scn, specie="CO2")].values).copy(),
            "CH4": np.asarray(f_base.concentration.loc[dict(scenario=self.scn, specie="CH4")].values).copy(),
            "N2O": np.asarray(f_base.concentration.loc[dict(scenario=self.scn, specie="N2O")].values).copy(),
        }
        return f_base, conc_base


    def get_all_pIRF(
        self,
        pulse_size: float = 1.0,
        co2_emis: Tuple[str, ...] = ("CO2 FFI",),
        pulse_split: str = "equal",
        eps: float = 1e-18,
    ):
        """
        Compute pIRF for CO2, CH4, N2O using:
          - ONE baseline run
          - ONE perturbed run per gas (3 runs)
        Total: 4 runs instead of 6 (≈33% faster).
        Returns dict: {"CO2": irf, "CH4": irf, "N2O": irf}
        Also stores in self.co2_f/ch4_f/n2o_f like wrappers.
        """
        # indices (same logic as above working code)
        pulse_year = self._get_pulse_year()
        tb_pulse = self._get_timebound_index(pulse_year)
        tp_pulse = self._get_timepoint_index(pulse_year)

        dt = 1.0 / float(self.ts_per_year)
        tp_axis_self = _as_float_array(self.f.timepoints)
        tb_axis_self = _as_float_array(self.f.timebounds)
        tp_pulse_time = float(tp_axis_self[tp_pulse])
        
        #tb_resp0 = _first_index_ge(tb_axis_self, tp_pulse_time - 1e-12)
        tb_resp0 = _first_index_ge(tb_axis_self, tp_pulse_time + dt - 1e-12)

        nH = self.H_max * self.ts_per_year + 1
        tb1 = tb_resp0 + nH
        if tb1 > tb_axis_self.size:
            raise ValueError("Horizon exceeds timebounds.")

        # baseline once
        f_base, conc_base = self._run_baseline_once()

        if self.debug:
            tb_axis = _as_float_array(f_base.timebounds)
            tp_axis = _as_float_array(f_base.timepoints)
            print(
                f"[pIRF-base] scn={self.scn} year_index={self.year_index} pulse_year={pulse_year} "
                f"tb_pulse={tb_pulse}(tb={tb_axis[tb_pulse]:.1f}) tp_pulse={tp_pulse}(tp={tp_axis[tp_pulse]:.3f})"
            )

        def _pulse_and_run(specie_conc: str, emis_spec: Union[str, Tuple[str, ...]]):
            # perturbed run (independent deepcopy, safe)
            f_pert = copy.deepcopy(self.f)

            # apply pulse
            if isinstance(emis_spec, str):
                emis_list = (emis_spec,)
            else:
                emis_list = tuple(emis_spec)

            if pulse_split == "equal":
                parts = [float(pulse_size) / len(emis_list)] * len(emis_list)
            else:
                parts = [float(pulse_size)] * len(emis_list)

            if self.debug:
                print(f"[pIRF-debug] gas={specie_conc} emis={','.join(emis_list)} tp_pulse={tp_pulse}(tp={tp_axis_self[tp_pulse]:.3f})")

            for sp_e, p_add in zip(emis_list, parts):
                em_slice = f_pert.emissions.loc[dict(scenario=self.scn, specie=sp_e)].copy()
                em_vals = np.asarray(em_slice.values)

                if em_vals.ndim == 1:
                    before = float(em_vals[tp_pulse])
                    em_vals[tp_pulse] = before + p_add
                    after = float(em_vals[tp_pulse])
                else:
                    before = float(np.nanmedian(em_vals[tp_pulse, ...]))
                    em_vals[tp_pulse, ...] = em_vals[tp_pulse, ...] + p_add
                    after = float(np.nanmedian(em_vals[tp_pulse, ...]))

                em_slice.values = em_vals
                f_pert.emissions.loc[dict(scenario=self.scn, specie=sp_e)] = em_slice

                if self.debug:
                    print(f"[SANITY] {sp_e} emis@tp_pulse median before/after: {before:.6g} -> {after:.6g} (add {p_add:.6g})")

            f_pert.run()

            conc_pert = np.asarray(
                f_pert.concentration.loc[dict(scenario=self.scn, specie=specie_conc)].values
            ).copy()

            dC = conc_pert - conc_base[specie_conc]
            dC_h = dC[tb_resp0:tb1, ...].copy()

            # IRF
            dC0 = np.asarray(dC_h[0, ...])
            small = np.abs(dC0) < eps
            dC0_safe = np.where(small, eps, dC0)
            irf = dC_h / dC0_safe
            irf[0, ...] = 1.0

            return dC_h, irf

        # run 3 perturbed cases
        dC_co2, irf_co2 = _pulse_and_run("CO2", co2_emis)
        dC_ch4, irf_ch4 = _pulse_and_run("CH4", "CH4")
        dC_n2o, irf_n2o = _pulse_and_run("N2O", "N2O")

        # store like wrappers
        self.co2_f["CO2_dC_ppm_t"] = dC_co2
        self.co2_f["CO2_IRF_fair_t"] = irf_co2
        self.ch4_f["CH4_dC_ppb_t"] = dC_ch4
        self.ch4_f["CH4_IRF_fair_t"] = irf_ch4
        self.n2o_f["N2O_dC_ppb_t"] = dC_n2o
        self.n2o_f["N2O_IRF_fair_t"] = irf_n2o

        return {"CO2": irf_co2, "CH4": irf_ch4, "N2O": irf_n2o}


    
def majorghg_get_irf(
    f,
    scn: str,
    H_max: int = 100,
    ts_per_year: int = 1,
    year_index: int = 269,
    debug: bool = False,
):
    return majorghg_get_f(
        f=f,
        scn=scn,
        H_max=H_max,
        ts_per_year=ts_per_year,
        year_index=year_index,
        debug=debug,
    )