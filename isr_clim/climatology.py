"""Pairing, diurnal climatology and scores, made as the FPI plots make them
(FPI_ForwardModel/Plotting/climatology_plots.py,
plot_seasonal_diurnal_climatology), so that the RISR-N ion-drift numbers
can be read beside the FPI neutral-wind ones. The same functions are
applied to both instruments.

Time. Everything is selected and paired in UT. Each observation takes the
WACCM-X SD step nearest to it if that step lies within 30 min
(fpi.model_at_times, the FM_V2 matching window) and is otherwise
unpaired. The model is never interpolated between its 3-hourly steps, so
only observations within 30 min of 00, 03, ..., 21 UT are paired.

Local time is used only to bin: the station's time zone (from the FPI
registry coordinates) at its standard offset, without daylight saving,
bin h holding h-30 to h+30 min (fpi.local_hour_bin). With the offset
fixed, the 3-hourly UT steps fill the same every third local hour in
every season (civil time would split spring and autumn across two
interleaved sets); statistics come from those bins only, and the gaps
are bridged for display only (plots.py, fpi.fill_hour_gaps).

Model line. The FPI runs are obs_mapped: a model step exists only where an
observation lies within 30 min of it, so the FPI model line is the mean of
those steps, each counted once. Here likewise: the steps (per bin) that at
least one RISR-N record was paired with (obs_mapped_steps).

Means. Observations: inverse-variance weighted, w = 1 / max(err, 0.5)^2,
with the error of the weighted mean sqrt(1 / sum w) and the unweighted
spread. Residuals (model - obs): weighted the same way. Scores: on every
paired observation (not on hourly means), as in the FPI wind_stats CSV.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import xarray as xr

from . import config, fpi

MIN_PAIRS = 5        # FPI: statistics only when more than 4 pairs
ERR_FLOOR = 0.5      # m/s, FPI: errors below this are raised to it


def in_band(df: pd.DataFrame, band=config.STATION_BAND) -> pd.DataFrame:
    return df[(df["cgm_lat"] >= band[0]) & (df["cgm_lat"] <= band[1])]


def pair(records: pd.DataFrame, model: pd.DataFrame, window=fpi.MATCH_WINDOW_MINUTES) -> pd.DataFrame:
    """``records`` with, per bin, the model step paired with each record:
    ``model_time`` and model_<value> for every value column of ``model``
    (NaN / NaT where no step lies within ``window`` minutes)."""
    values = [c for c in ("vn_geo", "ve_geo") if c in model]
    out = records.copy()
    out["model_time"] = pd.NaT
    for c in values:
        out[f"model_{c}"] = np.nan
    for lat, idx in out.groupby("cgm_lat").groups.items():
        m = model[np.isclose(model["cgm_lat"], lat)].drop_duplicates("time").sort_values("time")
        if m.empty:
            continue
        # Pair by step index with FPI's own rule, so the step itself is known.
        step = xr.DataArray(np.arange(len(m), dtype=float), coords={"time": m["time"].to_numpy()}, dims="time")
        k = fpi.model_at_times(step, out.loc[idx, "time"].to_numpy(), window)
        ok = np.isfinite(k)
        rows, k = np.asarray(idx)[ok], k[ok].astype(int)
        out.loc[rows, "model_time"] = m["time"].to_numpy()[k]
        for c in values:
            out.loc[rows, f"model_{c}"] = m[c].to_numpy()[k]
    return out


def obs_mapped_steps(paired: pd.DataFrame, value: str) -> tuple[np.ndarray, np.ndarray]:
    """(time, value) of every model step, per bin, paired with at least one
    record: what an obs_mapped forward-model run would have simulated."""
    s = paired.dropna(subset=["model_time"]).drop_duplicates(["cgm_lat", "model_time"])
    return s["model_time"].to_numpy("datetime64[ns]"), s[f"model_{value}"].to_numpy(float)


@dataclass
class Diurnal:
    """Season x local-hour climatology of one instrument and component.

    ``table``: one row per (season, hour 0-23) with obs_mean, obs_err,
    obs_std, obs_n, model_mean, model_std, model_n, res_mean, res_std,
    res_n. ``obs``: every observation (time, season, hour, obs, model,
    weight w, residual res; model and res NaN where unpaired). ``steps``:
    the model steps the line is drawn from. ``span``: (first, last) time of
    the observations.
    """
    table: pd.DataFrame
    obs: pd.DataFrame
    steps: pd.DataFrame
    span: tuple

    @property
    def pairs(self) -> pd.DataFrame:
        return self.obs[np.isfinite(self.obs["res"])]


def _weighted(g, value, weight):
    sw = g[weight].sum()
    mean = (g[value] * g[weight]).sum() / sw if sw > 0 else np.nan
    var = (g[value] ** 2 * g[weight]).sum() / sw - mean ** 2 if sw > 0 else np.nan
    return pd.Series({"mean": mean, "err": np.sqrt(1.0 / sw) if sw > 0 else np.nan,
                      "std": np.sqrt(max(var, 0.0)) if np.isfinite(var) else np.nan,
                      "n": len(g)})


def _grouped(df, value, keys):
    """_weighted per group, with its columns even when ``df`` is empty."""
    if df.empty:
        return pd.DataFrame(columns=["mean", "err", "std", "n"],
                            index=pd.MultiIndex.from_arrays([[], []], names=keys), dtype=float)
    return df.groupby(keys)[[value, "w"]].apply(_weighted, value, "w")


def seasonal_diurnal(time, obs, err, model_at_obs, model_time, model_value, tz: str) -> Diurnal:
    """Season x local-hour means of the observations, of the model steps
    (``model_time``, ``model_value``) and of the paired residuals."""
    o = pd.DataFrame({"time": np.asarray(time, "datetime64[ns]"), "obs": np.asarray(obs, float),
                      "err": np.asarray(err, float), "model": np.asarray(model_at_obs, float)})
    o = o[np.isfinite(o["obs"])]
    # FPI: safe_dv = dv.where(dv > 0.5, 0.5), which also gives a missing error 0.5.
    o["w"] = 1.0 / np.where(o["err"] > ERR_FLOOR, o["err"], ERR_FLOOR) ** 2
    o["season"] = pd.DatetimeIndex(o["time"]).month.map(config.SEASONS)
    o["hour"] = np.asarray(fpi.local_hour_bin(o["time"].to_numpy(), tz))
    o["res"] = o["model"] - o["obs"]

    m = pd.DataFrame({"time": np.asarray(model_time, "datetime64[ns]"), "model": np.asarray(model_value, float)})
    m = m[np.isfinite(m["model"])]
    m["season"] = pd.DatetimeIndex(m["time"]).month.map(config.SEASONS)
    m["hour"] = np.asarray(fpi.local_hour_bin(m["time"].to_numpy(), tz))

    p = o[np.isfinite(o["res"])]
    keys = ["season", "hour"]
    obs_t = _grouped(o, "obs", keys).add_prefix("obs_")
    obs_t["obs_std"] = o.groupby(keys)["obs"].std(ddof=0)
    res_t = _grouped(p, "res", keys).add_prefix("res_").drop(columns="res_err")
    mod_t = m.groupby(keys)["model"].agg(model_mean="mean", model_std=lambda s: s.std(ddof=0),
                                         model_n="count")
    seasons = [s for s in config.SEASON_ORDER if s in set(o["season"])]
    full = pd.MultiIndex.from_product([seasons, range(24)], names=keys)
    table = pd.concat([obs_t, mod_t, res_t], axis=1).reindex(full).reset_index()
    for c in ("obs_n", "model_n", "res_n"):
        table[c] = table[c].fillna(0).astype(int)
    span = (o["time"].min(), o["time"].max()) if len(o) else (pd.NaT, pd.NaT)
    return Diurnal(table, o.drop(columns="err").reset_index(drop=True), m.reset_index(drop=True), span)


def scores(d: Diurnal, labels: dict, hours: dict | None = None, hours_label: str = "all") -> pd.DataFrame:
    """FPI wind_stats columns per season (model - obs), plus the standard
    deviations behind the FPI Taylor diagrams (ddof=1) and their ratio.
    ``hours``: {season: set of local hours} to restrict to, e.g. the hours
    the other instrument observes, named ``hours_label`` in the Hours
    column; ``labels`` are leading columns."""
    rows = []
    for season in config.SEASON_ORDER:
        s = d.pairs[d.pairs["season"] == season]
        if hours is not None:
            s = s[s["hour"].isin(hours.get(season, ()))]
        if len(s) < MIN_PAIRS:
            continue
        o, m = s["obs"].to_numpy(), s["model"].to_numpy()
        res = m - o
        so, sm = np.std(o, ddof=1), np.std(m, ddof=1)
        rows.append({**labels, "Season": season, "Hours": hours_label,
                     "Start_Date": d.span[0].strftime("%Y-%m-%d"), "End_Date": d.span[1].strftime("%Y-%m-%d"),
                     "N": len(s), "Bias (m/s)": round(float(res.mean()), 1),
                     "RMSE (m/s)": round(float(np.sqrt(np.mean(res ** 2))), 1),
                     "MAE (m/s)": round(float(np.abs(res).mean()), 1),
                     "r": round(float(np.corrcoef(o, m)[0, 1]), 3) if so > 0 and sm > 0 else np.nan,
                     "Obs std (m/s)": round(float(so), 1), "Model std (m/s)": round(float(sm), 1),
                     "Std ratio": round(float(sm / so), 3) if so > 0 else np.nan})
    return pd.DataFrame(rows)


def paired_hours(d: Diurnal) -> dict:
    """{season: local hours holding at least one paired observation}."""
    return {s: set(g["hour"]) for s, g in d.pairs.groupby("season")}
