"""The two figures that test the FPI findings at the station.

* ``diurnal``: the FPI seasonal diurnal figure (Plotting/climatology_plots.py
  plot_seasonal_diurnal_climatology) drawn for the RISR-N ion drift instead
  of the FPI neutral wind: observed weighted hourly means with their error
  and spread, the WACCM-X SD line through its own 3-hourly bins (bridged
  for display), the weighted mean residual, and the residual histogram.
  Same layout, local-time axis (noon to noon, station's standard time),
  colours and style, so it can be laid beside the FPI figure.
* ``residual_comparison``: per season, the WACCM-X SD residual of the ion
  drift (vs RISR-N) and of the neutral wind (vs FPI) on one axis. If the
  neutral-wind errors come from the model's high-latitude convection, the
  two follow each other in local time and sign.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # batch jobs have no display
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from . import config, fpi  # noqa: E402  (fpi applies the FPI plots' style)
from .climatology import Diurnal  # noqa: E402

PLOT_DATES = pd.date_range("2000-01-01 12:00", periods=25, freq="h")
SEASON_NAMES = {"DJF": "Winter", "MAM": "Spring", "JJA": "Summer", "SON": "Autumn"}
GAP_BINS = fpi.display_gap_bins(fpi.SD_STEP_MINUTES)
MODEL = fpi.MODEL_DISPLAY[config.MODEL_KEY]
ION_COLOR, NEUTRAL_COLOR = "teal", fpi.SINGLE_MODEL_RESIDUAL_COLOR


def _noon(values) -> np.ndarray:
    """24 local-hour bins -> 25 values at 12:00 ... 12:00, as the FPI panels."""
    v = np.roll(np.asarray(values, float), -12)
    return np.append(v, v[0])


def _markers(values) -> list:
    return list(np.flatnonzero(np.isfinite(values)))


def _time_axis(ax, xlabel=None):
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:00"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=4))
    ax.set_xlim(PLOT_DATES[0], PLOT_DATES[-1])
    if xlabel:
        ax.set_xlabel(xlabel)


def _season(d: Diurnal, season: str) -> pd.DataFrame:
    return d.table[d.table["season"] == season].sort_values("hour")


def _span(d: Diurnal) -> str:
    return f"{d.span[0]:%Y-%m-%d} to {d.span[1]:%Y-%m-%d}"


def diurnal(d: Diurnal, component: str, station: fpi.Station, band, out: Path) -> Path:
    seasons = [s for s in config.SEASON_ORDER if (_season(d, s)["obs_n"] > 0).any()]
    res = (d.pairs["model"] - d.pairs["obs"]).to_numpy()
    res_bound = max(np.percentile(np.abs(res), 99) * 1.2, 150) if len(res) else 200
    hist_bins = np.linspace(-res_bound, res_bound, 21)
    vals = np.r_[d.obs["obs"].to_numpy(), d.steps["model"].to_numpy()]
    vals = vals[np.isfinite(vals)]
    v_min, v_max = (np.percentile(vals, 1), np.percentile(vals, 99)) if len(vals) else (-300, 300)
    pad = max((v_max - v_min) * 0.15, 50)

    fig, axes = plt.subplots(len(seasons), 3, figsize=(18, 3 * len(seasons)),
                             gridspec_kw={"width_ratios": [2.5, 1.5, 1]}, squeeze=False)
    for i, season in enumerate(seasons):
        t = _season(d, season)
        ax, ax_res, ax_hist = axes[i]
        last = i == len(seasons) - 1
        xlabel = fpi.local_time_label(station.tz) if last else None

        ax.axhline(0, **fpi.ZERO_LINE_KW)
        o, eo, so = _noon(t["obs_mean"]), _noon(t["obs_err"]), _noon(t["obs_std"])
        ax.errorbar(PLOT_DATES, o, yerr=eo, fmt="o", color="black", ecolor="black", elinewidth=1, capsize=3,
                    label="RISR-N")
        ax.fill_between(PLOT_DATES, o - so, o + so, color="black", alpha=0.2)
        m_plot = _noon(t["model_mean"])
        m_line = _noon(fpi.fill_hour_gaps(t["model_mean"], GAP_BINS))
        em_line = _noon(fpi.fill_hour_gaps(t["model_std"], GAP_BINS))
        ax.plot(PLOT_DATES, m_line, "-", linewidth=2, color=MODEL["color"], label=MODEL["label"],
                marker="o", markevery=_markers(m_plot))
        ax.fill_between(PLOT_DATES, m_line - em_line, m_line + em_line, color=MODEL["color"], alpha=0.2)
        ax.set_ylim(v_min - pad, v_max + pad)
        ax.set_title(f"Season: {SEASON_NAMES[season]}")
        ax.grid(True, linestyle="--")
        ax.set_ylabel(f"{component} Ion Drift (m/s)")
        _time_axis(ax, xlabel)
        if i == 0:
            ax.legend(loc="upper right")

        ax_res.axhline(0, **fpi.ZERO_LINE_KW)
        r_plot = _noon(t["res_mean"])
        r_line = _noon(fpi.fill_hour_gaps(t["res_mean"], GAP_BINS))
        r_std = _noon(fpi.fill_hour_gaps(t["res_std"], GAP_BINS))
        ax_res.plot(PLOT_DATES, r_line, "-", color=fpi.SINGLE_MODEL_RESIDUAL_COLOR, linewidth=1.5,
                    marker="o", markevery=_markers(r_plot), label="Weighted Mean Res")
        ax_res.fill_between(PLOT_DATES, r_line - r_std, r_line + r_std, color=fpi.SINGLE_MODEL_RESIDUAL_COLOR,
                            alpha=0.25, label="+/-1 s.d.")
        ax_res.set_ylim(-res_bound, res_bound)
        ax_res.grid(True, linestyle="--")
        ax_res.set_ylabel("Δ (m/s)")
        ax_res.legend(loc="upper right", fontsize="small")
        _time_axis(ax_res)

        r_s = (d.pairs.loc[d.pairs["season"] == season, "model"]
               - d.pairs.loc[d.pairs["season"] == season, "obs"]).to_numpy()
        ax_hist.hist(r_s, bins=hist_bins, color=fpi.SINGLE_MODEL_RESIDUAL_COLOR, edgecolor="black", alpha=0.7)
        ax_hist.text(0.05, 0.95, f"N = {len(r_s)}", transform=ax_hist.transAxes, fontsize=10,
                     verticalalignment="top", bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
        if len(r_s):
            ax_hist.axvline(r_s.mean(), color="red", linestyle="dashed", linewidth=2, label=f"Mean: {r_s.mean():.1f}")
            ax_hist.legend(fontsize="small")
        ax_hist.axvline(0, color="gray", linewidth=1.5)
        ax_hist.set_xlim(-res_bound, res_bound)
        ax_hist.grid(True, linestyle="--", axis="y")
        ax_hist.set_ylabel("Frequency")
        if last:
            ax_hist.set_xlabel(f"{MODEL['label']} - RISR-N (m/s)")

    fig.suptitle(f"Climatological Diurnal Variation of {component} Ion Drift\n{station.plot_name} | "
                 f"RISR-N CGM latitude {band[0]}-{band[1]}° | {_span(d)}")
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def residual_comparison(ion: Diurnal, neutral: Diurnal, ion_stats: pd.DataFrame, neutral_stats: pd.DataFrame,
                        component: str, station: fpi.Station, band, out: Path) -> Path:
    """Weighted mean residual (model - obs) +/- 1 s.d. per local hour, ion
    drift and neutral wind, one row per season either instrument has."""
    seasons = [s for s in config.SEASON_ORDER
               if (_season(ion, s)["res_n"] > 0).any() or (_season(neutral, s)["res_n"] > 0).any()]
    # Scaled to the plotted curves (mean +/- 1 s.d.), not the raw residuals,
    # whose ion-drift tails would flatten the neutral-wind curve.
    bound = 150.0
    for d in (ion, neutral):
        t = d.table[d.table["season"].isin(seasons)]
        edge = np.abs(np.r_[t["res_mean"] + t["res_std"], t["res_mean"] - t["res_std"]])
        if np.isfinite(edge).any():
            bound = max(bound, 1.1 * np.nanmax(edge))

    def label(name, stats, season):
        s = stats[stats["Season"] == season]
        if s.empty:
            return name
        s = s.iloc[0]
        return f"{name}: bias {s['Bias (m/s)']:.0f}, r {s['r']:.2f}, σ$_m$/σ$_o$ {s['Std ratio']:.2f} (N={s['N']})"

    fig, axes = plt.subplots(len(seasons), 1, figsize=(11, 3.2 * len(seasons)), squeeze=False, sharex=True)
    for ax, season in zip(axes[:, 0], seasons):
        ax.axhline(0, **fpi.ZERO_LINE_KW)
        for d, stats, name, color in (
                (ion, ion_stats, f"Ion drift, {MODEL['label']} - RISR-N", ION_COLOR),
                (neutral, neutral_stats, f"Neutral wind, {MODEL['label']} - FPI", NEUTRAL_COLOR)):
            t = _season(d, season)
            if t.empty or not (t["res_n"] > 0).any():
                continue
            r_plot = _noon(t["res_mean"])
            r_line = _noon(fpi.fill_hour_gaps(t["res_mean"], GAP_BINS))
            r_std = _noon(fpi.fill_hour_gaps(t["res_std"], GAP_BINS))
            ax.plot(PLOT_DATES, r_line, "-", color=color, linewidth=1.8, marker="o",
                    markevery=_markers(r_plot), label=label(name, stats, season))
            ax.fill_between(PLOT_DATES, r_line - r_std, r_line + r_std, color=color, alpha=0.15)
        ax.set_ylim(-bound, bound)
        ax.set_ylabel("Model - Obs (m/s)")
        ax.set_title(f"Season: {SEASON_NAMES[season]}")
        ax.grid(True, linestyle="--")
        ax.legend(loc="upper right", fontsize="small")
        _time_axis(ax)
    _time_axis(axes[-1, 0], fpi.local_time_label(station.tz))
    fig.suptitle(f"{MODEL['label']} {component.lower()} residuals at {station.plot_name}: ion drift vs neutral wind\n"
                 f"RISR-N CGM latitude {band[0]}-{band[1]}° ({_span(ion)}); "
                 f"FPI 45° looks ({_span(neutral)}); weighted mean ± 1 s.d.")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out
