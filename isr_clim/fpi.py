"""The parts of FPI_ForwardModel this package uses, imported from it rather
than copied, so the ion-drift comparison is made exactly as the FPI one:

* the site registry (FPI_ForwardModel/config/sites.yaml): the station's
  coordinates and display name, and the WACCM-X SD archive path;
* FM_V2: geometric height from Z3 (WGS84 gravity), the units check every
  provider makes, longitude-seam handling for horizontal interpolation, the
  SD archive's file naming and the model-observation matching window;
* Plotting/climatology_plots.py: the pairing rule (nearest model step
  within 30 min, never interpolated in time), the local-standard-time
  binning and its axis label, the
  3-hourly display bridging, and the model labels and colours. Importing it
  also applies the FPI plots' style.

It also loads the FPI observations and the forward-model run for the
station, so their residuals can be set beside the ion-drift residuals.
"""

from __future__ import annotations

import glob
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import xarray as xr

from . import config

if not (config.FPI_ROOT / "fpi_common").is_dir():
    raise ImportError(f"FPI_ForwardModel not found at {config.FPI_ROOT} (set ISR_FPI_ROOT)")
for _p in (config.FPI_ROOT, config.FPI_ROOT / "FM_V2", config.FPI_ROOT / "Plotting"):
    if str(_p) not in sys.path:
        sys.path.append(str(_p))

from climatology_plots import (  # noqa: E402
    DEFAULT_MATCH_WINDOW_MINUTES,
    MODEL_DISPLAY,
    SINGLE_MODEL_RESIDUAL_COLOR,
    ZERO_LINE_KW,
    display_gap_bins,
    extract_model_component,
    fill_hour_gaps,
    local_hour_bin,
    local_time_label,
    model_at_times,
    model_cadence,
)
from fpi_common.registry import load_registry  # noqa: E402
from src.constants import g0  # noqa: E402
from src.core.standard import TIME_TOLERANCE, require_units  # noqa: E402
from src.providers.waccmx import LAYOUTS  # noqa: E402
from src.utils.height import geometric_from_geopotential  # noqa: E402
from src.utils.longitude import prepare_lon_interp  # noqa: E402

# One window everywhere: the forward model's own (FM_V2), which the FPI
# plots use for pairing.
MATCH_WINDOW_MINUTES = int(TIME_TOLERANCE / np.timedelta64(1, "m"))
assert MATCH_WINDOW_MINUTES == DEFAULT_MATCH_WINDOW_MINUTES

SD_FILE_PATTERN = LAYOUTS["sd"]          # f.e20.FXSD.f19_f19.001.cam.h1.{date}-00000.nc
SD_STEP_MINUTES = 180                    # the SD files are 3-hourly, 00-21 UT

__all__ = [
    "MATCH_WINDOW_MINUTES", "MODEL_DISPLAY", "SD_FILE_PATTERN", "SD_STEP_MINUTES",
    "SINGLE_MODEL_RESIDUAL_COLOR", "Station", "ZERO_LINE_KW", "display_gap_bins", "fill_hour_gaps",
    "fpi_meridional_zonal", "g0", "geometric_from_geopotential", "local_hour_bin", "local_time_label",
    "model_at_times",
    "prepare_lon_interp", "require_units", "sd_dir", "station",
]


@dataclass(frozen=True)
class Station:
    name: str
    plot_name: str
    lat: float
    lon: float
    tz: str


@lru_cache(maxsize=None)
def _registry():
    return load_registry()


@lru_cache(maxsize=None)
def station(name: str = config.SITE) -> Station:
    """The FPI station from the registry, with the time zone the FPI plots
    bin by (Plotting/climatology_plots.py main: TimezoneFinder at the site)."""
    from timezonefinder import TimezoneFinder

    s = _registry().site(name)
    tz = TimezoneFinder().timezone_at(lng=s.lon_deg, lat=s.lat_deg)
    return Station(s.name, s.plot_name, s.lat_deg, s.lon_deg, tz)


def sd_dir() -> Path:
    """The WACCM-X SD archive (registry ``waccmx_sd_dir``, env FPI_WACCMX_SD_DIR)."""
    return _registry().model_dir("waccmx_sd_dir")


def load_fpi_outputs(site: str = config.SITE, combo: str = config.FPI_COMBO):
    """(filtered observations, forward-model run, look elevation) for
    ``site``, opened as climatology_plots.main opens them."""
    reg = _registry()
    obs_file = reg.filtered_obs_file(site)
    files = sorted(glob.glob(str(reg.simulation_dir(combo, site) / "data" / "out_*.nc")))
    if not obs_file.is_file() or not files:
        raise FileNotFoundError(f"no FPI observations ({obs_file}) or {combo} output for {site}")
    ds_obs = xr.open_dataset(obs_file).drop_vars("altitude", errors="ignore")
    sim = xr.open_mfdataset(files, combine="nested", concat_dim="time", parallel=False, join="outer",
                            coords="minimal", data_vars="minimal", compat="override")
    return ds_obs, sim, float(reg.site(site).elevation_deg)


def fpi_meridional_zonal(site: str = config.SITE, combo: str = config.FPI_COMBO, outputs=None) -> dict:
    """The FPI observations and the forward-model run for ``site``, per
    component, selected and signed exactly as climatology_plots does.

    Returns {component: dict(time, obs, err, model_at_obs, model_time,
    model_value)}: every observation at the site's elevation along the
    component's axis (azimuth within 2 deg, the opposite look negated so
    both read positive north / east), its uncertainty, the model paired
    with it (model_at_times), and the model series the FPI line is drawn
    from. ``outputs``: load_fpi_outputs() already loaded.
    """
    ds_obs, sim, elevation = outputs or load_fpi_outputs(site, combo)
    # As plot_seasonal_diurnal_climatology: observations within the run's span.
    ds_obs = ds_obs.where((ds_obs.timestamps >= sim.time.min()) & (ds_obs.timestamps <= sim.time.max()),
                          drop=True)
    _, window = model_cadence(sim)
    out = {}
    for component, (azm_pos, azm_neg) in (("Meridional", (0, 180)), ("Zonal", (90, 270))):
        pos = np.abs((ds_obs.azm - azm_pos + 180) % 360 - 180) <= 2.0
        neg = np.abs((ds_obs.azm - azm_neg + 180) % 360 - 180) <= 2.0
        sel = ds_obs.where((pos | neg) & (ds_obs.elm >= elevation - 2.0) & (ds_obs.elm <= elevation + 2.0),
                           drop=True)
        is_neg = np.abs((sel.azm - azm_neg + 180) % 360 - 180) <= 2.0
        v = xr.where(is_neg, -sel["v_horiz"], sel["v_horiz"]).values
        series = extract_model_component(sim, "horiz_wind", elevation, azm_pos, azm_neg)
        times = sel.timestamps.values
        out[component] = dict(
            time=times, obs=v, err=sel["dv_horiz"].values,
            model_at_obs=model_at_times(series, times, window),
            model_time=series.time.values, model_value=series.values,
        )
    return out
