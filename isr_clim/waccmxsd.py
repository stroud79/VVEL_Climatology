"""WACCM-X SD ion drifts at the RISR-N bin positions.

Each daily file (f.e20.FXSD.f19_f19.001.cam.h1.<date>-00000.nc, FPI's
LAYOUTS["sd"]) holds 8 instantaneous outputs, 00, 03, ..., 21 UT. The
3-hourly steps are kept as they are: nothing here or later interpolates
between them (see climatology.pair). Per step, following the FM_V2 WACCM-X
provider and integrator:

1. UI (zonal) and VI (meridional) ion drifts, geographic, from edynamo,
   and Z3 are refused unless their units are m/s and m (require_units).
2. Geometric height from Z3 with WGS84 gravity at each column's latitude
   (geometric_from_geopotential, EQUATIONS.md 1.2-1.3).
3. Linear interpolation in latitude and longitude to each bin's
   geographic position, the longitude seam closed (prepare_lon_interp).
4. Linear interpolation in geometric height onto a uniform MODEL_DZ_KM grid
   through MODEL_ALT_RANGE_KM, the layer the RISR-N vectors are fitted
   over, NaN above the model top; the mean over that grid is the
   height-weighted layer mean. The model top is at ~450-600 km, so the
   layer stops there (z_top_km). The SD drifts vary by a few per cent at
   most over that layer (E x B is nearly height-independent along the
   near-vertical field lines), so this choice is not critical: it changes
   the drifts by under 1 m/s from a plain mean over model levels.

The result is cached as one small CSV per day (CACHE/<date>.csv), so the
archive is read once. Days without UI/VI are skipped and listed in
CACHE/_skipped.csv; once the fields have been added, rerunning extracts
just those days. A step whose fields are all NaN (the last step of the
truncated 2008-01-11 and 2017-12-31 files) stays NaN and pairs with nothing.
"""

from __future__ import annotations

import logging
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from . import config, fpi, geometry

log = logging.getLogger(__name__)

FILE_GLOB = fpi.SD_FILE_PATTERN.format(date="*")
FILE_RE = re.compile(re.escape(fpi.SD_FILE_PATTERN).replace(re.escape("{date}"), r"(\d{4}-\d{2}-\d{2})") + "$")
UNITS = {"UI": ("m/s",), "VI": ("m/s",), "Z3": ("m",)}
CACHE_COLS = ["time", "cgm_lat", "ve_geo", "vn_geo", "z_top_km"]


def list_files(root: Path, start: str | None = None, end: str | None = None) -> list[tuple[str, Path]]:
    out = []
    for p in sorted(root.glob(FILE_GLOB)):
        m = FILE_RE.search(p.name)
        if m and (start is None or m.group(1) >= start) and (end is None or m.group(1) <= end):
            out.append((m.group(1), p))
    return out


def layer_mean(z_km: np.ndarray, values: np.ndarray, alt_range=config.MODEL_ALT_RANGE_KM,
               dz=config.MODEL_DZ_KM) -> np.ndarray:
    """Height-weighted mean of ``values`` over ``alt_range`` for columns
    (..., lev) with geometric heights ``z_km``: linear interpolation in
    height onto a uniform grid (NaN outside each column), then the mean of
    the finite grid values."""
    grid = np.arange(alt_range[0], alt_range[1] + 0.5 * dz, dz)
    zs = z_km.reshape(-1, z_km.shape[-1])
    vs = values.reshape(-1, values.shape[-1])
    out = np.full(len(zs), np.nan)
    for k, (z, v) in enumerate(zip(zs, vs)):
        ok = np.isfinite(z) & np.isfinite(v)
        if ok.sum() < 2:
            continue
        order = np.argsort(z[ok])
        g = np.interp(grid, z[ok][order], v[ok][order], left=np.nan, right=np.nan)
        if np.isfinite(g).any():
            out[k] = np.nanmean(g)
    return out.reshape(values.shape[:-1])


def extract_day(path: Path, cgm_lats: tuple, alt_range=config.MODEL_ALT_RANGE_KM,
                dz=config.MODEL_DZ_KM) -> pd.DataFrame:
    """Layer-mean UI/VI at each bin for every output time in ``path``."""
    with xr.open_dataset(path) as raw:
        missing = [v for v in UNITS if v not in raw.data_vars]
        if missing:
            raise KeyError(f"{path.name}: no {missing}")
        for v, units in UNITS.items():
            fpi.require_units(raw[v], units, path.name)
        times = pd.DatetimeIndex(raw["time"].values)
        pos = geometry.bin_positions(int(times[0].year), tuple(cgm_lats))
        plat = pos["glat"].to_numpy()
        # Only the rows around the bins: three grid rows of margin either side.
        dlat = float(np.abs(np.diff(raw["lat"].values[:2]))[0])
        ds = raw[list(UNITS)].sel(lat=slice(plat.min() - 3 * dlat, plat.max() + 3 * dlat))
        ds = ds.transpose("time", "lev", "lat", "lon").astype("float64")
        ds = ds.assign_coords(altitude=fpi.geometric_from_geopotential(ds["Z3"], ds["lat"], fpi.g0))
        ds, plon = fpi.prepare_lon_interp(ds, pos["glon"].to_numpy())
        h = ds.interp(lat=xr.DataArray(plat, dims="bin"), lon=xr.DataArray(plon, dims="bin"),
                      method="linear", kwargs={"fill_value": np.nan}).load()
    order = ("time", "bin", "lev")
    z = h["altitude"].transpose(*order).values / 1000.0
    ve = layer_mean(z, h["UI"].transpose(*order).values, alt_range, dz)
    vn = layer_mean(z, h["VI"].transpose(*order).values, alt_range, dz)
    nt, nb = ve.shape
    z_top = np.where(np.isfinite(ve) & np.isfinite(vn), np.nanmax(z, axis=-1), np.nan)
    return pd.DataFrame({
        "time": np.repeat(times.values, nb),
        "cgm_lat": np.tile(pos["cgm_lat"].to_numpy(), nt),
        "ve_geo": ve.ravel(), "vn_geo": vn.ravel(), "z_top_km": z_top.ravel(),
    })[CACHE_COLS]


def _cache_one(date: str, path: Path, cache: Path, cgm_lats: tuple) -> str:
    df = extract_day(path, cgm_lats)
    tmp = cache / f".{date}.csv.tmp"
    df.to_csv(tmp, index=False, float_format="%.6g")
    tmp.replace(cache / f"{date}.csv")
    return "extracted"


def run(root: Path, cache: Path, cgm_lats: tuple, workers: int = 1, start=None, end=None,
        overwrite: bool = False) -> int:
    """Extract every day not yet cached. Returns the number of failures."""
    cache.mkdir(parents=True, exist_ok=True)
    files = list_files(root, start, end)
    todo = [(d, p) for d, p in files if overwrite or not (cache / f"{d}.csv").exists()]
    log.info("WACCM-X SD: %d files in %s, %d to extract", len(files), root, len(todo))
    skipped = []
    with ProcessPoolExecutor(max_workers=max(1, workers)) as pool:
        jobs = {pool.submit(_cache_one, d, p, cache, tuple(cgm_lats)): (d, p) for d, p in todo}
        for n, job in enumerate(as_completed(jobs), 1):
            d, p = jobs[job]
            try:
                job.result()
            except Exception as err:   # one bad day must not stop the rest
                skipped.append(dict(date=d, path=str(p), reason=repr(err)))
            if n % 250 == 0 or n == len(todo):
                log.info("  [%d/%d] %d skipped so far", n, len(todo), len(skipped))
    pd.DataFrame(skipped, columns=["date", "path", "reason"]).to_csv(cache / "_skipped.csv", index=False)
    return len(skipped)


def load_cache(cache: Path, cgm_lats=None) -> pd.DataFrame:
    """Every cached day, optionally only the bins ``cgm_lats``."""
    parts = []
    for p in sorted(cache.glob("????-??-??.csv")):
        d = pd.read_csv(p, parse_dates=["time"])
        if list(d.columns) != CACHE_COLS:
            raise ValueError(f"{p} has columns {list(d.columns)}, not {CACHE_COLS}: it was made by an "
                             "earlier version; re-run the extract step with --overwrite")
        parts.append(d if cgm_lats is None else d[np.isin(d["cgm_lat"], cgm_lats)])
    if not parts:
        raise FileNotFoundError(f"no extracted WACCM-X SD days in {cache}")
    return pd.concat(parts, ignore_index=True).sort_values(["cgm_lat", "time"]).reset_index(drop=True)
