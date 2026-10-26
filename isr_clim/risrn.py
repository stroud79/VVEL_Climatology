"""RISR-N long-pulse vector velocities (Madrigal kindat 5960) -> quality-
controlled records, geographic east and north with their errors.

Archive layout (globalDownload.py --tree): ROOT/<YYYY>/ran/<exp>/ran*.hdf5.

Records are kept individually, as the FPI pipeline keeps its observations:
each is paired with the model step within 30 min of it and weighted by its
own error in the diurnal means (climatology.py). Three things about the
archive are handled here:

* Each experiment holds the same measurements several times, integrated
  over different times (30-500 s). One file per experiment is used, the
  first cadence available in config.PREFERRED_CADENCE_S; where an
  experiment has two files at that cadence, the one with more finite
  vectors is kept.
* Files come in two layouts, told apart by their columns and checked
  against the files' own direction angles. Most give vi1/vi2 (geographic
  east/north: atan2(vi1, vi2) equals gnangle, "CW angle from GEO North");
  some 2016-2018 files give vipe/vipn instead (perpendicular east/north in
  the local magnetic frame: atan2(vipe, vipn) equals nangle, "from GM
  North"). The latter are rotated to geographic with the IGRF declination
  at the bin, and their errors propagated (geometry.mag_to_geo_errors).
  vipn is perpendicular to B rather than horizontal; at Resolute
  (inclination ~87-88 deg) the two differ by ~0.1%, which is neglected.
* Quality control, on the file's own components: both horizontal errors
  < MAX_VELOCITY_ERROR and the speed < MAX_SPEED. Errors of exactly 1 m/s
  mark degenerate fits, not measurements (config.MIN_VELOCITY_ERROR), and
  are dropped: weighted by 1/error^2 they would decide whole hourly bins.
"""

from __future__ import annotations

import logging
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from . import config, geometry

log = logging.getLogger(__name__)

LAYOUTS = {  # frame: (east, north, east error, north error)
    "geo": ("vi1", "vi2", "dvi1", "dvi2"),
    "mag": ("vipe", "vipn", "dvipe", "dvipn"),
}
RECORD_COLS = ("time", "cgm_lat", "ve_geo", "vn_geo", "dve_geo", "dvn_geo", "frame", "experiment")


def _layout(names) -> str:
    for frame, cols in LAYOUTS.items():
        if set(cols) <= set(names):
            return frame
    raise ValueError(f"no known velocity columns in {sorted(names)}")


def scan(root: Path) -> pd.DataFrame:
    """One row per file: experiment, cadence, layout, record count, span."""
    rows = []
    for path in sorted(root.glob("*/ran/*/*.hdf5")):
        try:
            with h5py.File(path, "r") as h:
                tab = h["Data/Table Layout"]
                names = tab.dtype.names
                frame = _layout(names)
                e, n = LAYOUTS[frame][:2]
                a = tab.fields(["ut1_unix", "ut2_unix", "inttms", e, n])[:]
            rows.append(dict(
                path=str(path), experiment=path.parent.name, frame=frame,
                cadence=int(np.round(np.median(a["inttms"]))),
                n_records=len(a), n_finite=int(np.sum(np.isfinite(a[e]) & np.isfinite(a[n]))),
                t0=pd.to_datetime(a["ut1_unix"].min(), unit="s"),
                t1=pd.to_datetime(a["ut2_unix"].max(), unit="s"), error=""))
        except Exception as err:  # unreadable file: reported, not fatal
            rows.append(dict(path=str(path), experiment=path.parent.name, error=repr(err)))
    return pd.DataFrame(rows)


def select(files: pd.DataFrame) -> pd.DataFrame:
    """One file per experiment (see module docstring)."""
    ok = files[files["error"] == ""].copy()
    rank = {c: i for i, c in enumerate(config.PREFERRED_CADENCE_S)}
    ok["rank"] = ok["cadence"].map(rank).fillna(len(rank))
    ok = ok.sort_values(["experiment", "rank", "n_finite"], ascending=[True, True, False])
    return ok.groupby("experiment", as_index=False).first().drop(columns="rank")


def load_file(path: str) -> pd.DataFrame:
    """Quality-controlled records of one file, geographic components."""
    with h5py.File(path, "r") as h:
        a = h["Data/Table Layout"][:]
    frame = _layout(a.dtype.names)
    e, n, de, dn = LAYOUTS[frame]
    df = pd.DataFrame({
        "time": pd.to_datetime(0.5 * (a["ut1_unix"] + a["ut2_unix"]), unit="s"),
        "cgm_lat": np.round(a["cgm_lat"].astype(float), 4),
        "ve": a[e].astype(float), "vn": a[n].astype(float),
        "dve": a[de].astype(float), "dvn": a[dn].astype(float),
    })
    good = (np.isfinite(df.ve) & np.isfinite(df.vn) & np.isfinite(df.dve) & np.isfinite(df.dvn)
            & (df.dve < config.MAX_VELOCITY_ERROR) & (df.dvn < config.MAX_VELOCITY_ERROR)
            & (df.dve > config.MIN_VELOCITY_ERROR) & (df.dvn > config.MIN_VELOCITY_ERROR)
            & (np.hypot(df.ve, df.vn) < config.MAX_SPEED))
    df = df[good].reset_index(drop=True)
    if frame == "geo":
        df["ve_geo"], df["vn_geo"], df["dve_geo"], df["dvn_geo"] = df.ve, df.vn, df.dve, df.dvn
    else:
        decl = geometry.positions_for(df["time"], df["cgm_lat"])["decl"].to_numpy()
        df["ve_geo"], df["vn_geo"] = geometry.mag_to_geo(df.ve, df.vn, decl)
        df["dve_geo"], df["dvn_geo"] = geometry.mag_to_geo_errors(df.dve, df.dvn, decl)
    df["frame"] = frame
    return df.drop(columns=["ve", "vn", "dve", "dvn"])


def run(root: Path, out_dir: Path) -> pd.DataFrame:
    """Scan, select and quality-control the archive. Writes risrn_files.csv
    (every file and whether it was used) and risrn_records.csv.gz."""
    out_dir.mkdir(parents=True, exist_ok=True)
    files = scan(root)
    if files.empty:
        raise FileNotFoundError(f"no RISR-N files under {root}")
    chosen = select(files)
    files["used"] = files["path"].isin(chosen["path"])
    parts = []
    for _, f in chosen.iterrows():
        df = load_file(f["path"])
        files.loc[files["path"] == f["path"], "n_good"] = len(df)
        df["experiment"] = f["experiment"]
        parts.append(df)
    records = pd.concat(parts, ignore_index=True).sort_values(["time", "cgm_lat"])[list(RECORD_COLS)]
    files.to_csv(out_dir / "risrn_files.csv", index=False)
    records.to_csv(out_dir / "risrn_records.csv.gz", index=False, float_format="%.6g")
    bad = files[files["error"] != ""]
    log.info("RISR-N: %d files, %d experiments used, %d unreadable, %d records kept",
             len(files), len(chosen), len(bad), len(records))
    return records


def load_records(out_dir: Path) -> pd.DataFrame:
    path = out_dir / "risrn_records.csv.gz"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; run the 'risrn' step first")
    return pd.read_csv(path, parse_dates=["time"])
