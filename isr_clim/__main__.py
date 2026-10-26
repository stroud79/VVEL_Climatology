"""Command line for the RISR-N / WACCM-X SD ion-drift check of the FPI findings.

    python -m isr_clim risrn                 # RISR-N quality-controlled records
    python -m isr_clim extract --workers 8   # WACCM-X SD drifts at the bins (cached per day)
    python -m isr_clim climatology           # pairing, statistics and the validation figures
    python -m isr_clim all --workers 8       # all three in order

The station, its time zone and the WACCM-X SD archive come from
FPI_ForwardModel's registry (config/sites.yaml; isr_clim/fpi.py); the RISR-N
archive and the output folder from isr_clim/config.py. Each can be
overridden with --site, --sd-dir, --risrn-dir and --out.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import climatology, config, fpi, plots, risrn, waccmxsd

log = logging.getLogger("isr_clim")

CACHE = "waccmxsd_bins"   # Outputs/<CACHE>/<date>.csv


def step_risrn(args) -> int:
    risrn.run(args.risrn_dir, args.out)
    return 0


def step_extract(args) -> int:
    cgm = tuple(sorted(risrn.load_records(args.out)["cgm_lat"].unique()))
    failed = waccmxsd.run(args.sd_dir, args.out / CACHE, cgm, workers=args.workers,
                          start=args.start, end=args.end, overwrite=args.overwrite)
    if failed:
        log.warning("%d days skipped; see %s", failed, args.out / CACHE / "_skipped.csv")
    return 0


def _in_window(times, start, end) -> np.ndarray:
    """Mask of ``times`` from ``start`` to the end of day ``end`` (either may be None)."""
    t = np.asarray(times, "datetime64[ns]")
    keep = np.ones(t.shape, bool)
    if start:
        keep &= t >= np.datetime64(start)
    if end:
        keep &= t < np.datetime64(end) + np.timedelta64(1, "D")
    return keep


def step_climatology(args) -> int:
    station = fpi.station(args.site)
    band = tuple(args.band)
    tag = f"{station.name}_{band[0]}-{band[1]}"
    records = climatology.in_band(risrn.load_records(args.out), band)
    records = records[_in_window(records["time"], args.start, args.end)]
    model = waccmxsd.load_cache(args.out / CACHE, np.sort(records["cgm_lat"].unique()))
    paired = climatology.pair(records, model)
    paired.to_csv(args.out / f"matched_records_{tag}.csv.gz", index=False, float_format="%.6g")
    n_paired = int(paired["model_time"].notna().sum())
    log.info("%s: %d of %d RISR-N records in CGM %s-%s lie within %d min of a 3-hourly WACCM-X SD step "
             "(model %s to %s)", station.name, n_paired, len(paired), *band, fpi.MATCH_WINDOW_MINUTES,
             model["time"].min(), model["time"].max())

    try:
        neutral = fpi.fpi_meridional_zonal(station.name, args.fpi_combo)
    except FileNotFoundError as err:
        log.warning("no FPI comparison: %s", err)
        neutral = None

    label = fpi.MODEL_DISPLAY[config.MODEL_KEY]["label"]
    stats = []
    for comp, (value, error) in config.COMPONENTS.items():
        steps = climatology.obs_mapped_steps(paired, value)
        ion = climatology.seasonal_diurnal(paired["time"], paired[value], paired[error],
                                           paired[f"model_{value}"], *steps, station.tz)
        ion_labels = dict(Instrument="RISR-N", Quantity="ion drift", Model=label, Component=comp)
        ion_stats = climatology.scores(ion, ion_labels)
        stats.append(ion_stats)
        plots.diurnal(ion, comp, station, band, args.out / "plots" / f"ion_diurnal_{tag}_{comp.lower()}.png")
        table = ion.table.set_index(["season", "hour"]).add_prefix("ion_")
        if neutral is not None:
            n = neutral[comp]
            k = _in_window(n["time"], args.start, args.end)
            km = _in_window(n["model_time"], args.start, args.end)
            neu = climatology.seasonal_diurnal(n["time"][k], n["obs"][k], n["err"][k], n["model_at_obs"][k],
                                               n["model_time"][km], n["model_value"][km], station.tz)
            neu_stats = climatology.scores(neu, dict(Instrument="FPI", Quantity="neutral wind", Model=label,
                                                     Component=comp))
            # The ion drift scored only in the (season, local hour) bins the FPI observes in.
            stats += [neu_stats, climatology.scores(ion, ion_labels, hours=climatology.paired_hours(neu),
                                                    hours_label="FPI-observed")]
            plots.residual_comparison(ion, neu, ion_stats, neu_stats, comp, station, band,
                                      args.out / "plots" / f"residual_comparison_{tag}_{comp.lower()}.png")
            table = table.join(neu.table.set_index(["season", "hour"]).add_prefix("fpi_"), how="outer")
        table.reset_index().to_csv(args.out / f"diurnal_{tag}_{comp.lower()}.csv", index=False,
                                   float_format="%.4g")
    pd.concat(stats, ignore_index=True).to_csv(args.out / f"validation_stats_{tag}.csv", index=False)
    log.info("outputs in %s", args.out)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=("risrn", "extract", "climatology", "all"))
    ap.add_argument("--risrn-dir", type=Path, default=config.RISRN_DIR)
    ap.add_argument("--sd-dir", type=Path, default=None,
                    help="WACCM-X SD archive (default: waccmx_sd_dir in the FPI registry)")
    ap.add_argument("--out", type=Path, default=config.OUTPUT_DIR)
    ap.add_argument("--site", default=config.SITE, help="FPI registry station to validate (default: %(default)s)")
    ap.add_argument("--fpi-combo", default=config.FPI_COMBO,
                    help="FPI forward-model run to compare with (default: %(default)s)")
    ap.add_argument("--start", help="first date, YYYY-MM-DD (default: all)")
    ap.add_argument("--end", help="last date, YYYY-MM-DD (default: all)")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--overwrite", action="store_true", help="re-extract days already cached")
    ap.add_argument("--band", type=float, nargs=2, default=config.STATION_BAND,
                    metavar=("LO", "HI"), help="CGM latitude band of RISR-N bins over the station")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", stream=sys.stdout)
    args.sd_dir = args.sd_dir or fpi.sd_dir()
    args.out.mkdir(parents=True, exist_ok=True)
    steps = {"risrn": [step_risrn], "extract": [step_extract], "climatology": [step_climatology],
             "all": [step_risrn, step_extract, step_climatology]}[args.step]
    for f in steps:
        rc = f(args)
        if rc:
            return rc
    return 0


if __name__ == "__main__":
    sys.exit(main())
