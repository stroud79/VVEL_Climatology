# ISR_Climatology

A secondary check of the FPI_ForwardModel findings at Resolute Bay. The
FPI comparison shows WACCM-X SD's thermospheric neutral wind there erring
in a strongly local-time-dependent way (too far south around local
midnight, too far north in the morning, too far west after midnight). At
Resolute the neutral wind is driven largely by ion drag, so this asks
whether WACCM-X SD's ion drift, the convection that drag comes from, errs
in the same way. The ion drift is measured by RISR-N long-pulse vector
velocities (Madrigal kindat 5960).

FPI_ForwardModel is the reference throughout, and its code is imported, not
copied (`isr_clim/fpi.py`):

| What | From FPI_ForwardModel |
|---|---|
| Station coordinates, display name, WACCM-X SD path | `config/sites.yaml` via `fpi_common.registry` (rule 1) |
| Units check, Z3 → geometric height (WGS84 gravity), longitude seam | `FM_V2/src` (`require_units`, `geometric_from_geopotential`, `prepare_lon_interp`) |
| Pairing: nearest model step within 30 min, never interpolated between the 3-hourly steps | `Plotting/climatology_plots.model_at_times`, FM_V2 `TIME_TOLERANCE` |
| Local-time bins (station standard time, no daylight saving, h ± 30 min), 3-hourly display bridging, labels, colours, style | `Plotting/climatology_plots` |
| Inverse-variance weighted means, residuals, scores | the same rules as `plot_seasonal_diurnal_climatology`, applied to both instruments |
| Cluster environment | `config/env.sh` |

`tests/test_core.py::test_fpi_statistics_are_reproduced` checks that this
package, given the FPI data, reproduces the statistics of FPI's current
`plot_seasonal_diurnal_climatology` exactly. Compare against FPI plots
made with that code: `wind_stats_*.csv` files made before it (before
14 Oct 2026) interpolated the 3-hourly model linearly in time.

## Steps

```
python -m isr_clim risrn                 # RISR-N records -> Outputs/risrn_records.csv.gz
python -m isr_clim extract --workers 8   # WACCM-X SD drifts at the bins -> Outputs/waccmxsd_bins/
python -m isr_clim climatology           # pairing, statistics, figures
python -m isr_clim all --workers 8       # all three
```

On BlueBEAR, from this folder: `sbatch slurm/submit_climatology.sh all`
(or a single step; see the script header). FPI_ForwardModel is expected
beside this folder on the RDS share (`ISR_FPI_ROOT` otherwise).

## What each step does

**risrn** (`isr_clim/risrn.py`)
- One file per experiment: every experiment holds the same measurements at
  several integration times, so the first available of 180, 300, 120, 60,
  500, 30 s is used (`Outputs/risrn_files.csv` records which).
- Quality control: both horizontal error estimates < 150 m/s, speed < 3 km/s.
  Records whose error is exactly 1 m/s are degenerate fits (e.g.
  −2537 ± 1 m/s; genuine errors start at ~6 m/s) and are dropped, since
  under inverse-variance weighting a few of them would decide whole
  hourly means.
- Two file layouts, confirmed against the files' own direction angles:
  `vi1`/`vi2` are geographic east/north (`gnangle`, from geographic north);
  `vipe`/`vipn` (91 files, 2016-2018) are local-magnetic (`nangle`), and are
  rotated to geographic with the IGRF declination at the bin, errors
  propagated.
- Records are kept individually, as FPI keeps its observations.

**extract** (`isr_clim/waccmxsd.py`)
- UI and VI (edynamo ion drifts, geographic, m/s), 8 instantaneous steps
  per day at 00, 03, ..., 21 UT, kept as they are.
- Interpolated linearly to each bin's position (bins on the radar's AACGM
  meridian at 300 km), then linearly in geometric height and averaged over
  175-700 km, the layer the radar vectors are fitted over (`altb`/`alte`).
  The SD model top is at ~450-600 km, so the layer stops there
  (`z_top_km`). The drifts change by a few per cent at most over that
  layer.
- One CSV per day in `Outputs/waccmxsd_bins/`; rerunning extracts only
  missing days.

**climatology** (`isr_clim/climatology.py`, `isr_clim/plots.py`), for the
bins of CGM latitude 82.5-83.5° over Resolute:
- Each record is paired with the WACCM-X SD step within 30 min of it, so
  only records near 00, 03, ..., 21 UT are paired. The model line is the
  mean over the paired steps, each once, as in FPI's `obs_mapped` runs.
- `plots/ion_diurnal_<tag>_<component>.png`: FPI's seasonal diurnal
  figure for the ion drift (observed vs model, weighted mean residual,
  residual histogram), same layout and local-time axis as FPI's
  `wind_diurnal_Resolute_Bay_<component>_45elm.png`.
- `plots/residual_comparison_<tag>_<component>.png`: per season, the
  WACCM-X SD residual of the ion drift (vs RISR-N) and of the neutral wind
  (vs FPI) against local time, with each season's bias, r and std ratio.
- `validation_stats_<tag>.csv`: FPI's `wind_stats` columns plus observed
  and model standard deviations and their ratio, for RISR-N (all hours, and
  only the season/local-hour bins the FPI observes) and for the FPI.
- `diurnal_<tag>_<component>.csv`: the hourly values behind the figures.
- `matched_records_<tag>.csv.gz`: every RISR-N record with its model step.

Components are geographic: Meridional (positive north) and Zonal
(positive east), as the FPI horizontal winds at azimuths 0 and 90.

## Tests

`python -m pytest tests`
