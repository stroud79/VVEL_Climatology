# ISR_Climatology — a first-principles explanation of the pipeline

> **Superseded.** This describes the version run as SLURM job 57469373
> (hourly means, linear time interpolation between the 3-hourly model
> outputs, solar/magnetic local time). The package has since been aligned
> with FPI_ForwardModel: records paired with the model step within 30 min
> and never interpolated in time, the station's civil local time, WGS84
> geometric heights, inverse-variance weighted means. See `README.md`;
> the physical background in §4.0 still applies.

This document explains the whole `isr_clim` package in the order in which the
pipeline actually executes: from the RISR-N HDF5 files and WACCM-X SD netCDF
files on disk, through every transformation, to the CSV tables and PNG figures
in `Outputs/`. It explains the physics and mathematics needed to understand each
step, shows how each idea is expressed in the code, and checks the numbers
against the outputs of a real run (SLURM job 57469373, 2026-10-14, whose log is
`slurm-isr_clim-57469373.out` in the repository root).

**How claims are labelled.** Most statements are read directly from the code
and need no label. Where it matters, the following tags are used:

| Tag | Meaning |
|---|---|
| **[code]** | Verified by reading the source in this repository. |
| **[output]** | Verified by recomputing from the files in `Outputs/` (with a stand-alone script; the pipeline itself was not re-run). |
| **[background]** | Standard physics or knowledge about external tools and data sets. It is *not* implemented or recorded in this repository and is given only so the code can be understood. |
| **[interpretation]** | A reasoned inference by the author of this document that the repository does not state. |
| **[uncertain]** | Something that cannot be settled from the repository. |

Line numbers refer to the files as they were when this document was written;
function names are given everywhere so the references survive edits.

---

## Contents

1. [Overview and purpose](#1-overview-and-purpose)
2. [Repository architecture](#2-repository-architecture)
3. [End-to-end pipeline](#3-end-to-end-pipeline)
4. [Detailed step-by-step execution walkthrough](#4-detailed-step-by-step-execution-walkthrough)
   - [4.0 Physical background used throughout](#40-physical-background-used-throughout)
   - [4.1 Stage 0 — entry point, configuration and step dispatch](#41-stage-0--entry-point-configuration-and-step-dispatch)
   - [4.2 Stage 1 — `risrn`: radar archive to hourly bin means](#42-stage-1--risrn-radar-archive-to-hourly-bin-means)
   - [4.3 Stage 2 — `extract`: WACCM-X SD drifts at the radar bins](#43-stage-2--extract-waccm-x-sd-drifts-at-the-radar-bins)
   - [4.4 Stage 3 — `climatology`: matching, binning, scoring, plotting](#44-stage-3--climatology-matching-binning-scoring-plotting)
5. [Final outputs and their interpretation](#5-final-outputs-and-their-interpretation)
6. [Numerical methods, assumptions and limitations](#6-numerical-methods-assumptions-and-limitations)
7. [Validation, testing and reproducibility](#7-validation-testing-and-reproducibility)
8. [Glossary and symbol reference](#8-glossary-and-symbol-reference)
9. [End-to-end worked example](#9-end-to-end-worked-example)
10. [Appendix — discrepancies between documentation and code, and open questions](#10-appendix--discrepancies-between-documentation-and-code-and-open-questions)

---

## 1. Overview and purpose

### 1.1 What the repository does

`isr_clim` builds a **climatology of high-latitude ionospheric ion drift**
(plasma convection velocity) above **Resolute Bay, Nunavut** (geographic
74.73° N, 94.91° W) and uses it to **evaluate a whole-atmosphere model**:

* **Observations** — vector ion velocities measured by **RISR-N** (the Resolute
  Bay Incoherent Scatter Radar – North), as distributed by the Madrigal database
  under kind-of-data code 5960 ("long-pulse vector velocities", per `README.md`).
* **Model** — the ion drifts `UI` (eastward) and `VI` (northward) from
  **WACCM-X SD** (Whole Atmosphere Community Climate Model with thermosphere and
  ionosphere eXtension, run in Specified-Dynamics mode), daily history files
  `f.e20.FXSD.f19_f19.001.cam.h1.<date>-00000.nc`.

The central design idea, stated in `isr_clim/__init__.py` and the module
docstrings, is that **the model is sampled exactly where and when the radar
measured**: at the geographic position of every radar latitude bin and at the
centre time of every radar hourly mean. Only then are both reduced to
season × local-time climatologies and scored. A second, unrestricted model
climatology ("model at all times") is produced alongside to show how much the
radar's intermittent campaign sampling changes the model's own mean.

### 1.2 The problem it solves

An incoherent scatter radar does not run continuously; RISR-N runs in
campaigns ("experiments") of hours to days. A naïve comparison of a radar
climatology with a model climatology built from every model time step mixes
two effects: model error and the radar's uneven sampling in season, local time,
solar cycle and geomagnetic activity. Sampling the model at the radar's times
and positions removes the second effect for the scored comparison.

There are also several practical problems the code solves:

* each radar experiment is stored several times at different integration times
  (cadences) — only one copy may be used or records would be double-counted;
* two different file layouts exist (geographic vs local-magnetic velocity
  components);
* the radar resolves the flow in bins of **magnetic** latitude, whereas the
  model lives on a **geographic** latitude–longitude grid;
* the model is output every 3 h, the radar averages are hourly;
* the model archive is large (3653 daily files), so it must be read only once.

### 1.3 Input-to-output workflow in one paragraph

`python -m isr_clim all` runs three steps. **`risrn`** scans the radar
archive, picks one file per experiment, quality-controls each velocity record,
expresses it in both geographic and local-magnetic east/north components,
averages records into hourly means per magnetic-latitude bin, and attaches each
bin's geographic position, solar local time and magnetic local time
(`Outputs/risrn_hourly.csv.gz`). **`extract`** reads every WACCM-X SD day once,
averages `UI`/`VI` over model levels between 175 and 700 km, interpolates them
bilinearly to the bin positions, rotates them into the magnetic frame, and
caches one CSV per day (`Outputs/waccmxsd_cache/`). **`climatology`**
interpolates the cached model linearly in time to each radar hourly mean,
writes the matched pairs, computes skill scores in a magnetic-latitude band
over Resolute, builds season × local-hour and season × latitude × local-hour
tables, and draws the figures.

### 1.4 Main physical and computational concepts

| Concept | Where it matters |
|---|---|
| E×B drift of F-region ions | What the radar measures and what `UI`/`VI` represent (§4.0) |
| Polar-cap convection (antisunward flow) | Why the diurnal curves look like sinusoids (§4.0, §5) |
| Geographic vs local-magnetic east/north frames; magnetic declination | Frame rotation (§4.2.5) |
| AACGM magnetic coordinates, magnetic local time | Bin placement and MLT binning (§4.2.4, §4.2.7) |
| IGRF geomagnetic field model | Declination at each bin (§4.2.4) |
| Geopotential vs geometric height | Selecting model levels (§4.3.4) |
| Vertical averaging, bilinear spatial interpolation, linear time interpolation | Model sampling (§4.3, §4.4.2) |
| Group-by aggregation (hourly means, climatologies) | §4.2.6, §4.4.4–4.4.6 |
| Verification statistics (bias, RMSE, MAE, Pearson r, std ratio) | §4.4.4 |

---

## 2. Repository architecture

### 2.1 Files

| Path | Role |
|---|---|
| `isr_clim/__main__.py` | Command-line entry point; defines the three steps and their order. |
| `isr_clim/config.py` | Paths, radar location, reference altitude, QC thresholds, model altitude range, time-gap limit, station band, season map. |
| `isr_clim/risrn.py` | Stage 1: scan, select, load and QC radar files; hourly bin means. |
| `isr_clim/geometry.py` | Shared geometry: bin positions (AACGM→geographic), IGRF declination, frame rotations, solar and magnetic local time. |
| `isr_clim/waccmxsd.py` | Stage 2: WACCM-X SD extraction, vertical average, bilinear interpolation, per-day cache. |
| `isr_clim/climatology.py` | Stage 3: time matching, all-times model series, climatology tables, skill scores. |
| `isr_clim/plots.py` | Stage 3: figures. |
| `isr_clim/__init__.py` | Package docstring only. |
| `tests/test_core.py` | Five unit tests of the data-free pieces. |
| `slurm/submit_climatology.sh` | Batch script for the BlueBEAR cluster. |
| `requirements.txt` | Unpinned dependency list. |
| `README.md` | Short user guide. |
| `slurm-isr_clim-57469373.out/.stats` | Log and accounting of the run whose outputs are in `Outputs/`. |
| `Outputs/` | Results of that run (ignored by git, see `.gitignore`). |

### 2.2 Module dependencies

```
__main__ ──► config
   │
   ├──► risrn ───────► config, geometry
   ├──► waccmxsd ────► config, geometry
   ├──► climatology ─► config, geometry
   └──► plots ───────► config
geometry ──► config, aacgmv2, ppigrf
```

`geometry` is the only module that touches magnetic coordinates and the
geomagnetic field, and it is shared by all three stages. This guarantees that
the radar and the model use *identical* bin positions and declinations (the
same function, `geometry.bin_positions`, with the same arguments).

### 2.3 External dependencies (`requirements.txt`)

| Package | Used for |
|---|---|
| `numpy`, `pandas` | All array and table work. |
| `h5py` | Reading RISR-N HDF5 files (`risrn.scan`, `risrn.load_file`). |
| `netCDF4`, `cftime` | Reading WACCM-X SD files and decoding their time axis (`waccmxsd.extract_day`). |
| `aacgmv2` | Geographic ↔ AACGM conversions and magnetic local time (`geometry`). |
| `ppigrf` | IGRF magnetic field vector → declination (`geometry.bin_positions`). |
| `matplotlib` | Figures (`plots`), using the non-interactive `Agg` backend. |
| `pytest` | Tests. |

No versions are pinned (see §7.4 for why that matters). The logged run used
Python 3.10.19 in the conda environment `WACCMX1`.

### 2.4 Entry points

* `python -m isr_clim <step> [options]` → `isr_clim/__main__.py:main`.
* `sbatch slurm/submit_climatology.sh <step> [options]` → finds the project
  folder, activates conda, checks imports, then runs the same command with
  `time` (§4.1.4).

### 2.5 Configuration and data locations

`isr_clim/config.py` computes `PROJECT_ROOT` (the folder holding `isr_clim/`)
and `SHARE_ROOT` (its parent). Defaults:

| Setting | Default | Override |
|---|---|---|
| `RISRN_DIR` | `SHARE_ROOT/RISRN_data` | env `ISR_RISRN_DIR`, CLI `--risrn-dir` |
| `WACCMX_SD_DIR` | `SHARE_ROOT/WACCMX_SD/cam.h1` | env `ISR_WACCMX_SD_DIR`, CLI `--sd-dir` |
| `OUTPUT_DIR` | `PROJECT_ROOT/Outputs` | env `ISR_OUTPUT_DIR`, CLI `--out` |

Environment variables are read once at import time; command-line options take
precedence because their defaults *are* the `config` values.

---

## 3. End-to-end pipeline

### 3.1 Flow diagram

```
                         python -m isr_clim all   (__main__.main)
                                     │
 ┌───────────────────────────────────┼───────────────────────────────────────────┐
 │ STEP 1  risrn.run                 ▼                                           │
 │  RISRN_data/<YYYY>/ran/<exp>/ran*.hdf5                                         │
 │        │ scan()        one row per file: frame, cadence, counts, span        │
 │        ▼                                                                     │
 │  select()          one file per experiment (cadence preference, n_finite)    │
 │        │──────────────────────────────────────► Outputs/risrn_files.csv       │
 │        ▼                                                                     │
 │  load_file()  per chosen file:                                               │
 │     read Table Layout → time = mid-integration, cgm_lat, ve, vn, dve, dvn     │
 │     QC (finite, errors < 150 m/s, speed < 3000 m/s)                          │
 │     geometry.positions_for → bin_positions(year) [AACGM→geo, IGRF decl]      │
 │     geometry.geo_to_mag / mag_to_geo → ve_geo, vn_geo, ve_mag, vn_mag        │
 │        ▼                                                                     │
 │  hourly()   group (hour, cgm_lat): mean, std, n ≥ 3; stamp hh:30;            │
 │             attach glat, glon, mlon, decl, slt, mlt                          │
 │        │──────────────────────────────────────► Outputs/risrn_hourly.csv.gz   │
 └────────┼─────────────────────────────────────────────────────────────────────┘
          │ unique cgm_lat values (23 bins)
 ┌────────▼─────────────────────────────────────────────────────────────────────┐
 │ STEP 2  waccmxsd.run  (ProcessPoolExecutor, one task per day)                │
 │  WACCMX_SD/cam.h1/f.e20.FXSD.f19_f19.001.cam.h1.<date>-00000.nc               │
 │     extract_day():                                                           │
 │       read UI, VI, Z3 north of (southernmost bin − 4°)                        │
 │       Z3 (geopotential) → geometric height z                                 │
 │       mean of UI, VI over levels with 175 ≤ z ≤ 700 km                        │
 │       _bilinear → value at each bin's (glat, glon)                           │
 │       geometry.geo_to_mag with the same declination as the radar             │
 │        │──────────────────────────────────────► Outputs/waccmxsd_cache/<date>.csv
 │        │                                       Outputs/waccmxsd_cache/_skipped.csv
 └────────┼─────────────────────────────────────────────────────────────────────┘
          │ load_cache(): all days concatenated
 ┌────────▼─────────────────────────────────────────────────────────────────────┐
 │ STEP 3  __main__.step_climatology                                            │
 │  climatology.match(obs, model)  linear time interpolation, ≤ 3.5 h brackets  │
 │        │──────────────────────────────────────► Outputs/matched_hourly.csv.gz │
 │  climatology.model_all_times(model)  model at every hh:30, with slt, mlt     │
 │  climatology.scores(matched, band)   ─────────► stats_waccmxsd_<band>.csv     │
 │  for component in --components, for lt in --lt:                              │
 │     climatology.diurnal  ──────────────────────► diurnal_<c>_<lt>_<band>.csv  │
 │     plots.diurnal        ──────────────────────► plots/diurnal_…png           │
 │     climatology.lat_lt   ──────────────────────► lat_lt_<c>_<lt>.csv          │
 │     plots.lat_lt         ──────────────────────► plots/lat_lt_…png            │
 └──────────────────────────────────────────────────────────────────────────────┘
```

### 3.2 Stage summary

| Stage | Function | Consumes | Produces | Why the next stage needs it |
|---|---|---|---|---|
| 0 | `__main__.main` | CLI args, env, `config` | `args`, `Outputs/` folder | Supplies paths, band, components, worker count. |
| 1 | `risrn.run` | RISR-N HDF5 archive | `risrn_files.csv`, `risrn_hourly.csv.gz` | The set of bins defines where the model is sampled; the hourly times define when. |
| 2 | `waccmxsd.run` | netCDF archive, bin list from stage 1 | `waccmxsd_cache/<date>.csv` | A compact time series of model drift at each bin. |
| 3 | `step_climatology` | Stage 1 hourly means, stage 2 cache | Matched table, scores, climatologies, plots | Final products. |

Stage 2 depends on stage 1 only through the **list of bin latitudes**
(`step_extract`, `__main__.py:40`). Stage 3 depends on both outputs.

Logged timings of the reference run: stage 1 ≈ 52 s; stage 2 ≈ 11.5 min for
3653 files on 8 workers; stage 3 ≈ 1.5 min; total 13 min 57 s.

---

## 4. Detailed step-by-step execution walkthrough

### 4.0 Physical background used throughout

This section introduces the physics once so that later sections can refer
back to it. Everything here is **[background]**: the repository relies on these
ideas but does not implement them.

#### 4.0.1 The F-region ionosphere and E×B drift

Above about 150 km the upper atmosphere is partly ionised (the ionosphere); the
F region (≈150–600 km) is dominated by O⁺ ions and electrons. Consider one ion
species of mass $m_i$ and charge $q$ moving with velocity $\mathbf v$ through
neutral gas moving with velocity $\mathbf u$, in electric field $\mathbf E$ and
magnetic field $\mathbf B$. Its momentum equation (neglecting pressure
gradients and gravity, which matter mainly along $\mathbf B$) is

$$
m_i\frac{d\mathbf v}{dt} = q\left(\mathbf E + \mathbf v\times\mathbf B\right) - m_i\nu_{in}\left(\mathbf v-\mathbf u\right),
$$

where $\nu_{in}$ (s⁻¹) is the ion–neutral collision frequency. Define the ion
gyrofrequency $\Omega_i = qB/m_i$ (rad s⁻¹). In steady state and when
$\nu_{in}\ll\Omega_i$ the collision term is negligible and
$\mathbf E + \mathbf v\times\mathbf B = 0$ for the components perpendicular to
$\mathbf B$. Taking the cross product with $\mathbf B$:

$$
\mathbf E\times\mathbf B = -(\mathbf v\times\mathbf B)\times\mathbf B = B^2\mathbf v_\perp - (\mathbf v\cdot\mathbf B)\mathbf B
\quad\Rightarrow\quad
\boxed{\mathbf v_\perp = \frac{\mathbf E\times\mathbf B}{B^2}}
$$

Units: (V m⁻¹)(T)/(T²) = V m⁻¹ T⁻¹ = m s⁻¹. ✔

Order of magnitude: for O⁺ ($m_i = 16\times1.66\times10^{-27}$ kg) in
$B\approx5\times10^{-5}$ T, $\Omega_i \approx 300$ s⁻¹, while $\nu_{in}$ at
300 km is of order 0.1–1 s⁻¹, so $\nu_{in}/\Omega_i\sim10^{-3}$ and the drift
is pure E×B. The two become comparable near 120–130 km. An electric field of
20 mV m⁻¹ gives $|\mathbf v_\perp| = 0.02/5\times10^{-5} = 400$ m s⁻¹, the
same order as the speeds in the outputs.

Two consequences used implicitly by the pipeline:

1. **The drift is nearly independent of height in the F region.** Magnetic
   field lines are very good conductors, so the electric potential is almost
   constant along them, and $\mathbf E_\perp$ (hence $\mathbf v_\perp$) changes
   only through field-line geometry. Averaging the model drift over
   175–700 km (§4.3.5) therefore averages values that should be similar, and
   the radar's single "F-region" vector is comparable with that average.
2. **At Resolute the field is almost vertical.** The magnetic inclination is
   ≈ 87–88° (the risrn docstring says "~87 deg"), so the perpendicular plane
   is almost horizontal and the horizontal components of $\mathbf v_\perp$ are
   essentially the full drift. A field line tilted by $90^\circ-I\approx2$–$3^\circ$
   from the vertical moves only ≈ 17–26 km horizontally per 500 km of height,
   so placing all bins at one reference altitude (300 km, §4.2.4) is benign.

#### 4.0.2 Polar-cap convection

The solar wind's interaction with Earth's magnetosphere (the Dungey cycle)
imposes an electric potential across the polar cap, producing a two-cell
convection pattern: flow is **antisunward** (from the dayside to the
nightside) across the central polar cap and returns sunward at lower latitudes
in the auroral zone. Resolute Bay sits at AACGM latitude ≈ 82.2–82.8° at
300 km (comment above `config.STATION_BAND`), i.e. normally inside the polar
cap.

A simple model shows why the outputs look like sinusoids in local time. Look
down on the pole with the Sun in the $+x$ direction. A station at local time
$\tau$ lies at azimuth $\varphi = 2\pi(\tau-12\,\mathrm h)/24\,\mathrm h$ from the
noon meridian, and its local "northward" unit vector points toward the pole:
$\hat{\mathbf n} = -(\cos\varphi, \sin\varphi)$. For uniform antisunward flow
$\mathbf v = -V\hat{\mathbf x}$,

$$
v_n = \mathbf v\cdot\hat{\mathbf n} = V\cos\varphi ,
$$

which is $+V$ (northward) at local noon, $-V$ (southward) at midnight, and zero
at 06 and 18 local time. The observed diurnal curves (§5.2) have this
shape with $V\approx300$ m s⁻¹, with zero crossings shifted (≈ 04:30–05:30
and 17:00–18:30 solar local time in the band), as expected because the real pattern is
organised about the *magnetic* pole and skewed by the interplanetary magnetic
field **[interpretation]**.

#### 4.0.3 How the radar obtains a vector velocity

An incoherent scatter radar transmits pulses and receives the weak echo
scattered from thermal fluctuations of the plasma. The Doppler shift of the
echo gives the **line-of-sight** ion velocity in each range gate of each
beam. RISR-N is a phased array that interleaves many beam directions. A vector
velocity is obtained by assuming that, within one bin of magnetic latitude,
the perpendicular drift is uniform, and solving the overdetermined set of
line-of-sight equations (with error estimates) for the vector. The product
used here gives such vectors per **corrected geomagnetic latitude** bin
(`cgm_lat`), integrated over a fixed time (`inttms`). **This inversion is done
by the data provider before the files reach the repository**; the code starts
from the finished vectors and their error estimates `dvi1`, `dvi2` /
`dvipe`, `dvipn`.

#### 4.0.4 WACCM-X SD and its ion drift

WACCM-X is a configuration of the NCAR Community Earth System Model that
extends from the ground to the upper thermosphere (≈ 500–700 km depending on
solar activity). In **SD** ("specified dynamics") mode the lower and middle
atmosphere are nudged towards a meteorological reanalysis, so the model
follows the real weather of each day. Its electrodynamics module ("edynamo",
named in the `waccmxsd` docstring) computes electric fields and the resulting
E×B ion drifts, output as `UI` (zonal) and `VI` (meridional) in geographic
components. At high latitudes the convection electric field in WACCM-X is
normally imposed by an empirical model driven by solar-wind/geomagnetic
indices **[background; which one this run used is [uncertain] — it is not
recorded in the repository]**. `f19_f19` denotes the ≈ 1.9° latitude × 2.5°
longitude finite-volume grid. The `h1` files hold 8 instantaneous outputs per
day at 3-hour spacing (docstring of `waccmxsd`; confirmed **[output]**: 3651 of
3653 cached days have exactly 8 times 00, 03, …, 21 UT).

---

### 4.1 Stage 0 — entry point, configuration and step dispatch

**Purpose.** Parse the command line, set up logging and the output folder,
and run the requested steps in order.

**Inputs.** Command-line arguments, environment variables
(`ISR_RISRN_DIR`, `ISR_WACCMX_SD_DIR`, `ISR_OUTPUT_DIR`, `SLURM_CPUS_PER_TASK`),
and constants in `config.py`.

#### 4.1.1 Argument parsing (`__main__.main`, lines 79–98)

| Option | Default | Used by | Effect |
|---|---|---|---|
| `step` | (required) `risrn`/`extract`/`climatology`/`all` | dispatcher | Which steps run. |
| `--risrn-dir` | `config.RISRN_DIR` | stage 1 | Radar archive root. |
| `--sd-dir` | `config.WACCMX_SD_DIR` | stage 2 | Model archive folder. |
| `--out` | `config.OUTPUT_DIR` | all | Output folder (created if missing, line 101). |
| `--start`, `--end` | none | stages 2, 3 | Inclusive date range `YYYY-MM-DD`. **Not used by stage 1**, which always processes the whole radar archive. |
| `--workers` | `$SLURM_CPUS_PER_TASK` or 1 | stage 2 | Process-pool size. |
| `--overwrite` | off | stage 2 | Re-extract days already in the cache. |
| `--band LO HI` | `config.STATION_BAND` = (82.5, 83.5) | stage 3 | Magnetic-latitude band for scores and diurnal tables. |
| `--components` | `vn_geo vn_mag` | stage 3 | Which components get diurnal and lat–LT tables/plots (scores always include all four). |
| `--lt` | `slt mlt` | stage 3 | Local-time coordinate(s) to bin by. |
| `--lt-start` | 12 | stage 3 plots | Left edge of the SLT diurnal plot axis (noon-to-noon). MLT plots always start at 0 (line 71). |
| `-v` | off | logging | DEBUG level. |

#### 4.1.2 Dispatch (lines 99–108)

Logging goes to stdout with timestamps. `steps` maps the chosen step to a list
of functions; `all` is `[step_risrn, step_extract, step_climatology]`. Each
returns 0 on success; a non-zero return would stop the loop, though in practice
failures surface as exceptions, which end the program with a traceback and a
non-zero exit status.

#### 4.1.3 Configuration constants (`config.py`)

| Constant | Value | Units | Consumed in |
|---|---|---|---|
| `RADAR_LAT`, `RADAR_LON` | 74.72955, 265.09424 | deg (geographic; lon east 0–360) | `geometry.bin_positions` |
| `REF_ALT_KM` | 300.0 | km | AACGM conversions and IGRF evaluation |
| `PREFERRED_CADENCE_S` | (180, 300, 120, 60, 500, 30) | s | `risrn.select` |
| `MAX_VELOCITY_ERROR` | 150.0 | m s⁻¹ | `risrn.load_file` QC |
| `MAX_SPEED` | 3000.0 | m s⁻¹ | `risrn.load_file` QC |
| `MIN_SAMPLES_PER_HOUR` | 3 | records | `risrn.hourly` |
| `MODEL_ALT_RANGE_KM` | (175.0, 700.0) | km | `waccmxsd.extract_day` |
| `MODEL_MAX_GAP_H` | 3.5 | h | `climatology.match` |
| `STATION_BAND` | (82.5, 83.5) | deg AACGM lat | CLI default for `--band` |
| `SEASONS`, `SEASON_ORDER` | DJF/MAM/JJA/SON by calendar month | — | `climatology.add_time_bins`, `scores`, plots |

The comments in `config.py` give the reasons: 175–700 km is "the range the
RISR-N vector velocities are resolved over"; the band 82.5–83.5° brackets
Resolute's AACGM latitude (≈ 82.2–82.8° at 300 km in 2009–2018) and starts at
82.5 because "the files' lowest bin, 82.25 deg, never holds a vector"
(consistent with the outputs: the lowest bin present in `risrn_hourly` is 82.5
**[output]**).

#### 4.1.4 The SLURM wrapper (`slurm/submit_climatology.sh`)

Requests 1 node, 8 CPUs, 32 GB, 12 h. It walks up from the submit directory to
find the folder containing `isr_clim/` (lines 33–44), activates conda
(`ISR_CONDA_BASE`, `ISR_CONDA_ENV`, defaulting to the author's miniconda and
`WACCMX1`), **imports every dependency and aborts with the real import error
if any fails** (lines 64–76), prints a header and runs
`time python -m isr_clim <step> <args>`. Because `SLURM_CPUS_PER_TASK=8`,
`--workers` defaults to 8. The accounting file shows the reference run used
1.2 GB of the 32 GB requested and 31 % CPU efficiency (stages 1 and 3 are
single-process).

**Outputs of stage 0.** The parsed `args` namespace and an existing `--out`
folder.

---

### 4.2 Stage 1 — `risrn`: radar archive to hourly bin means

Entry: `__main__.step_risrn` → `risrn.run(args.risrn_dir, args.out)`
(`risrn.py:134–157`). `run` calls, in order, `scan` → `select` → `load_file`
(per chosen file) → `hourly`, then writes two files.

#### 4.2.1 Inputs: the archive and the file format

* **Layout on disk** (docstring): `ROOT/<YYYY>/ran/<experiment>/ran*.hdf5`,
  as produced by Madrigal's `globalDownload.py --tree`. `<experiment>` is a
  Madrigal experiment name such as `24aug09b`.
* **Inside each file**: a compound (record) array at
  `Data/Table Layout`. Each row is one **(integration period, magnetic
  latitude bin)** pair. Fields used:

| Field | Meaning | Units |
|---|---|---|
| `ut1_unix`, `ut2_unix` | Start and end of the integration period | s since 1970-01-01 UTC |
| `inttms` | Integration time | s |
| `cgm_lat` | Corrected geomagnetic latitude of the bin | deg |
| `vi1`, `vi2` | Eastward, northward ion velocity, **geographic** frame ("geo" layout) | m s⁻¹ |
| `dvi1`, `dvi2` | Their error estimates | m s⁻¹ |
| `vipe`, `vipn` | Perpendicular-east, perpendicular-north ion velocity, **local magnetic** frame ("mag" layout) | m s⁻¹ |
| `dvipe`, `dvipn` | Their error estimates | m s⁻¹ |

The archive used for the reference run contained 615 files in 290 experiments
(2009–2018), all readable **[output, log]**.

#### 4.2.2 `scan(root)` — inventory every file (`risrn.py:50–69`)

**Purpose.** Learn, cheaply, what each file contains so that one file per
experiment can be chosen.

**Step by step.** For each path matching `*/ran/*/*.hdf5` (sorted):

1. Open the table, get the field names, and decide the layout with
   `_layout(names)` (lines 43–47): `"geo"` if all of `vi1, vi2, dvi1, dvi2`
   exist, else `"mag"` if all of `vipe, vipn, dvipe, dvipn` exist, else
   `ValueError`.
2. Read only five columns: `ut1_unix`, `ut2_unix`, `inttms`, and the east and
   north velocity.
3. Record: `experiment` (the parent folder name), `frame`, `cadence` =
   `round(median(inttms))`, `n_records`, `n_finite` (rows where *both*
   components are finite), `t0`/`t1` (earliest start, latest end), `error=""`.
4. Any exception (corrupt file, unknown layout) is caught and stored as a row
   with `error=repr(err)` — **reported, not fatal**.

**Output.** A DataFrame, one row per file.

#### 4.2.3 `select(files)` — one file per experiment (`risrn.py:72–78`)

**Why.** The docstring explains that each experiment holds *the same
measurements* several times, integrated over different times (e.g. 60, 180,
300, 500 s). Using more than one would count the same plasma several times.

**Algorithm.**

1. Keep rows with `error == ""`.
2. `rank` = position of the file's cadence in `PREFERRED_CADENCE_S`
   (180 → 0, 300 → 1, 120 → 2, 60 → 3, 500 → 4, 30 → 5); any other cadence gets
   rank 6 (still usable, least preferred).
3. Sort by `experiment` ↑, `rank` ↑, `n_finite` ↓ and take the first row per
   experiment (`groupby(...).first()`; since all columns of readable rows are
   non-null this is the first row).

So the preferred cadence wins; among files of the same cadence the one with
more finite vectors wins (the docstring's stated tie-break).

**Design rationale.** The code states the preference order but not why 180 s
comes first. **[interpretation]** A 3-minute integration is a compromise:
long enough for low statistical error, short enough that ≥ 3 samples fit in an
hour many times over, and it is the most common cadence in the archive (275 of
615 files; 266 of the 290 chosen files are 180 s, 16 are 60 s, 6 are 300 s,
2 are 120 s **[output]**).

`run` then marks `files["used"]`, and later fills `n_good` for used files, and
writes everything to **`Outputs/risrn_files.csv`** (columns `path, experiment,
frame, cadence, n_records, n_finite, t0, t1, error, used, n_good`). This is
the audit trail of which file represented each experiment.

#### 4.2.4 `load_file(path)` — records, quality control, positions (`risrn.py:81–107`)

For each chosen file, in `select` order:

**(a) Read and tabulate** (lines 83–92). The whole table is read. A DataFrame
is built with:

* `time` = midpoint of the integration period,
  $t = \tfrac12(t_{\mathrm{ut1}}+t_{\mathrm{ut2}})$, converted to a naïve UTC
  timestamp;
* `cgm_lat` rounded to 4 decimals (so it can be used as an exact join key);
* `ve`, `vn`, `dve`, `dvn` from the layout's four columns, as float64.

**(b) Quality control** (lines 93–97). A record is kept only if

$$
\text{finite}(v_e, v_n, \delta v_e, \delta v_n)\;\wedge\;
\delta v_e < 150\ \mathrm{m\,s^{-1}}\;\wedge\;
\delta v_n < 150\ \mathrm{m\,s^{-1}}\;\wedge\;
\sqrt{v_e^2+v_n^2} < 3000\ \mathrm{m\,s^{-1}} .
$$

The error threshold removes poorly constrained vectors (low signal, poor beam
geometry); the speed threshold removes physically implausible ones
(polar-cap drifts rarely exceed ~2 km s⁻¹ **[background]**). The speed test is
applied in the file's own frame; since the rotation in (d) preserves length,
the result is frame-independent. Note that the error thresholds are applied to
the *file's* components (geographic for "geo" files, magnetic for "mag"
files) — the two are not exactly equivalent tests, though for comparable
errors the difference is small **[interpretation]**. In the used files
1,856,573 of 3,071,089 finite vectors (60.5 %) pass **[output]**. The counts
`n_finite` and `n_good` are attached to `df.attrs`.

**(c) Bin positions** — `geometry.positions_for(df["time"], df["cgm_lat"])`.
This is the first use of the geometry module and it is the same machinery the
model extraction uses, so it is explained here in full.

*First principles: AACGM coordinates.* **[background]** Earth's field is
roughly, but not exactly, a tilted dipole. *Altitude-Adjusted Corrected
Geomagnetic* (AACGM) coordinates label a point by the dipole field line that
reaches the same apex as the real (IGRF) field line through the point: trace
the real field line from the point to its apex at radius $r_a$; the AACGM
latitude $\lambda$ satisfies $\cos^2\lambda = (R_E+h)/r_a$ for a point at
altitude $h$, and the AACGM longitude is the dipole longitude of the apex.
Points on one field line share coordinates, which is why the convection
pattern (which maps along field lines) is ordered by them. The `aacgmv2`
library evaluates this with fitted coefficients, and the transformation drifts
from year to year as the internal field changes.

*Implementation* — `geometry.bin_positions(year, cgm_lats)` (lines 23–37),
memoised with `lru_cache` so each (year, bin set) is computed once per
process:

1. `when = 1 July of year` — one epoch per year (mid-year).
2. Radar AACGM longitude: `aacgmv2.convert_latlon(RADAR_LAT, RADAR_LON, 300 km,
   when, "G2A")` → `radar_mlon`.
3. Each bin is placed on that **magnetic meridian**:
   `convert_latlon_arr(cgm_lats, radar_mlon, 300 km, when, "A2G")` →
   geographic `glat`, `glon`. Longitude is wrapped to (−180, 180]:
   `(glon + 180) % 360 − 180`.
4. IGRF field at each bin at 300 km: `ppigrf.igrf(glon, glat, 300, when)` →
   $B_e, B_n, B_u$ (nT; east, north, up). The **declination** — the angle of
   the horizontal field east of geographic north — is

   $$
   D = \operatorname{atan2}(B_e, B_n)\quad\text{(degrees, east positive)}.
   $$

5. Return a table `cgm_lat, glat, glon, mlon (=radar_mlon), decl`.

Real values for 2012 **[output]**: bin 83.0 → (75.14° N, −94.67°), decl
−20.2°; bin 85.0 → (77.35° N, −92.88°), decl −31.6°; and for 2013 (the 88° bin
has no 2012 data) bin 88.0 → (80.73° N, −88.73°), decl −50.7°. `radar_mlon` drifts from −36.05° (2009) to
−32.57° (2018). The large, rapidly varying declination reflects proximity to
the north magnetic pole; using the *bin's* declination rather than the radar's
is important for the higher bins.

*Design rationale.* The radar's product is resolved in magnetic latitude only,
so a longitude must be chosen to give each bin a geographic location. The code
chooses the radar's own AACGM meridian (module docstring of `geometry`): the
bins of a near-zenith/meridional radar field of view lie approximately along
it, and it makes the bins' positions a function of `cgm_lat` and year only
**[the "approximately" is an interpretation; the actual beam footprints are
not in the repository]**.

`positions_for(times, cgm_lat)` (lines 40–49) maps every row to its year's
table: for each distinct year it calls `bin_positions(year, all_bins_present)`
and re-indexes by the row's rounded `cgm_lat`, then restores the original row
order. Output: a DataFrame aligned to the input rows with `glat, glon, mlon,
decl`.

**(d) Frame rotation** (lines 99–104). Each record is expressed in both frames
using that row's declination.

*Derivation.* Let $\hat{\mathbf e}, \hat{\mathbf n}$ be geographic east and
north unit vectors. Magnetic north lies at azimuth $D$ east of geographic north,
so the local magnetic unit vectors are

$$
\hat{\mathbf n}_m = \sin D\,\hat{\mathbf e} + \cos D\,\hat{\mathbf n},\qquad
\hat{\mathbf e}_m = \cos D\,\hat{\mathbf e} - \sin D\,\hat{\mathbf n}.
$$

Writing $\mathbf v = v_{e,m}\hat{\mathbf e}_m + v_{n,m}\hat{\mathbf n}_m$ and
collecting terms gives **magnetic → geographic** (`geometry.mag_to_geo`, lines
52–57):

$$
v_{e,g} = v_{e,m}\cos D + v_{n,m}\sin D,\qquad
v_{n,g} = v_{n,m}\cos D - v_{e,m}\sin D .
$$

Its inverse (the transpose of the rotation matrix) is **geographic →
magnetic** (`geometry.geo_to_mag`, lines 60–63):

$$
v_{e,m} = v_{e,g}\cos D - v_{n,g}\sin D,\qquad
v_{n,m} = v_{n,g}\cos D + v_{e,g}\sin D .
$$

This is a pure rotation: it preserves speed, and with $D<0$ (west, as at
Resolute) magnetic north points west of geographic north.

* "geo" files: `ve_geo, vn_geo` = file values; `ve_mag, vn_mag` from
  `geo_to_mag`.
* "mag" files: `ve_mag, vn_mag` = file values; `ve_geo, vn_geo` from
  `mag_to_geo`.

*Assumption* (stated in the `risrn` docstring): the "mag" layout's components
are perpendicular to $\mathbf B$, not horizontal. The perpendicular-north
direction is tilted from horizontal by $90^\circ - I$, so its horizontal
projection is $v_{pn}\sin I$; the neglected fraction is $1-\sin I$. The
docstring says this is "under 0.1%" at $I\approx87^\circ$; strictly
$1-\sin87^\circ = 0.14\,\%$, and it falls below 0.1 % for $I>87.4^\circ$.
Either way it is far below the measurement error **[interpretation]**.

*Why both frames?* Geographic components are what the model outputs and what
geographic local time refers to; magnetic components are the natural frame for
convection organised by the geomagnetic field. Carrying both lets every
downstream table be produced in either frame without re-reading data.

In the reference run 236 used files were "geo" and 54 "mag" (all 2016–2018;
91 "mag" files exist in total, most of them alternative cadences)
**[output]**.

**Outputs of `load_file`.** Columns `time, cgm_lat, dve, dvn, ve_geo, vn_geo,
ve_mag, vn_mag, frame`; `run` adds `experiment` and concatenates all files into
`records` (1,856,573 rows in the reference run).

#### 4.2.5 Why the rotation is done before averaging

Because $D$ is constant for a given (bin, year), and rotation is linear, the
mean of rotated records equals the rotation of the mean records. Doing it
record by record is nevertheless necessary for the **standard deviations** of
the converted frame (they are not related by a simple rotation of the two
standard deviations) and for mixing "geo" and "mag" files in one hourly bin.

#### 4.2.6 `hourly(records)` — hourly bin means (`risrn.py:113–131`)

**Purpose.** Reduce irregular 0.5–8-min records to a regular hourly series per
bin, matching the time resolution at which the model will be compared.

**Step by step.**

1. `hour` = record time floored to the hour.
2. Group by `(hour, cgm_lat)` and aggregate: mean and sample standard
   deviation (pandas, $N-1$ denominator) of each of `vn_geo, ve_geo, vn_mag,
   ve_mag`; mean of `dvn` and `dve`; the first `experiment`; and the record
   count.

   $$
   \bar v = \frac1n\sum_{k=1}^{n} v_k,\qquad s_v = \sqrt{\frac{1}{n-1}\sum_{k=1}^{n}(v_k-\bar v)^2}.
   $$

3. Flatten column names: means keep the plain name (`vn_geo`), standard
   deviations become `vn_geo_std`, the count becomes `n`.
4. Drop groups with `n < MIN_SAMPLES_PER_HOUR` (3).
5. Stamp each mean at the **centre of the hour**: `time = hour + 30 min`.
6. Attach `glat, glon, mlon, decl` from `positions_for` at the new times.
7. Attach solar and magnetic local time (below).

Notes:

* The mean is unweighted (not error-weighted) and `dvn`/`dve` are the mean of
  the per-record error estimates, **not** the standard error of the hourly
  mean.
* The hh:30 stamp assumes the records are spread across the hour; an hour
  with records only in its first 10 minutes is still stamped at hh:30
  **[interpretation]**.
* At 180 s cadence an hour holds at most 20 records; the distribution of `n`
  peaks at 16–19, and values up to 59 come from 60 s experiments
  **[output]**.

#### 4.2.7 Local times (`geometry.solar_local_time`, `geometry.magnetic_local_time`)

*Solar local time* (lines 66–70). Earth rotates 360° in 24 h, i.e. 15° h⁻¹,
so the mean solar time at longitude $\phi$ (deg east) is

$$
\mathrm{SLT} = \left(\mathrm{UT} + \frac{\phi}{15^\circ\,\mathrm h^{-1}}\right)\bmod 24\ \mathrm h .
$$

The code uses the bin's `glon`, not the radar's. This is *mean* solar time:
the equation of time (± 16 min) is not applied **[code]**.

*Magnetic local time* (lines 73–81). **[background]** MLT is the angle,
in hours, between the point's magnetic longitude and the magnetic longitude of
the subsolar point: $\mathrm{MLT} = 12 + (\Lambda - \Lambda_{\odot})/15$ (mod 24),
where $\Lambda_\odot$ depends on time. `aacgmv2.convert_mlt(mlon, time,
m2a=False)` evaluates it. The code groups rows by distinct time so the
expensive solar computation is done once per time. Because every bin has the
same `mlon` (the radar's), all bins at one time share the same MLT.

**Outputs of stage 1.** `run` writes:

* **`Outputs/risrn_files.csv`** — see §4.2.3.
* **`Outputs/risrn_hourly.csv.gz`** (`%.6g` floats), one row per (hour, bin):

| Column | Meaning | Units |
|---|---|---|
| `cgm_lat` | AACGM latitude of the bin | deg |
| `vn_geo`, `ve_geo`, `vn_mag`, `ve_mag` | Hourly mean drift components | m s⁻¹ |
| `*_std` | Sample std of the records in the hour | m s⁻¹ |
| `dvn`, `dve` | Mean reported error of the file's north/east components | m s⁻¹ |
| `experiment` | Experiment that supplied the (first) record | — |
| `n` | Records averaged (≥ 3) | — |
| `time` | Hour centre, UTC | — |
| `glat`, `glon` | Geographic position of the bin at 300 km | deg |
| `mlon` | AACGM longitude of the radar (and so of every bin) | deg |
| `decl` | IGRF declination at the bin | deg, east positive |
| `slt`, `mlt` | Solar and magnetic local time | h |

The reference run produced **112,081** hourly bin means from 2009-08-24 to
2018-10-11, in 23 bins from 82.5° to 88.0° every 0.25° **[output, log]**.
Bin populations are very uneven: 20 rows at 82.5°, 1,618 at 82.75°, ≈ 8,000
at 83.5–84.5°, 25 at 88.0°.

**Error handling in stage 1.** Unreadable files are listed, not fatal;
an empty archive raises `FileNotFoundError` (`run`, line 141); `load_file`
does not catch errors, so a file that `scan` could read but `load_file` cannot
(e.g. missing `cgm_lat`) would stop the step **[code]**.

---

### 4.3 Stage 2 — `extract`: WACCM-X SD drifts at the radar bins

Entry: `__main__.step_extract` (lines 39–45). It reads
`risrn_hourly.csv.gz` (raising `FileNotFoundError` with a hint if stage 1 has
not run), takes the **sorted unique `cgm_lat` values** as a tuple and calls
`waccmxsd.run(sd_dir, out/"waccmxsd_cache", cgm, workers, start, end,
overwrite)`.

#### 4.3.1 Purpose

Turn ≈ 3 653 large 3-D netCDF files into a small table "model drift at each
radar bin, every 3 h", read the archive only once, and make the step
resumable.

#### 4.3.2 Inputs

Daily files whose names match
`f\.e20\.FXSD\.f19_f19\.001\.cam\.h1\.(\d{4}-\d{2}-\d{2})-00000\.nc$`
(`FILE_RE`, line 33). Required variables (`NEEDED`, line 35):

| Variable | Dimensions | Meaning | Units |
|---|---|---|---|
| `time` | (time) | Output times, CF-encoded | e.g. days since … |
| `lev` | (lev) | Model level (hybrid pressure) | — |
| `lat`, `lon` | (lat), (lon) | Grid, lat ascending, lon 0–360 | deg |
| `UI`, `VI` | (time, lev, lat, lon) | Zonal / meridional ion drift (geographic) | m s⁻¹ (cm s⁻¹ handled) |
| `Z3` | (time, lev, lat, lon) | Geopotential height | m |

`lev` is checked for presence but not otherwise used **[code]**.

#### 4.3.3 Orchestration — `run` (lines 118–137)

1. `list_files(root, start, end)` (lines 38–44): glob, regex-match the date,
   and keep dates within `[start, end]` by *string* comparison (valid because
   the format is ISO `YYYY-MM-DD`).
2. `todo` = days whose `cache/<date>.csv` does not exist (all days with
   `--overwrite`). **This is the resume mechanism**: the SLURM header notes a
   run that hits the time limit can be resubmitted unchanged.
3. A `ProcessPoolExecutor(max_workers=workers)` runs `_cache_one(date, path,
   cache, cgm_lats, alt_range)` for each day. Separate processes sidestep
   Python's GIL and the netCDF library's thread-safety limits.
4. Results are collected with `as_completed`. **Any exception in a day is
   caught** and appended to `skipped` with its `repr` — "one bad day must not
   stop the rest". Progress is logged every 50 days.
5. `cache/_skipped.csv` (`date, path, reason`) is **rewritten** with this run's
   failures; the number of failures is returned and `step_extract` logs a
   warning if it is non-zero.

`_cache_one` (lines 110–115) writes to a hidden `.<date>.csv.tmp` and then
atomically renames it to `<date>.csv`, so an interrupted job never leaves a
half-written day that the resume logic would mistake for a finished one.

#### 4.3.4 `extract_day(path, cgm_lats)` — reading and height conversion (lines 67–88)

**(a) Checks and units** (lines 69–76). Missing variables → `KeyError("no
[...]")` (these are the "files without UI/VI" the docstring mentions; they end
up in `_skipped.csv`). If `UI`/`VI` declare units `cm/s` or `cm s-1` they are
multiplied by 0.01; anything else (including no units attribute) is taken as
m s⁻¹.

**(b) Times** (lines 77–79). `cftime.num2date` decodes the CF time axis using
the file's calendar (default `standard`), then each date is converted to ISO
text and parsed by pandas. Going through `cftime` makes non-standard model
calendars decodable. **[output]** The cache contains 29 February in 2008, 2012
and 2016, so this run uses a calendar with leap days.

**(c) Horizontal subset** (lines 80–87). Bin positions for this file's year
are computed with the same `geometry.bin_positions` used for the radar. Only
grid latitudes $\ge \min(\text{glat}) - 4^\circ$ are read (≈ 71° N and
northward), because the bins lie between ≈ 75° N and 81° N. The 4° margin
guarantees at least one grid row (spacing ≈ 1.9°) south of the southernmost
bin for interpolation, while reading only ≈ 10 % of the global field. All
longitudes are kept. `_as_float` (lines 47–49) converts masked arrays to
float64 with masked points → NaN. `Z3` is converted from m to km.

**(d) Geopotential height → geometric height** (line 88).

*First principles.* Geopotential $\Phi(z)$ is the work per unit mass needed to
lift a parcel from sea level to height $z$ against gravity. For a spherical
Earth of radius $R$, gravity falls off as
$g(z) = g_0\,R^2/(R+z)^2$, so

$$
\Phi(z) = \int_0^z g(z')\,dz' = g_0R^2\left[\frac1R-\frac1{R+z}\right] = \frac{g_0Rz}{R+z}.
$$

Geopotential height is $Z = \Phi/g_0 = Rz/(R+z)$. Solving for $z$:

$$
\boxed{z = \frac{R\,Z}{R - Z}}
\qquad R = 6371\ \text{km (`EARTH_RADIUS_KM`)}.
$$

Units: km·km/km = km ✔. The difference is substantial at these heights:
$Z=300$ km → $z=314.8$ km; $Z=600$ km → $z=662.4$ km. Without the
conversion the 175–700 km window would sit ≈ 5–10 % too low.

*Assumptions.* Spherical Earth, $R$ = 6371 km, and that the model's `Z3` is a
geopotential height defined with constant reference gravity $g_0$ relative to
sea level. **[uncertain]** Whether WACCM-X's `Z3` is defined exactly this way
(and above sea level rather than the surface) is not documented in the
repository; any discrepancy would shift the window by a few km at most.

#### 4.3.5 Vertical average over 175–700 km (lines 89–93)

For each (time, lat, lon) column, a level $\ell$ is "in range" if
$175 \le z_\ell \le 700$ km and both $U_\ell$ and $V_\ell$ are finite. With
$N_{\mathrm{lev}}$ in-range levels:

$$
\bar U = \frac{1}{N_{\mathrm{lev}}}\sum_{\ell\ \mathrm{in\ range}} U_\ell ,\qquad
\bar V = \frac{1}{N_{\mathrm{lev}}}\sum_{\ell\ \mathrm{in\ range}} V_\ell ,
$$

and NaN where $N_{\mathrm{lev}} = 0$. The `np.where(in_range, ui, 0.0)` /
`nansum` construction implements the masked sum without Python loops; the
`errstate` guard silences the 0/0 warning.

*Rationale.* The code says the range is "the range the RISR-N vector
velocities are resolved over" (`config.py`). Physically (§4.0.1) the E×B drift
is nearly height-independent there, so the average mostly suppresses small
level-to-level differences.

*It is an average over levels, not over height.* **[interpretation]** WACCM-X's
upper levels are spaced roughly uniformly in log-pressure (a fraction of a
scale height), not in km, so each level represents a different thickness in
km, and the mean is weighted per scale height rather than per km. Because the
drift varies little with height this is a minor effect. The number of levels
in range varies between 30 and 37 (mostly 33–36) **[output]**, because the
thermosphere expands and contracts with solar activity and season, moving the
heights of pressure levels; the top of WACCM-X may lie below 700 km at
solar minimum, in which case the window is truncated at the model top.

#### 4.3.6 Bilinear interpolation to the bins — `_bilinear` (lines 52–64)

**Concept.** On a regular grid, a value at an arbitrary point is estimated
from the four surrounding grid values, weighted by the fractional distances
along each axis. With grid latitudes $\theta_i$, longitudes $\phi_j$ and target
$(\theta,\phi)$:

$$
w_x = \frac{\phi-\phi_j}{\phi_{j+1}-\phi_j},\qquad
w_y = \frac{\theta-\theta_i}{\theta_{i+1}-\theta_i},
$$

$$
f(\theta,\phi) = (1-w_y)(1-w_x)f_{i,j} + (1-w_y)w_x f_{i,j+1} + w_y(1-w_x)f_{i+1,j} + w_y w_x f_{i+1,j+1}.
$$

It is exact for fields linear in latitude and longitude separately, continuous,
and never overshoots the four corner values.

**Implementation details.**

* Longitude periodicity: the first longitude column is appended at
  $\phi_0+360^\circ$ (lines 56–57) and target longitudes are wrapped to
  [0, 360) (line 58), so a point between 357.5° and 360° interpolates between
  the last and first columns. Bin longitudes are ≈ −95° → 265°.
* Cell index: `searchsorted(..., side="right") - 1`, clipped to valid cells.
* Leading dimensions (`...`) are carried through, so one call interpolates all
  8 times at once; output shape `(time, bin)`.
* Applied to $\bar U$, $\bar V$ and, as a diagnostic, $N_{\mathrm{lev}}$
  (hence non-integer `nlev` values where neighbouring columns differ).

**Assumptions and caveats.**

* Interpolation is in degree space, not great-circle distance; fine away from
  the pole (bins are south of ≈ 81° N).
* East and north components are interpolated separately. Near the pole the
  local east/north basis rotates by the longitude difference between grid
  columns (2.5°), which is ignored **[interpretation]**; the resulting error is
  at most a few per cent of the vector magnitude.
* If any of the four corners is NaN the result is NaN (NumPy arithmetic).
* The model grid (cells ≈ 210 km north–south × 60 km east–west at 77° N) is
  far coarser than the radar bins (0.25° ≈ 28 km), so adjacent bins can share
  the same four grid points.

#### 4.3.7 Magnetic frame and cache row layout (lines 99–107)

The `(time, bin)` arrays are flattened into a long table: `time` repeated per
bin, `cgm_lat` tiled per time, `ve_geo` (from `UI`), `vn_geo` (from `VI`),
`nlev`. The magnetic components are computed with `geometry.geo_to_mag` using
**the same per-bin, per-year declination the radar used**, so radar and model
magnetic components are defined identically.

**Output** — `Outputs/waccmxsd_cache/<YYYY-MM-DD>.csv`, `%.6g` floats, 8 times
× 23 bins = 184 rows per normal day:

| Column | Meaning | Units |
|---|---|---|
| `time` | Model output time (instantaneous) | UTC |
| `cgm_lat` | Radar bin | deg |
| `ve_geo`, `vn_geo` | 175–700 km mean of `UI`, `VI` at the bin | m s⁻¹ |
| `nlev` | Interpolated number of model levels in the window | — |
| `ve_mag`, `vn_mag` | Rotated into the local magnetic frame | m s⁻¹ |

The reference run extracted 3,653 days (2008-01-01 to 2017-12-31) with no
failures, so `_skipped.csv` holds only its header **[output]**.

**Validation and error handling.** Missing variables raise `KeyError`,
collected per day. **Not detected:** days that read successfully but contain
fewer outputs or all-NaN fields. Two such days exist in the cache
**[output]**: `2008-01-11.csv` has only 3 output times, the third (06:00) with
`nlev = 0` and empty velocities; `2017-12-31.csv` has 7 times, the last (18:00)
likewise empty. In both, the final output is all-NaN, which suggests truncated
files **[interpretation]**. These gaps are later handled by the time-matching
rules (§4.4.2), not by stage 2. (This is also why the log reports the model
period as ending at "2017-12-31 18:00:00".)

`load_cache(cache)` (lines 140–144) reads every `????-??-??.csv` (so not
`_skipped.csv` or `.tmp` files), concatenates them and sorts by
`(cgm_lat, time)`. It raises `FileNotFoundError` if the cache is empty.

**Important coupling.** The cache files are keyed only by date. If stage 1 is
re-run and produces a different set of bins, existing cache days are *not*
re-extracted unless `--overwrite` is given; bins absent from the cache simply
get no model value in stage 3 (`match` skips them) **[code]**.

---

### 4.4 Stage 3 — `climatology`: matching, binning, scoring, plotting

Entry: `__main__.step_climatology` (lines 48–76). Order of operations:

1. `obs = _hourly(out)` — read `risrn_hourly.csv.gz`.
2. `model = waccmxsd.load_cache(out/"waccmxsd_cache")`.
3. Optional date filtering of both by `--start` / `--end` (end inclusive of the
   whole day: `< end + 1 day`).
4. `matched = climatology.match(obs, model)` → write `matched_hourly.csv.gz`;
   log how many rows have a model value for the first component.
5. `model_all = climatology.model_all_times(model)`.
6. `climatology.scores(matched, band)` → `stats_waccmxsd_<LO>-<HI>.csv`.
7. For each component in `--components` and each `lt` in `--lt`: `diurnal`
   table + plot, `lat_lt` table + plot.

#### 4.4.1 Inputs

`obs` (112,081 rows, schema §4.2.7) and `model` (≈ 672,000 rows, schema §4.3.7).

#### 4.4.2 `match(obs, model)` — linear interpolation in time (`climatology.py:27–47`)

**Purpose.** Give every radar hourly mean a model value at the same bin and the
same instant (hh:30), so they can be paired.

**First principles.** Given samples $(t_0, f_0)$ and $(t_1, f_1)$ with
$t_0 \le t < t_1$, linear interpolation is

$$
f(t) = f_0 + (f_1 - f_0)\,\frac{t - t_0}{t_1 - t_0}.
$$

With 3-hourly model output and hh:30 radar stamps, $t$ is always 0.5, 1.5 or
2.5 h after $t_0$, so the weights on $(f_0, f_1)$ are $(5/6, 1/6)$, $(1/2, 1/2)$ or
$(1/6, 5/6)$.

**Algorithm**, per bin (`groupby("cgm_lat")`):

1. Select the model rows for that bin (`np.isclose` on latitude) sorted by time;
   skip the bin if there are none.
2. Convert times to `datetime64[s]`; `k = searchsorted(model_times, obs_times)`
   gives, for each observation, the index of the first model time ≥ it. The
   bracketing pair is `lo = k−1`, `hi = k` (clipped).
3. A value is accepted (`inside`) only if the observation lies strictly within
   the model record ($0 < k < N$) **and** the bracketing outputs are at most
   `MODEL_MAX_GAP_H` = 3.5 h apart.
4. For each of the four components, `np.interp` (on times as float seconds)
   gives the interpolated value; it is kept only if `inside` and both
   bracketing values are finite; otherwise NaN.
5. Results are written to new columns `model_vn_geo`, `model_ve_geo`,
   `model_vn_mag`, `model_ve_mag`.

**Why 3.5 h?** The normal spacing is 3 h; 3.5 h tolerates small irregularities
but refuses to bridge a missing output (6 h) or a missing day (≥ 24 h), where a
straight line would be fiction **[interpretation]**. The finiteness check
handles the all-NaN outputs at 2008-01-11 06:00 and 2017-12-31 18:00.

**Validation.** `tests/test_core.py::test_match_interpolates_in_time_and_refuses_gaps`
checks the midpoint of 0 and 300 m s⁻¹ is 150 m s⁻¹ and that an observation in
a 24 h gap gets NaN.

**Output.** `matched` = `obs` + four `model_*` columns, written to
**`Outputs/matched_hourly.csv.gz`**. In the reference run 108,352 of 112,081
rows got model values; all 3,729 unmatched rows are from 2018, after the end of
the model archive **[output]**.

**Numerical side-effect — damping of the diurnal cycle.** **[interpretation]**
Linear interpolation of a sinusoid of period $P$ between samples $h$ apart
damps it. For the 24 h harmonic ($\omega = 2\pi/24$ h⁻¹) and the three hh:30
offsets, the amplitude factors are

$$
\left|\tfrac56 e^{-i\omega\cdot0.5} + \tfrac16 e^{i\omega\cdot2.5}\right| \approx 0.958,\qquad
\cos(\omega\cdot1.5) = \cos(\pi/8) \approx 0.924,\qquad 0.958,
$$

an average ≈ 0.95. The model's diurnal amplitude is therefore biased low by
≈ 5 % relative to the instantaneous model, and faster variations (e.g.
substorm-scale changes) are smoothed much more. This affects both `model` and
`model_all` equally.

#### 4.4.3 `model_all_times(model)` — the model at every hour (`climatology.py:57–73`)

**Purpose.** A model climatology over the entire model period, independent of
when the radar ran. The docstring explains why the model is interpolated to
every hour rather than used at its native 3-hourly times: at a fixed longitude,
outputs at 00, 03, …, 21 UT fall at only 8 distinct local hours, leaving two
thirds of the local-hour bins empty.

**Steps.**

1. Build an hourly grid from `floor(min time)` to `ceil(max time)`, shifted to
   hh:30 (matching the radar stamps), × every bin in the cache.
2. Reuse `match(grid, model)` — the same interpolation and gap rules.
3. Rename `model_*` columns to plain component names; drop rows where all four
   are NaN.
4. Attach `glon`, `mlon` via `positions_for`, then `slt`, `mlt` exactly as in
   stage 1.

For 2008–2017 this is ≈ 87,600 hours × 23 bins ≈ 2 million rows.

*Consistency check* **[output]**. In summer (no gaps) every SLT hour of
`diurnal_*_slt` has `model_all_n` = 92 days × 10 years × 5 band bins = 4,600
exactly. Winter has 903 DJF days, but the two short days of §4.3.7 remove
hours: on 2008-01-11 UT 03:30–23:30 (the 06:00 output is NaN and the next
valid output is 18 h later), and on 2017-12-31 UT 15:30–23:30. With
SLT ≈ UT − 6 h in the band, the counts are 903 × 5 = 4,515 for SLT 18–20
(UT 00–02, never lost), 902 × 5 = 4,510 for SLT 21–08 (lost once) and
901 × 5 = 4,505 for SLT 09–17 (lost on both days) — exactly the values in the
table.

#### 4.4.4 `scores(matched, band)` — skill statistics (`climatology.py:106–127`)

**Purpose.** Quantify how well the model reproduces the hourly observations in
the station band.

**Steps.** Select rows with $82.5 \le$ `cgm_lat` $\le 83.5$ (inclusive at both
ends); add `season` (by UT calendar month: DJF, MAM, JJA, SON); for each of the
four components and each season plus `ALL`, drop rows lacking either value;
skip groups with fewer than 3 pairs. With observations $o_i$, model values
$m_i$, $d_i = m_i - o_i$, $i=1\ldots N$:

| Column | Definition | Meaning |
|---|---|---|
| `N` | number of pairs | |
| `Bias (m/s)` | $\bar d = \frac1N\sum d_i$ | Mean model-minus-observation; positive = model too northward/eastward. |
| `RMSE (m/s)` | $\sqrt{\frac1N\sum d_i^2}$ | Typical error size, sensitive to large errors. |
| `MAE (m/s)` | $\frac1N\sum\lvert d_i\rvert$ | Typical error size, robust to outliers. |
| `r` | $\dfrac{\sum(o_i-\bar o)(m_i-\bar m)}{\sqrt{\sum(o_i-\bar o)^2\sum(m_i-\bar m)^2}}$ | Pearson correlation of hourly values. |
| `Obs std`, `Model std` | $\sigma = \sqrt{\frac1N\sum(x_i-\bar x)^2}$ | Variability; **population** std (NumPy default, $N$ denominator). |
| `Std ratio` | $\sigma_m/\sigma_o$ | < 1 means the model is too smooth. |

plus `Model, Season, Component, Band, Start, End` (first/last date in the
group). Values are rounded to 3 decimals.

These statistics are linked by the exact identity (with population stds)

$$
\mathrm{RMSE}^2 = \bar d^{\,2} + \sigma_o^2 + \sigma_m^2 - 2r\sigma_o\sigma_m ,
$$

which is a useful consistency check on the table (§5.3). Note: the `_std`
columns elsewhere (hourly means, diurnal tables) use the sample ($N-1$)
definition; `scores` uses $N$. For thousands of pairs the difference is
negligible.

The docstring says the table follows "the style of the FPI wind_stats tables";
those tables are not in this repository.

#### 4.4.5 `diurnal(...)` — season × local-hour climatology (`climatology.py:82–95`)

**Purpose.** Show the mean daily cycle over Resolute in each season for the
radar, the model at radar times, and the model at all times.

**Binning** — `add_time_bins(df, lt_col)` (lines 50–54):
`season` from the UT month; `lt_hour = floor(lt) mod 24`, so hour $h$ covers
$[h, h+1)$ local time.

**Steps.**

1. Select rows in the band.
2. `mt`: matched rows with **both** observation and model values for the
   component. Consequently the "obs" climatology also excludes 2018, which has
   no model.
3. `ma`: `model_all` rows in the band with a value.
4. Group each by `(season, lt_hour)` and summarise with `_summary`
   (lines 76–79): mean, sample std, median, 25th and 75th percentiles, count.
5. Concatenate as `obs_*`, `model_*`, `model_all_*` columns.

**Output** — `Outputs/diurnal_<component>_<slt|mlt>_<LO>-<HI>.csv` (`%.4g`):
`season, lt_hour`, then for each of `obs`, `model`, `model_all`:
`mean, std, median, q25, q75, n` (m s⁻¹ except `n`).

**How to read `obs` vs `model` vs `model_all`.**

* `obs` vs `model`: like-for-like, same hours and bins — differences are model
  error.
* `model` vs `model_all`: the same model, different sampling — the README says
  this "shows the effect of the radar's campaign sampling". **[interpretation]**
  Strictly it mixes three sampling effects: (i) which days the radar ran;
  (ii) the period — `model_all` covers 2008-01 to 2017-12, the radar
  2009-08 to 2017-11, so `model_all` includes ≈ 1.6 years the radar never saw,
  near solar minimum; (iii) **latitude weighting** — `model_all` weights the 5
  band bins equally, whereas the matched pairs are dominated by 83.25–83.5°
  (in the reference run: 6, 1,260, 4,399, 5,579, 7,862 pairs at 82.5, 82.75,
  83.0, 83.25, 83.5°). The model varies smoothly with latitude in the band, so
  (iii) is probably small.

#### 4.4.6 `lat_lt(...)` — season × latitude × local-hour (`climatology.py:98–103`)

Using all matched rows (all 23 bins) with both values, group by
`(season, cgm_lat, lt_hour)` and record the observed mean, model mean and
count. **Output** — `Outputs/lat_lt_<component>_<slt|mlt>.csv` (`%.4g`):
`season, cgm_lat, lt_hour, obs, model, n`. The difference (model − obs) is
computed only in the plot, not stored, and the CSV includes cells with any
`n ≥ 1`.

#### 4.4.7 Figures (`plots.py`)

Matplotlib with the `Agg` backend (no display needed on a cluster), serif font,
200 dpi PNGs in `Outputs/plots/`.

**`plots.diurnal`** (lines 30–73) — one row per season present:

* x positions: `_roll(lt_hour, start) = ((h + 0.5 − start) mod 24) + start`
  places each hour at its bin centre on a `start … start+24` axis
  (noon-to-noon for SLT, 0–24 for MLT); ticks every 3 h labelled modulo 24.
* Dotted red: `model_all_mean`. Solid red with markers: `model_mean`, red
  shading ±1 `model_std`. Black dots: `obs_mean`, grey shading ±1 `obs_std`.
  Dashed zero line.
* Title per panel: season and `N` = Σ `obs_n`. Shared symmetric y-limits =
  1.05 × the largest |mean ± std| across panels.
* The shading is the **spread of hourly means** (day-to-day variability), not
  the uncertainty of the climatological mean (which would be ≈ std/√n).

**`plots.lat_lt`** (lines 76–110) — rows = seasons, columns = RISR-N,
WACCM-X SD, WACCM-X SD minus RISR-N:

* Cells with `n < 3` (`min_n`) are hidden.
* `pcolormesh` with x edges 0…24 (h) and y edges at bin centres ± half the
  median spacing.
* Diverging `RdBu_r` colour map, symmetric limits: 98th percentile of
  |obs| and |model| for the first two columns, of |model − obs| for the third
  (robust to a few extreme cells).

**Outputs of stage 3**: see §5.

---

## 5. Final outputs and their interpretation

### 5.1 Inventory (reference run)

| File | Produced by | Rows / content |
|---|---|---|
| `risrn_files.csv` | stage 1 | 615 files (290 used) |
| `risrn_hourly.csv.gz` | stage 1 | 112,081 hourly bin means |
| `waccmxsd_cache/<date>.csv` ×3653, `_skipped.csv` | stage 2 | 184 rows per normal day; no failures |
| `matched_hourly.csv.gz` | stage 3 | 112,081 rows; 108,352 with model |
| `stats_waccmxsd_82.5-83.5.csv` | stage 3 | 4 components × 5 seasons |
| `diurnal_{vn_geo,vn_mag}_{slt,mlt}_82.5-83.5.csv/.png` | stage 3 | 4 seasons × 24 hours |
| `lat_lt_{vn_geo,vn_mag}_{slt,mlt}.csv/.png` | stage 3 | season × bin × hour |

### 5.2 Diurnal tables and figures

**What they represent.** For one velocity component, the distribution, over
all matched hourly means in the band, of the drift as a function of season and
local hour. Units m s⁻¹; positive = northward (or eastward). `lt_hour` = $h$
means local time in $[h, h+1)$; figures plot it at $h+0.5$.

**Trace back.** `obs_*` ← `vn_geo` hourly means (§4.2.6) ← QC'd records
(§4.2.4) ← radar vectors. `model_*` ← `match` (§4.4.2) ← cache (§4.3) ←
`UI`/`VI`. `model_all_*` ← `model_all_times` (§4.4.3). `lt_hour` ← `slt`
(bin longitude) or `mlt` (radar AACGM longitude) (§4.2.7).

**What the reference run shows** (`diurnal_vn_geo_slt_82.5-83.5.png`):
in every season the northward geographic drift follows the antisunward
pattern of §4.0.2 — ≈ +250 to +350 m s⁻¹ in the late morning, ≈ −300 to
−400 m s⁻¹ near local midnight, crossing zero near 04:30–05:30 and
17:00–18:30 SLT. The
model reproduces the phase and amplitude well in the night sector; around
noon it is more northward than observed (DJF 12 h: obs 171.8, model
299.2 m s⁻¹ **[output]**). The model's spread (red band) is much narrower than
the observed spread (grey band). The dotted all-times curve lies close to the
at-radar-times curve, i.e. campaign sampling changes the model's mean little
in this band.

### 5.3 Statistics table

`stats_waccmxsd_82.5-83.5.csv` (reference run, `ALL` rows):

| Component | N | Bias | RMSE | MAE | r | Obs std | Model std | Std ratio |
|---|---|---|---|---|---|---|---|---|
| `vn_geo` | 19106 | 25.9 | 272.2 | 201.8 | 0.654 | 353.1 | 276.4 | 0.78 |
| `ve_geo` | 19106 | 21.7 | 278.3 | 204.0 | 0.570 | 337.4 | 180.1 | 0.53 |
| `vn_mag` | 19106 | 16.3 | 272.8 | 202.1 | 0.646 | 353.8 | 262.3 | 0.74 |
| `ve_mag` | 19106 | 30.0 | 277.7 | 204.8 | 0.572 | 336.6 | 200.1 | 0.59 |

**Interpretation.** Biases are small (≈ 15–30 m s⁻¹, < 10 % of the
variability); the model captures roughly 33–43 % of the hourly variance
($r^2$); and it is too smooth (std ratio 0.5–0.8), especially in the east
component. **[interpretation]** Part of the RMSE is irreducible by
construction: radar hourly means include real sub-3-hour variability and
measurement noise, while the model value is an interpolated 3-hourly snapshot
from a model whose high-latitude driving is itself an average pattern
(§4.0.4).

**Plausibility checks you can do.**

* The identity of §4.4.4 for `vn_geo` ALL:
  $25.937^2 + 353.124^2 + 276.443^2 - 2(0.654)(353.124)(276.443)$
  $= 672.7 + 124{,}696.8 + 76{,}420.5 - 127{,}676 \approx 74{,}114$, and
  $\sqrt{74{,}114} = 272.2$ = the tabulated RMSE ✔.
* Recomputing the `vn_geo` ALL row directly from `matched_hourly.csv.gz`
  (band rows with both values) reproduces N = 19,106, bias 25.937, RMSE
  272.242, MAE 201.785, r 0.654, stds 353.124 / 276.443 exactly
  **[output]**.
* Seasonal N sum to ALL: 4770 + 6032 + 3070 + 5234 = 19,106 ✔.
* `vn_geo` and `vn_mag` have nearly equal RMSE because they differ only by a
  ≈ 20° rotation.

### 5.4 Latitude × local-time tables and figures

`lat_lt_<c>_<lt>.csv` gives the observed and model means in each
(season, bin, hour) cell; the figure adds model − obs. They show how the
agreement varies across the radar field of view (82.5–88°). In the reference
run (`lat_lt_vn_geo_mlt.png`) the model pattern is almost independent of
latitude (vertical stripes in MLT), the observations show more latitudinal
structure, and the difference panels show a model excess of northward flow on
the morning/noon side and a deficit on the afternoon/evening side
**[output, visual reading]**. Cells at the highest and lowest bins rest on few
hourly means (`n` column) and should be read with care.

### 5.5 Matched table

`matched_hourly.csv.gz` is the most useful file for further analysis: every
column of `risrn_hourly` plus the four `model_*` columns. A row with empty
`model_*` had no model value (outside the model period, or across a gap).

### 5.6 Which upstream assumptions affect which output

| Assumption / choice | Hourly obs | Model cache | Matched/scores | Diurnal | Lat–LT |
|---|---|---|---|---|---|
| QC thresholds (150, 3000 m s⁻¹) | ✔ | | ✔ | ✔ | ✔ |
| Cadence preference | ✔ (n, noise) | | ✔ | ✔ | ✔ |
| Bins on radar AACGM meridian at 300 km, mid-year epoch | ✔ (position, decl, SLT) | ✔ (sampling point) | ✔ | ✔ | ✔ |
| IGRF declination | ✔ (frame rotation) | ✔ | ✔ (mag comps) | ✔ | ✔ |
| 175–700 km level mean; Z3 conversion | | ✔ | ✔ | ✔ | ✔ |
| Bilinear interpolation | | ✔ | ✔ | ✔ | ✔ |
| Linear time interpolation, 3.5 h gap rule | | | ✔ | ✔ (both model series) | ✔ |
| hh:30 stamping | ✔ | | ✔ | ✔ | ✔ |
| Band 82.5–83.5 | | | ✔ (scores) | ✔ | |
| Season by UT month, lt floor binning | | | ✔ (season) | ✔ | ✔ |

---

## 6. Numerical methods, assumptions and limitations

### 6.1 Numerical methods summary

| Method | Where | Order / properties | Main error source |
|---|---|---|---|
| Rotation by declination | `geometry.mag_to_geo/geo_to_mag` | Exact, norm-preserving | Declination model error; mid-year epoch |
| AACGM ↔ geographic | `aacgmv2` | Fitted-coefficient approximation of field-line tracing | Library coefficient accuracy (small) |
| Geopotential → geometric height | `waccmxsd.extract_day` | Exact for spherical Earth, $g\propto r^{-2}$ | Model `Z3` definition [uncertain] |
| Masked arithmetic mean over levels | `extract_day` | Per-level weighting | Not thickness-weighted |
| Bilinear interpolation | `waccmxsd._bilinear` | Exact for bilinear fields; second-order in grid spacing | Coarse model grid; vector-basis rotation ignored |
| Linear time interpolation | `climatology.match` | Second-order in Δt | Damps the diurnal cycle by ≈ 5 % (§4.4.2) |
| Group means / quantiles | pandas | Exact | Sampling noise in sparse cells |

### 6.2 Physical assumptions

1. Radar vector velocities are horizontal (perpendicular ≈ horizontal at
   $I\approx87$–$88^\circ$).
2. The radar product represents the F-region drift over ≈ 175–700 km, and the
   model level mean over the same range is comparable.
3. Each magnetic-latitude bin can be represented by a single point on the
   radar's AACGM meridian at 300 km.
4. Magnetic coordinates and declination are fixed for each calendar year
   (evaluated on 1 July).
5. A 3-hourly instantaneous model state, linearly interpolated, represents the
   model at the radar's hour centre.
6. Solar local time is mean solar time; seasons are meteorological, by UT
   month, and are the same for all years.

### 6.3 Limitations and potential failure modes

* **Representativeness mismatch.** Radar = hourly average of ~3-min vectors;
  model = interpolated 3-hourly snapshots on a ~2° grid. This depresses
  correlation and the model/obs std ratio irrespective of model physics.
* **Silent partial days.** Stage 2 does not flag days with missing or all-NaN
  outputs (e.g. `2008-01-11`); they are only handled implicitly by the gap
  rule in stage 3.
* **Cache/bin coupling.** The cache does not record which bins or altitude
  range it was built with; changing `MODEL_ALT_RANGE_KM` or the bin set
  requires `--overwrite`, or the cache will silently mix configurations.
* **`_skipped.csv` scope.** It is rewritten each run, so after a run with
  `--start/--end` it lists only failures in that range.
* **Overlapping experiments.** Records from different experiments falling in
  the same (hour, bin) would be averaged together; `experiment` would then
  record only the first. Whether this occurs was not checked **[uncertain]**.
* **Stage 1 ignores `--start/--end`** and always processes the full archive.
* **Uneven sampling within the band.** Scores and diurnal "obs/model" are
  dominated by bins 83.25–83.5°; the 82.5° bin contributes 6 pairs.
* **2018 radar data are unused** in all comparisons because the model archive
  ends on 2017-12-31.
* **Small cells.** Diurnal hours with tens of hourly means, and lat–LT cells
  with 3–10, have large sampling error; the plots show spread (±1 std), not
  the standard error of the mean.
* **CSV precision.** Intermediate files are written with 6 significant digits
  (`%.6g`), and tables with 4 (`%.4g`); harmless for m s⁻¹-level work but
  not bit-reproducible with an in-memory pipeline.

### 6.4 Sensitivities worth knowing

* **QC thresholds** remove ≈ 40 % of finite vectors; looser thresholds would
  add noisier, likely larger-magnitude vectors, increasing obs std and RMSE.
* **Band choice** controls which bins enter scores; the band starts at 82.5
  because 82.25 is never populated (config comment).
* **Gap limit** (3.5 h): raising it to > 24 h would interpolate across missing
  days and fabricate model values.

---

## 7. Validation, testing and reproducibility

### 7.1 Unit tests (`tests/test_core.py`)

Run with `python -m pytest tests`. The tests need no data archive but do need
`aacgmv2` and `ppigrf` (imported by `geometry`).

| Test | What it establishes |
|---|---|
| `test_frame_rotation_round_trip` | `mag_to_geo(geo_to_mag(e, n, D), D)` returns the input for 50 random vectors and declinations in [−60°, 60°] — the two rotations are inverses. |
| `test_magnetic_north_points_along_declination` | A unit magnetic-north vector maps to geographic azimuth $D$ ($D=-20^\circ$): fixes the **sign convention** of the rotation. |
| `test_bilinear_reproduces_linear_fields_and_wraps_longitude` | `_bilinear` reproduces fields linear in latitude exactly, including at longitudes 359° and −1° (wrap-around), and a field equal to longitude at an interior point. |
| `test_match_interpolates_in_time_and_refuses_gaps` | Midpoint interpolation (0 → 300 gives 150) and NaN across a 24 h gap. |
| `test_bin_positions_place_lowest_bin_over_resolute` | For 2012, bin 82.5° lies within 1° latitude / 2° longitude of Resolute (74.7° N, 95° W), and geographic latitude increases with magnetic latitude. |

**What is not tested:** file scanning/selection, QC, hourly aggregation, the
height conversion and vertical averaging, unit scaling, `model_all_times`,
`diurnal`, `lat_lt`, `scores`, plotting, and the CLI. These were instead
cross-checked for this document against the reference outputs (§5.3, §9).

**Test status.** The tests were **not run** while writing this document: the
machine used had no Python scientific stack. Their results are therefore not
reported here.

### 7.2 Independent checks performed for this document

Using only the Python standard library on the files in `Outputs/`:

* stats `vn_geo` ALL row reproduced exactly from `matched_hourly.csv.gz`;
* diurnal `DJF` SLT hours 0 and 12 (`obs_mean`, `obs_std`, `obs_n`,
  `model_mean`) reproduced exactly;
* `lat_lt` cell (DJF, 82.75°, 0 h) reproduced exactly;
* time interpolation, frame rotation and SLT for the worked example (§9)
  reproduced to the printed precision;
* counts in the log (615 files, 290 experiments, 1,856,573 records, 112,081
  bins, 108,352 matched) agree with the CSVs.

### 7.3 Running the pipeline

```bash
pip install -r requirements.txt            # in the target environment
python -m isr_clim risrn                   # → Outputs/risrn_hourly.csv.gz
python -m isr_clim extract --workers 8     # → Outputs/waccmxsd_cache/
python -m isr_clim climatology             # → tables, scores, figures
# or
python -m isr_clim all --workers 8
# on BlueBEAR
sbatch slurm/submit_climatology.sh all
```

Input requirements: the RISR-N tree under `--risrn-dir` and the WACCM-X SD
`cam.h1` files under `--sd-dir`, with the file name pattern of §4.3.2.

### 7.4 Reproducibility considerations

* **Deterministic given inputs and library versions.** There is no randomness
  outside the tests (which use a fixed seed). Parallel extraction writes one
  file per day, so worker order does not matter.
* **Unpinned dependencies.** `aacgmv2` coefficient sets and the IGRF generation
  shipped with `ppigrf` change between releases; both alter bin positions,
  declinations, MLT and therefore every downstream number slightly. Recording
  the versions (e.g. `pip freeze`) alongside outputs is advisable.
* **Cache persistence.** Stage 2 is incremental; outputs depend on what is
  already in the cache (§6.3).
* **Configuration lives in code.** Thresholds are constants in `config.py`;
  changing them is a code change, not a run-time option (except band,
  components, local-time coordinate and dates).

---

## 8. Glossary and symbol reference

### 8.1 Terms

| Term | Meaning |
|---|---|
| AACGM | Altitude-Adjusted Corrected Geomagnetic coordinates (§4.2.4). `cgm_lat` in the radar files is an AACGM-type latitude. |
| Bin | One magnetic-latitude cell of the radar vector-velocity product; 0.25° wide, 82.5–88.0° here. |
| Cadence | Integration time of a radar file (`inttms`), s. |
| Declination $D$ | Angle of the horizontal geomagnetic field east of geographic north. |
| DJF/MAM/JJA/SON | Dec–Feb, Mar–May, Jun–Aug, Sep–Nov. |
| E×B drift | Plasma drift perpendicular to both E and B (§4.0.1). |
| edynamo | WACCM-X electrodynamics module producing `UI`/`VI`. |
| Experiment | A Madrigal radar run (folder under `ran/`). |
| F region | Ionosphere ≈ 150–600 km, where ions E×B drift. |
| Geo / mag frame | Local east/north components relative to geographic / magnetic north. |
| Geopotential height | $Z=\Phi/g_0$ (§4.3.4). |
| IGRF | International Geomagnetic Reference Field. |
| Inclination $I$ | Angle of B below horizontal. |
| ISR | Incoherent Scatter Radar. |
| Kindat 5960 | Madrigal data-kind code of the RISR-N long-pulse vector-velocity product. |
| Madrigal | Database distributing ISR data. |
| MLT | Magnetic local time (§4.2.7). |
| RISR-N | Resolute Bay Incoherent Scatter Radar – North. |
| SD | Specified Dynamics (nudged to reanalysis). |
| SLT | Solar (mean) local time. |
| Station band | Magnetic-latitude range over Resolute used for scores/diurnal tables. |
| WACCM-X | Whole Atmosphere Community Climate Model with thermosphere–ionosphere eXtension. |

### 8.2 Symbols

| Symbol | Meaning | Units | Code |
|---|---|---|---|
| $\mathbf E$, $\mathbf B$ | Electric, magnetic field | V m⁻¹, T | — |
| $\mathbf v_\perp$ | E×B drift | m s⁻¹ | — |
| $\Omega_i$, $\nu_{in}$ | Ion gyrofrequency, ion–neutral collision frequency | s⁻¹ | — |
| $v_{e,g}, v_{n,g}$ | Geographic east, north components | m s⁻¹ | `ve_geo`, `vn_geo` |
| $v_{e,m}, v_{n,m}$ | Magnetic east, north components | m s⁻¹ | `ve_mag`, `vn_mag` |
| $\delta v_e, \delta v_n$ | Radar error estimates | m s⁻¹ | `dve`, `dvn` |
| $D$ | Declination | deg | `decl` |
| $I$ | Inclination | deg | — |
| $B_e, B_n, B_u$ | IGRF field components | nT | `be`, `bn` |
| $\lambda$ | AACGM latitude | deg | `cgm_lat` |
| $\Lambda$ | AACGM longitude | deg | `mlon` |
| $\theta, \phi$ | Geographic latitude, longitude | deg | `glat`, `glon` |
| $Z$, $z$ | Geopotential, geometric height | km | `z3`, `z` |
| $R$ | Earth radius, 6371 | km | `EARTH_RADIUS_KM` |
| $U, V$ | Model zonal, meridional ion drift | m s⁻¹ | `UI`, `VI` |
| $N_{\mathrm{lev}}$ | Model levels in 175–700 km | — | `nlev` |
| $w_x, w_y$ | Bilinear weights | — | `wx`, `wy` |
| $n$ | Records per hourly mean | — | `n` |
| $o_i, m_i, d_i$ | Observed, modelled value, difference | m s⁻¹ | `comp`, `model_comp`, `diff` |
| $N$ | Number of matched pairs | — | `N` |
| $r$ | Pearson correlation | — | `r` |
| $\sigma_o, \sigma_m$ | Population stds | m s⁻¹ | `Obs std`, `Model std` |

---

## 9. End-to-end worked example

This traces **one real hourly mean** through the pipeline using values from the
reference run's outputs: bin `cgm_lat = 83.5`, hour 2009-08-24 16:00–17:00 UT
(the first row of `risrn_hourly.csv.gz`). All numbers are **[output]** unless
marked hypothetical.

**Step 1 — file selection.** Experiment `24aug09b` has two files
(`risrn_files.csv`):

| File | Cadence | n_records | n_finite | Used |
|---|---|---|---|---|
| `ran090824.003.hdf5` | 180 s | 960 | 715 | **True** (n_good 424) |
| `ran090824.004.hdf5` | 500 s | 336 | 263 | False |

180 s has rank 0 and 500 s rank 4, so `.003` is chosen. The experiment spans
16:14:31–18:22:24 UT. 960 records over ≈ 43 integration periods means ≈ 22
bins per period. QC keeps 424 of 715 finite vectors (59 %).

**Step 2 — records in the hour.** Records for bin 83.5 whose mid-times fall in
16:00–16:59 are grouped. The experiment starts at 16:14:31, leaving ≈ 15
three-minute periods in that hour; 14 pass QC (`n = 14`). The file is in the
"geo" layout, so `vn_geo`/`ve_geo` are the file's values.

**Step 3 — hourly mean and stamp.** Mean values: `vn_geo` = 205.071,
`ve_geo` = −45.5714 m s⁻¹ (std 121.2 and 86.4). Mean reported errors:
`dvn` = 74.6, `dve` = 33.6 m s⁻¹. Stamped `time = 16:30:00`.

**Step 4 — bin position (2009 table).** `mlon` = −36.0487° (radar's AACGM
longitude, 2009-07-01); the bin is placed at `glat` = 75.4781° N,
`glon` = −94.3759°; IGRF declination `decl` = −23.3777°.

**Step 5 — magnetic components** (each record was rotated; check on the means,
valid because $D$ is fixed): $\cos D = 0.91791$, $\sin D = -0.39678$.

$$
v_{e,m} = (-45.5714)(0.91791) - (205.071)(-0.39678) = -41.83 + 81.37 = 39.54
$$
$$
v_{n,m} = (205.071)(0.91791) + (-45.5714)(-0.39678) = 188.24 + 18.08 = 206.32
$$

File: `ve_mag` = 39.5399, `vn_mag` = 206.319 ✔.

**Step 6 — local times.**
$\mathrm{SLT} = 16.5 + (-94.3759)/15 = 16.5 - 6.2917 = 10.2083$ h ✔ (`slt`).
`mlt` = 9.028 h from `aacgmv2` (MLT ≈ SLT − 1.2 h here).

**Step 7 — model extraction** (`waccmxsd_cache/2009-08-24.csv`, bin 83.5):

| Time | `ve_geo` | `vn_geo` | `nlev` | `ve_mag` | `vn_mag` |
|---|---|---|---|---|---|
| 15:00 | −47.4299 | 203.010 | 33 | 37.0159 | 205.164 |
| 18:00 | 31.2116 | 195.008 | 33 | 106.027 | 166.615 |

Each is a mean over 33 model levels with geometric height in 175–700 km,
bilinearly interpolated to (75.48° N, 265.62° E), and rotated with the same
$D = -23.3777°$ (check at 15:00:
$-47.4299\times0.91791 - 203.010\times(-0.39678) = 37.02$ ✔).

**Step 8 — time matching.** 16:30 lies 1.5 h after 15:00 in a 3 h interval
(≤ 3.5 h, both ends finite), so weights are ½, ½:

| Component | Model | Observed | Model − obs |
|---|---|---|---|
| `vn_geo` | (203.010 + 195.008)/2 = **199.009** | 205.071 | −6.06 |
| `ve_geo` | (−47.4299 + 31.2116)/2 = **−8.109** | −45.571 | +37.46 |
| `vn_mag` | (205.164 + 166.615)/2 = **185.889** | 206.319 | −20.43 |
| `ve_mag` | (37.0159 + 106.027)/2 = **71.521** | 39.540 | +31.98 |

These are exactly the `model_*` values in the first row of
`matched_hourly.csv.gz` ✔.

**Step 9 — climatology binning.** August → `JJA`; `lt_hour` = ⌊10.2083⌋ = 10
(SLT) or ⌊9.028⌋ = 9 (MLT); 83.5 is inside [82.5, 83.5].

* It is one of the 115 hourly means behind `diurnal_vn_geo_slt_82.5-83.5.csv`
  row `JJA, 10` (obs mean 304.9, model 362.0, all-times model 363.0 m s⁻¹);
  this example (obs 205, model 199) lies below both means and well within one
  std (205.1).
* It is one of 43 hourly means in `lat_lt_vn_geo_slt.csv` cell
  `JJA, 83.5, 10` (obs 309.6, model 362.2).
* It is one of the 3,070 JJA pairs in `stats_waccmxsd_82.5-83.5.csv`, adding
  $d=-6.06$ m s⁻¹ to the bias sum and $36.7$ (m s⁻¹)² to the squared-error
  sum for `vn_geo`.

The flow shows the whole chain: one file of a two-hour campaign → 14 QC'd
3-minute vectors → one hourly mean at a computed geographic point → two model
snapshots at that same point → one interpolated model value → one pair in
three different summaries.

---

## 10. Appendix — discrepancies between documentation and code, and open questions

### 10.1 Comments or documentation that differ from the code

| Location | Says | Code does |
|---|---|---|
| `climatology.py` module docstring | "Times more than MODEL_MAX_GAP_H from the bracketing model outputs … get no model value" | The test is that the **bracketing outputs are ≤ 3.5 h apart** (`(mt[hi] - mt[lo]) <= max_gap`), not the observation's distance from them. `README.md` ("no value across gaps longer than 3.5 h") matches the code. |
| `risrn.py` `LAYOUTS` comment | "(east, north, east error, north error, frame)" | Tuples have four entries; the frame is the dictionary key. |
| `risrn.py` docstring | Mentions direction fields `gnangle` / `nangle` | Not read or used. |
| `risrn.py` docstring | Cadences "60/180/300/500 s" | `PREFERRED_CADENCE_S` also includes 120 and 30 s (and both occur in the archive). |
| `risrn.py` docstring | Perpendicular vs horizontal differ "by under 0.1%" at ~87° | $1-\sin 87^\circ = 0.14\,\%$; under 0.1 % needs $I > 87.4^\circ$. Negligible either way. |
| `README.md` | "the model at all times (shows the effect of the radar's campaign sampling)" | Also differs in period (2008–mid-2009 included) and latitude weighting (§4.4.5). |
| `README.md` | lat-LT products include "the model-minus-radar difference" | The difference is in the figure only, not in the CSV. |
| `waccmxsd.py` docstring | "8 three-hourly instantaneous outputs" per file | True for 3651 of 3653 days; `2008-01-11` (3 outputs) and `2017-12-31` (7 outputs), each ending in an all-NaN output, are exceptions that are not flagged. |

### 10.2 Open questions (cannot be settled from the repository)

1. How the WACCM-X SD run drives high-latitude convection (empirical model and
   its inputs), which bears directly on the std ratios.
2. The exact definition of `Z3` in this WACCM-X version (reference gravity,
   above sea level vs surface) and the model top height relative to 700 km.
3. The exact definition of `cgm_lat` and the inversion assumptions in the
   Madrigal kindat 5960 product (provider-side).
4. Whether different experiments ever overlap in time.
5. Why the 2008-01-11 and 2017-12-31 model files are shorter than 8 outputs.
6. The versions of `aacgmv2` and `ppigrf` used in the reference run.
