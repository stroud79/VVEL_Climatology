"""Default paths and analysis settings.

The project folder sits beside the data archives and FPI_ForwardModel on
the group RDS share (/rds/projects/r/reidb-waccm-x), so those paths are
taken relative to it. Each can be overridden with an environment variable
or on the command line.

FPI_ForwardModel is the reference for this comparison: the station's
coordinates, display name and time zone, the WACCM-X SD archive path, the
geopotential-to-geometric height conversion, the model-observation pairing
rule and the diurnal binning all come from it (isr_clim/fpi.py), never
from a copy here.
"""

from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SHARE_ROOT = PROJECT_ROOT.parent

FPI_ROOT = Path(os.environ.get("ISR_FPI_ROOT", SHARE_ROOT / "FPI_ForwardModel"))
RISRN_DIR = Path(os.environ.get("ISR_RISRN_DIR", SHARE_ROOT / "RISRN_data"))
OUTPUT_DIR = Path(os.environ.get("ISR_OUTPUT_DIR", PROJECT_ROOT / "Outputs"))

# The FPI station whose findings this validates (its name in FPI_ForwardModel's
# config/sites.yaml) and the FPI forward-model run they come from.
SITE = "Resolute"
FPI_COMBO = "waccmxsd_only"
# The FM_V2 model type behind that run; selects its label and colour in
# Plotting/climatology_plots.MODEL_DISPLAY.
MODEL_KEY = "waccmx_sd"

# RISR-N location (Madrigal experiment metadata). The radar is not an FPI
# site, so it is not in the FPI registry.
RADAR_LAT = 74.72955
RADAR_LON = 265.09424
# Altitude at which bin positions and magnetic coordinates are evaluated.
REF_ALT_KM = 300.0

# RISR-N vector-velocity selection and quality control.
PREFERRED_CADENCE_S = (180, 300, 120, 60, 500, 30)  # first available is used per experiment
MAX_VELOCITY_ERROR = 150.0   # m/s, applied to both horizontal components
# Degenerate fits in the archive report an error of exactly 1 m/s (e.g.
# -2537 m/s +/- 1, or an east component of exactly 0 +/- 1); the smallest
# genuine errors are ~6 m/s. Under inverse-variance weighting a handful of
# them would outweigh whole hours of data, so errors must exceed this.
MIN_VELOCITY_ERROR = 1.0     # m/s
MAX_SPEED = 3000.0           # m/s

# Model sampling: UI/VI averaged over geometric height through the layer the
# RISR-N vectors are fitted over (altb..alte of every kindat 5960 file), on
# a uniform height grid of MODEL_DZ_KM. The SD model top lies at ~450-600 km
# (solar minimum to maximum), so in practice the layer is truncated there
# (z_top_km in the cache).
MODEL_ALT_RANGE_KM = (175.0, 700.0)
MODEL_DZ_KM = 5.0

# Magnetic-latitude band over Resolute Bay (AACGM ~82.2-82.8 deg at 300 km,
# 2009-2018). The files' lowest bin, 82.25 deg, never holds a vector, so the
# band starts at 82.5.
STATION_BAND = (82.5, 83.5)

# Components compared, named as in the FPI plots: (east or north value,
# its error), both geographic. Positive north / positive east, as the FPI
# horizontal winds at azimuths 0 and 90.
COMPONENTS = {"Meridional": ("vn_geo", "dvn_geo"), "Zonal": ("ve_geo", "dve_geo")}

SEASONS = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
           6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}
SEASON_ORDER = ("DJF", "MAM", "JJA", "SON")
