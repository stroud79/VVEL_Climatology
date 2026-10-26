"""Geographic positions of the RISR-N magnetic-latitude bins, local magnetic
declination, and the rotation of magnetic-frame vectors to geographic.

The vector-velocity product resolves the flow in bins of corrected
geomagnetic latitude. Each bin is placed on the radar's magnetic meridian
(AACGM longitude of the radar at REF_ALT_KM) and converted to geographic
coordinates for the year in question, which is where the model is sampled.
"""

from __future__ import annotations

import datetime as dt
from functools import lru_cache

import aacgmv2
import numpy as np
import pandas as pd
import ppigrf

from . import config


@lru_cache(maxsize=None)
def bin_positions(year: int, cgm_lats: tuple) -> pd.DataFrame:
    """Geographic latitude/longitude (deg, lon in -180..180) and IGRF
    declination (deg, east positive) of each bin, for mid-``year``."""
    when = dt.datetime(year, 7, 1)
    _, radar_mlon, _ = aacgmv2.convert_latlon(
        config.RADAR_LAT, config.RADAR_LON, config.REF_ALT_KM, when, method_code="G2A")
    glat, glon, _ = aacgmv2.convert_latlon_arr(
        np.asarray(cgm_lats, float), np.full(len(cgm_lats), radar_mlon),
        config.REF_ALT_KM, when, method_code="A2G")
    glon = (np.asarray(glon) + 180.0) % 360.0 - 180.0
    be, bn, _ = ppigrf.igrf(glon, np.asarray(glat), config.REF_ALT_KM, when)
    decl = np.degrees(np.arctan2(np.ravel(be), np.ravel(bn)))
    return pd.DataFrame({"cgm_lat": np.asarray(cgm_lats, float), "glat": np.asarray(glat),
                         "glon": glon, "mlon": radar_mlon, "decl": decl})


def positions_for(times: pd.Series, cgm_lat: pd.Series) -> pd.DataFrame:
    """Bin positions matched row by row to (time, cgm_lat) pairs."""
    years = pd.DatetimeIndex(times).year
    out = []
    lats = tuple(sorted(set(np.round(cgm_lat.astype(float), 4))))
    for y in np.unique(years):
        pos = bin_positions(int(y), lats).set_index("cgm_lat")
        sel = years == y
        out.append(pos.reindex(np.round(cgm_lat[sel].astype(float), 4)).set_index(cgm_lat.index[sel]))
    return pd.concat(out).reindex(cgm_lat.index)


def mag_to_geo(e_mag, n_mag, decl_deg):
    """Rotate (east, north) components from the local magnetic frame to the
    geographic frame. Magnetic north lies at azimuth ``decl`` (east of
    geographic north)."""
    d = np.radians(decl_deg)
    return e_mag * np.cos(d) + n_mag * np.sin(d), n_mag * np.cos(d) - e_mag * np.sin(d)


def mag_to_geo_errors(de_mag, dn_mag, decl_deg):
    """1-sigma errors of :func:`mag_to_geo`'s components, propagated from
    independent magnetic-frame errors. The files give no east-north
    covariance, so it is taken as zero."""
    c2, s2 = np.cos(np.radians(decl_deg)) ** 2, np.sin(np.radians(decl_deg)) ** 2
    return np.sqrt(c2 * de_mag ** 2 + s2 * dn_mag ** 2), np.sqrt(c2 * dn_mag ** 2 + s2 * de_mag ** 2)
