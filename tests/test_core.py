"""Checks of the pieces that need no data archive, and one that the FPI
statistics are reproduced exactly when FPI_ForwardModel's outputs exist."""

import numpy as np
import pandas as pd
import pytest

from isr_clim import climatology, config, fpi


def test_pairing_takes_the_step_within_30_min_and_never_interpolates():
    model = pd.DataFrame({"time": pd.to_datetime(["2012-01-01 00:00", "2012-01-01 03:00"]),
                          "cgm_lat": 83.0, "vn_geo": [0.0, 300.0], "ve_geo": [10.0, 20.0]})
    rec = pd.DataFrame({"time": pd.to_datetime(["2012-01-01 00:20", "2012-01-01 01:30", "2012-01-01 02:31",
                                                "2012-01-01 07:00"]), "cgm_lat": 83.0})
    p = climatology.pair(rec, model)
    np.testing.assert_array_equal(p["model_vn_geo"].to_numpy(), [0.0, np.nan, 300.0, np.nan])
    assert p.loc[2, "model_time"] == pd.Timestamp("2012-01-01 03:00")
    assert p["model_time"].isna().tolist() == [False, True, False, True]


def test_pairing_is_per_bin():
    model = pd.DataFrame({"time": pd.to_datetime(["2012-01-01 00:00"] * 2), "cgm_lat": [83.0, 83.25],
                          "vn_geo": [1.0, 2.0], "ve_geo": [0.0, 0.0]})
    rec = pd.DataFrame({"time": pd.to_datetime(["2012-01-01 00:05"] * 2), "cgm_lat": [83.25, 83.0]})
    assert climatology.pair(rec, model)["model_vn_geo"].tolist() == [2.0, 1.0]


def test_obs_mapped_steps_count_each_step_once():
    t0 = pd.Timestamp("2012-01-01 00:00")
    p = pd.DataFrame({"cgm_lat": [83.0, 83.0, 83.25, 83.0], "model_time": [t0, t0, t0, pd.NaT],
                      "model_vn_geo": [5.0, 5.0, 7.0, np.nan]})
    _, v = climatology.obs_mapped_steps(p, "vn_geo")
    assert sorted(v) == [5.0, 7.0]


def test_diurnal_weights_and_local_hour_follow_fpi():
    tz = "America/Resolute"                     # CST (UTC-6) in January
    t = pd.to_datetime(["2012-01-10 06:10", "2012-01-10 06:20", "2012-01-10 06:25"])
    d = climatology.seasonal_diurnal(t, [100.0, 200.0, 50.0], [10.0, 20.0, np.nan], [150.0, 150.0, np.nan],
                                     pd.to_datetime(["2012-01-10 06:00"]), [150.0], tz)
    row = d.table[(d.table["season"] == "DJF") & (d.table["hour"] == 0)].iloc[0]
    w = np.array([1 / 100, 1 / 400, 1 / 0.25])   # a missing error takes the 0.5 m/s floor, as in FPI
    assert np.isclose(row["obs_mean"], np.sum(w * [100, 200, 50]) / w.sum())
    assert np.isclose(row["obs_err"], np.sqrt(1 / w.sum()))
    assert np.isclose(row["res_mean"], (50 * w[0] - 50 * w[1]) / (w[0] + w[1]))
    assert row["res_n"] == 2 and row["model_n"] == 1 and row["model_mean"] == 150.0


def test_scores_match_definitions_and_need_five_pairs():
    t = pd.date_range("2012-01-01 06:00", periods=6, freq="min")
    obs, mod = np.array([1.0, 2, 3, 4, 5, 6]), np.array([2.0, 2, 5, 4, 7, 9])
    d = climatology.seasonal_diurnal(t, obs, np.ones(6), mod, t[:1], mod[:1], "UTC")
    s = climatology.scores(d, {}).iloc[0]
    assert s["N"] == 6 and np.isclose(s["Bias (m/s)"], round((mod - obs).mean(), 1))
    assert np.isclose(s["r"], round(np.corrcoef(obs, mod)[0, 1], 3))
    assert climatology.scores(climatology.seasonal_diurnal(t[:4], obs[:4], np.ones(4), mod[:4], t[:1], mod[:1],
                                                           "UTC"), {}).empty


def test_layer_mean_is_height_weighted_and_stops_at_model_top():
    waccmxsd = pytest.importorskip("isr_clim.waccmxsd")    # needs aacgmv2 (via geometry)
    z = np.array([100.0, 200.0, 300.0, 500.0])             # uneven levels, top below 700 km
    v = 2.0 * z                                            # linear in height
    out = waccmxsd.layer_mean(z, v, (175.0, 700.0), 5.0)
    assert np.isclose(out, 2.0 * np.mean(np.arange(175.0, 500.0 + 1e-9, 5.0)))
    assert np.isclose(waccmxsd.layer_mean(z, np.full(4, 7.0)), 7.0)


def test_frame_rotation_and_error_propagation():
    geometry = pytest.importorskip("isr_clim.geometry")
    e, n = geometry.mag_to_geo(0.0, 1.0, -20.0)            # unit magnetic-north vector
    assert np.isclose(np.degrees(np.arctan2(e, n)), -20.0)
    de, dn = geometry.mag_to_geo_errors(np.array([10.0]), np.array([30.0]), np.array([90.0]))
    assert np.allclose([de[0], dn[0]], [30.0, 10.0])        # 90 deg swaps the axes
    de, dn = geometry.mag_to_geo_errors(10.0, 30.0, 37.0)
    assert np.isclose(de ** 2 + dn ** 2, 10.0 ** 2 + 30.0 ** 2)


def test_risrn_qc_drops_degenerate_one_metre_per_second_errors(tmp_path):
    risrn = pytest.importorskip("isr_clim.risrn")          # needs aacgmv2 (via geometry)
    import h5py

    names = ("ut1_unix", "ut2_unix", "cgm_lat", "vi1", "vi2", "dvi1", "dvi2")
    rows = [(0.0, 180.0, 83.0, 100.0, 200.0, 30.0, 40.0),     # kept
            (0.0, 180.0, 83.0, 0.0, -2537.0, 1.0, 1.0),       # degenerate fit
            (0.0, 180.0, 83.0, 50.0, 60.0, 160.0, 40.0),      # error too large
            (0.0, 180.0, 83.0, 3000.0, 100.0, 30.0, 40.0)]    # too fast
    with h5py.File(tmp_path / "x.hdf5", "w") as h:
        h["Data/Table Layout"] = np.array(rows, dtype=[(n, "f8") for n in names])
    df = risrn.load_file(str(tmp_path / "x.hdf5"))
    assert df[["ve_geo", "vn_geo", "dve_geo", "dvn_geo"]].values.tolist() == [[100.0, 200.0, 30.0, 40.0]]


def test_fpi_statistics_are_reproduced():
    """Given the FPI data, this package's pairing and scores give exactly the
    statistics of FPI's own plot_seasonal_diurnal_climatology (the current
    code, not a stored wind_stats CSV, which may predate it)."""
    try:
        outputs = fpi.load_fpi_outputs()
    except FileNotFoundError:
        pytest.skip("FPI outputs for the station not available")
    import climatology_plots as cp          # FPI_ForwardModel/Plotting, on the path via isr_clim.fpi
    import matplotlib.pyplot as plt

    ds_obs, sim, elevation = outputs
    st = fpi.station()
    models = cp.build_models([(config.MODEL_KEY, sim, config.FPI_COMBO)])
    for comp, n in fpi.fpi_meridional_zonal(outputs=outputs).items():
        fig, fig_nr, want = cp.plot_seasonal_diurnal_climatology(ds_obs, models, comp, elevation,
                                                                 st.plot_name, st.lat, st.tz)
        plt.close(fig)
        plt.close(fig_nr)
        d = climatology.seasonal_diurnal(n["time"], n["obs"], n["err"], n["model_at_obs"],
                                         n["model_time"], n["model_value"], st.tz)
        got = climatology.scores(d, {}).set_index("Season")
        want = want.set_index("Season")
        assert set(got.index) == set(want.index)
        for col in ("N", "Bias (m/s)", "RMSE (m/s)", "MAE (m/s)", "r"):
            np.testing.assert_allclose(got.loc[want.index, col].astype(float), want[col].astype(float),
                                       atol=1e-9, err_msg=f"{comp} {col}")
