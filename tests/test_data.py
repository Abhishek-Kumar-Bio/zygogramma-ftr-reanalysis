"""Checks on the data files and on the pipeline outputs (run `python run_all.py` first)."""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_digitised_table_is_complete():
    d = pd.read_csv(ROOT / "data" / "digitised_means.csv")
    assert len(d) == 60
    assert d.groupby("marker").size().eq(12).all()
    assert d.groupby(["marker", "sex", "phase"]).regime.apply(lambda s: sorted(s) == ["cold", "control", "heat"]).all()
    assert (d["mean"] > 0).all() and (d["se"] > 0).all()


def test_calibrated_anova_reproduces_published_table():
    path = ROOT / "results" / "tables" / "anova_vs_published.csv"
    if not path.exists():
        pytest.skip("run the pipeline first")
    c = pd.read_csv(path)
    assert len(c) == 20
    main = c[c.used_for_calibration]
    assert main.ratio_calibrated_to_published.between(0.8, 1.25).all()
    held_out = c[~c.used_for_calibration]
    assert held_out.ratio_calibrated_to_published.between(0.7, 1.45).all()
