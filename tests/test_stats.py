"""Checks that the summary-statistics functions agree with raw-data methods and known values."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zbico.stats import (log_response_ratio, hedges_g, factorial_anova_from_summary,  # noqa: E402
                         chi_square, risk_ratio, mantel_haenszel)


def simulate(seed=1, n=10):
    rng = np.random.default_rng(seed)
    rows = []
    for a in ["cold", "control", "heat"]:
        for b in ["stress", "recovery"]:
            for c in ["male", "female"]:
                mu = 10 + {"cold": -2, "control": 0, "heat": 3}[a] + (1.5 if b == "recovery" else 0) \
                     + (0.8 if c == "female" else 0) + (2 if (a == "heat" and b == "recovery") else 0)
                for y in rng.normal(mu, 1.3, n):
                    rows.append({"regime": a, "phase": b, "sex": c, "y": y})
    raw = pd.DataFrame(rows)
    cells = raw.groupby(["regime", "phase", "sex"])["y"].agg(mean="mean", sd="std", n="count").reset_index()
    return raw, cells


def test_lnrr_known_value():
    lnrr, se, lo, hi = log_response_ratio(120.0, 6.0, 100.0, 5.0)
    assert lnrr == pytest.approx(np.log(1.2))
    assert se == pytest.approx(np.sqrt(0.05 ** 2 + 0.05 ** 2))
    assert lo < lnrr < hi
    # identical groups -> no effect
    assert log_response_ratio(5, 1, 5, 1)[0] == 0


def test_hedges_g_direction_and_correction():
    g, se = hedges_g(12, 2, 10, 10, 2, 10)
    d = 1.0                                    # (12 - 10) / 2
    assert g == pytest.approx(d * (1 - 3 / (4 * 18 - 1)))
    assert se > 0


def test_anova_one_factor_matches_scipy():
    rng = np.random.default_rng(3)
    groups = [rng.normal(m, 1, 12) for m in (0, 0.5, 1.5)]
    cells = pd.DataFrame({"g": ["a", "b", "c"], "mean": [x.mean() for x in groups],
                          "sd": [x.std(ddof=1) for x in groups], "n": 12})
    aov = factorial_anova_from_summary(cells, ["g"])
    f, p = stats.f_oneway(*groups)
    assert aov.loc["g", "F"] == pytest.approx(f)
    assert aov.loc["g", "p"] == pytest.approx(p)


def test_anova_sums_of_squares_add_up():
    raw, cells = simulate()
    aov = factorial_anova_from_summary(cells, ["regime", "phase", "sex"])
    ss_total = ((raw.y - raw.y.mean()) ** 2).sum()
    assert aov["SS"].sum() == pytest.approx(ss_total)
    assert aov.loc["error", "df"] == 108
    assert aov["df"].sum() == 119


def test_anova_three_factor_matches_statsmodels():
    smf = pytest.importorskip("statsmodels.formula.api")
    sm = pytest.importorskip("statsmodels.api")
    raw, cells = simulate(seed=7)
    aov = factorial_anova_from_summary(cells, ["regime", "phase", "sex"])
    ref = sm.stats.anova_lm(smf.ols("y ~ C(regime) * C(phase) * C(sex)", data=raw).fit(), typ=2)
    pairs = {"regime": "C(regime)", "phase": "C(phase)", "sex": "C(sex)", "regime x phase": "C(regime):C(phase)",
             "regime x sex": "C(regime):C(sex)", "phase x sex": "C(phase):C(sex)",
             "regime x phase x sex": "C(regime):C(phase):C(sex)"}
    for mine, theirs in pairs.items():
        assert aov.loc[mine, "F"] == pytest.approx(ref.loc[theirs, "F"], rel=1e-8)


def test_chi_square_matches_scipy():
    table = [[135, 100, 22, 15, 235, 200], [365, 400, 128, 135, 265, 300]]
    chi2, df, p = chi_square(table)
    ref = stats.chi2_contingency(table, correction=False)
    assert chi2 == pytest.approx(ref[0])
    assert df == ref[2] == 5
    assert p == pytest.approx(ref[1])


def test_published_chi_square_is_reproduced_from_percentages():
    pct = [27, 20, 15, 10, 47, 40]
    chi2, df, _ = chi_square([pct, [100 - v for v in pct]])
    assert round(chi2, 2) == 53.88 and df == 5


def test_risk_ratio():
    rr, lo, hi = risk_ratio(200, 500, 15, 150)
    assert rr == pytest.approx(4.0)
    assert lo < rr < hi


def test_mantel_haenszel_single_stratum():
    a, b, c, d = 135, 365, 100, 400
    or_mh, lo, hi, chi2, p = mantel_haenszel([(a, b, c, d)])
    assert or_mh == pytest.approx(a * d / (b * c))
    n = a + b + c + d
    pearson = chi_square([[a, b], [c, d]])[0]
    assert chi2 == pytest.approx(pearson * (n - 1) / n)


def test_mantel_haenszel_matches_statsmodels():
    sm = pytest.importorskip("statsmodels.api")
    strata = [(135, 365, 100, 400), (22, 128, 15, 135), (235, 265, 200, 300)]
    or_mh, lo, hi, chi2, p = mantel_haenszel(strata)
    st = sm.stats.StratifiedTable([np.array([[a, b], [c, d]]) for a, b, c, d in strata])
    assert or_mh == pytest.approx(st.oddsratio_pooled)
    assert chi2 == pytest.approx(st.test_null_odds(correction=False).statistic)
    ref_lo, ref_hi = st.oddsratio_pooled_confint()
    assert lo == pytest.approx(ref_lo, rel=1e-6) and hi == pytest.approx(ref_hi, rel=1e-6)
