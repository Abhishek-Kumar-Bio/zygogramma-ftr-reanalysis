"""Statistics that work from summary data (group means, standard errors, counts).

Everything here is deliberately written out with numpy/scipy so each formula can be
read and checked. The tests in tests/ compare these functions against raw-data
implementations and textbook values.
"""
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------- effect sizes
def log_response_ratio(mean_t, se_t, mean_c, se_c, level=0.95):
    """Log response ratio lnRR = ln(treatment / control) with a delta-method SE.

    var(lnRR) = SE_t^2 / mean_t^2 + SE_c^2 / mean_c^2     (Hedges, Gurevitch & Curtis 1999)
    Returns (lnRR, SE, lower, upper). exp(lnRR) - 1 is the proportional change.
    """
    lnrr = np.log(mean_t / mean_c)
    se = np.sqrt((se_t / mean_t) ** 2 + (se_c / mean_c) ** 2)
    z = stats.norm.ppf(0.5 + level / 2)
    return lnrr, se, lnrr - z * se, lnrr + z * se


def hedges_g(mean_t, sd_t, n_t, mean_c, sd_c, n_c):
    """Standardised mean difference with the small-sample correction J. Returns (g, SE)."""
    df = n_t + n_c - 2
    sd_pooled = np.sqrt(((n_t - 1) * sd_t ** 2 + (n_c - 1) * sd_c ** 2) / df)
    d = (mean_t - mean_c) / sd_pooled
    j = 1 - 3 / (4 * df - 1)
    g = j * d
    se = np.sqrt((n_t + n_c) / (n_t * n_c) + g ** 2 / (2 * (n_t + n_c)))
    return g, se


# ---------------------------------------------------- factorial ANOVA from cells
def factorial_anova_from_summary(cells, factors, mean="mean", sd="sd", n="n"):
    """Between-subjects factorial ANOVA computed from cell means, SDs and a common n.

    `cells` has one row per cell of a fully crossed, balanced design.
    Sums of squares for each effect come from the cell means (inclusion-exclusion over
    marginal means); the error mean square is the average within-cell variance.
    Returns a DataFrame indexed by effect with SS, df, MS, F, p.
    """
    ns = cells[n].unique()
    if len(ns) != 1:
        raise ValueError("balanced design required (same n in every cell)")
    n_per = int(ns[0])
    levels = {f: cells[f].nunique() for f in factors}
    if len(cells) != int(np.prod(list(levels.values()))):
        raise ValueError("design is not fully crossed")
    grand = cells[mean].mean()
    n_cells = len(cells)

    def ss_of_margin(subset):
        """n_obs_per_margin_cell * sum((marginal mean - grand)^2) for a set of factors."""
        marg = cells.groupby(list(subset))[mean].mean()
        obs_per = n_per * n_cells / len(marg)
        return obs_per * float(((marg - grand) ** 2).sum())

    rows = {}
    ss_effect = {}
    for k in range(1, len(factors) + 1):
        for subset in combinations(factors, k):
            ss = ss_of_margin(subset)
            for j in range(1, k):                       # remove all lower-order effects
                for lower in combinations(subset, j):
                    ss -= ss_effect[lower]
            ss_effect[subset] = ss
            df = int(np.prod([levels[f] - 1 for f in subset]))
            rows[" x ".join(subset)] = {"SS": ss, "df": df}
    df_error = n_cells * (n_per - 1)
    ms_error = float((cells[sd] ** 2).mean())
    out = pd.DataFrame(rows).T
    out["MS"] = out["SS"] / out["df"]
    out["F"] = out["MS"] / ms_error
    out["p"] = stats.f.sf(out["F"], out["df"], df_error)
    out.loc["error"] = {"SS": ms_error * df_error, "df": df_error, "MS": ms_error, "F": np.nan, "p": np.nan}
    out["df"] = out["df"].astype(int)
    return out


# ------------------------------------------------------------------ count data
def chi_square(table):
    """Pearson chi-square for an r x c table of counts (no continuity correction).
    Returns (chi2, df, p)."""
    t = np.asarray(table, float)
    expected = t.sum(axis=1, keepdims=True) * t.sum(axis=0, keepdims=True) / t.sum()
    chi2 = float(((t - expected) ** 2 / expected).sum())
    df = (t.shape[0] - 1) * (t.shape[1] - 1)
    return chi2, df, float(stats.chi2.sf(chi2, df))


def risk_ratio(dead_t, n_t, dead_c, n_c, level=0.95):
    """Risk ratio of death, treatment vs control, with a log-scale Wald interval."""
    rr = (dead_t / n_t) / (dead_c / n_c)
    se = np.sqrt(1 / dead_t - 1 / n_t + 1 / dead_c - 1 / n_c)
    z = stats.norm.ppf(0.5 + level / 2)
    return rr, float(np.exp(np.log(rr) - z * se)), float(np.exp(np.log(rr) + z * se))


def mantel_haenszel(strata, level=0.95):
    """Common odds ratio across 2x2 strata and the Cochran-Mantel-Haenszel test.

    Each stratum is (a, b, c, d) = (exposed dead, exposed alive, unexposed dead, unexposed alive).
    Returns (OR_MH, lower, upper, chi2_CMH, p). Interval: Robins, Breslow & Greenland (1986).
    """
    s = np.asarray(strata, float)
    a, b, c, d = s.T
    n = a + b + c + d
    or_mh = float((a * d / n).sum() / (b * c / n).sum())
    expected = (a + b) * (a + c) / n
    var = (a + b) * (c + d) * (a + c) * (b + d) / (n ** 2 * (n - 1))
    chi2 = float((a.sum() - expected.sum()) ** 2 / var.sum())      # no continuity correction
    p_, q_ = (a + d) / n, (b + c) / n
    r_, s_ = a * d / n, b * c / n
    var_log = float((p_ * r_).sum() / (2 * r_.sum() ** 2)
                    + (p_ * s_ + q_ * r_).sum() / (2 * r_.sum() * s_.sum())
                    + (q_ * s_).sum() / (2 * s_.sum() ** 2))
    z = stats.norm.ppf(0.5 + level / 2)
    lo, hi = np.exp(np.log(or_mh) - z * np.sqrt(var_log)), np.exp(np.log(or_mh) + z * np.sqrt(var_log))
    return or_mh, float(lo), float(hi), chi2, float(stats.chi2.sf(chi2, 1))
