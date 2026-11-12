#!/usr/bin/env python3
"""Step 4. Mortality: reproduce the published chi-square, then analyse the counts.

Fig. 4 gives mortality as a percentage for six groups (3 regimes x 2 sexes); the group
sizes are in the Methods (500 per sex in each stress regime, 150 per sex in the control).

  a) Which whole-number percentages, consistent with the digitised bars, reproduce the
     published statistic (chi2 = 53.88, df = 5)?  -> pins down the percentages exactly.
  b) Chi-square on the actual numbers of beetles.
  c) Risk ratios against the control, and the male/female difference pooled over regimes
     (Mantel-Haenszel).

Outputs: results/tables/mortality_counts.csv, mortality_tests.csv
"""
from itertools import product
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zbico.stats import chi_square, risk_ratio, mantel_haenszel  # noqa: E402

TABLES = ROOT / "results" / "tables"
PUBLISHED_CHI2 = 53.88


def wilson(k, n, z=1.959964):
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(ROOT / "data" / "mortality.csv")

    # a) every whole-number combination within the digitisation error of each bar
    options = [(int(np.floor(v)), int(np.ceil(v))) for v in m.mortality_pct_digitised]
    matches, tried = [], 0
    for combo in product(*options):
        tried += 1
        chi2, _, _ = chi_square([combo, [100 - c for c in combo]])
        if round(chi2, 2) == PUBLISHED_CHI2:
            matches.append(combo)
    print(f"a) {tried} candidate sets of whole-number percentages; {len(matches)} reproduce chi2 = {PUBLISHED_CHI2}: {matches}")
    if len(matches) != 1:
        raise SystemExit("percentages are not uniquely identified - check the digitisation")
    m["mortality_pct"] = matches[0]
    chi_pct = chi_square([m.mortality_pct, 100 - m.mortality_pct])

    # b) counts. 15 % of 150 is 22.5, so the control-male percentage must be rounded: try 22 and 23
    m["dead_exact"] = m.mortality_pct / 100 * m.n_start
    tests = [{"test": "chi-square on percentages (n = 100 per group)", "chi2": chi_pct[0], "df": chi_pct[1], "p": chi_pct[2]}]
    for label, rounding in [("non-integer count rounded down", np.floor), ("non-integer count rounded up", np.ceil)]:
        dead = rounding(m.dead_exact.round(6)).astype(int)
        chi = chi_square([dead, m.n_start - dead])
        tests.append({"test": f"chi-square on counts ({label})", "chi2": chi[0], "df": chi[1], "p": chi[2]})
    m["dead"] = np.floor(m.dead_exact.round(6)).astype(int)        # conservative for the control
    m["alive"] = m.n_start - m.dead
    ci = [wilson(k, n) for k, n in zip(m.dead, m.n_start)]
    m["mortality_lo"] = [100 * c[0] for c in ci]
    m["mortality_hi"] = [100 * c[1] for c in ci]

    # c) risk ratios vs control, per sex
    g = m.set_index(["regime", "sex"])
    for regime in ["cold", "heat"]:
        for sex in ["male", "female"]:
            rr, lo, hi = risk_ratio(g.loc[(regime, sex), "dead"], g.loc[(regime, sex), "n_start"],
                                    g.loc[("control", sex), "dead"], g.loc[("control", sex), "n_start"])
            tests.append({"test": f"risk ratio {regime} vs control ({sex})", "estimate": rr, "lo": lo, "hi": hi})
    # heat vs cold directly
    for sex in ["male", "female"]:
        rr, lo, hi = risk_ratio(g.loc[("heat", sex), "dead"], g.loc[("heat", sex), "n_start"],
                                g.loc[("cold", sex), "dead"], g.loc[("cold", sex), "n_start"])
        tests.append({"test": f"risk ratio heat vs cold ({sex})", "estimate": rr, "lo": lo, "hi": hi})
    # male vs female, stratified by regime
    strata = [(g.loc[(r, "male"), "dead"], g.loc[(r, "male"), "alive"],
               g.loc[(r, "female"), "dead"], g.loc[(r, "female"), "alive"]) for r in ["cold", "control", "heat"]]
    or_mh, lo, hi, chi2, p = mantel_haenszel(strata)
    tests.append({"test": "male vs female, Mantel-Haenszel odds ratio across regimes",
                  "estimate": or_mh, "lo": lo, "hi": hi, "chi2": chi2, "df": 1, "p": p})

    m.drop(columns=["dead_exact"]).round(3).to_csv(TABLES / "mortality_counts.csv", index=False)
    t = pd.DataFrame(tests)
    t.to_csv(TABLES / "mortality_tests.csv", index=False)
    print(m[["regime", "sex", "n_start", "mortality_pct", "dead", "mortality_lo", "mortality_hi"]].round(1).to_string(index=False))
    pd.set_option("display.width", 200)
    print(t.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
