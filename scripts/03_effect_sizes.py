#!/usr/bin/env python3
"""Step 3. Effect sizes: every stressed group against the control sampled on the SAME day.

The paper compares raw values across days (day 10 vs day 20). But the control beetles
themselves change a lot between day 10 and day 20, so a raw before/after comparison
mixes two things: the after-effect of stress and ordinary ageing. Dividing each stressed
group by its same-day control removes the ageing part.

Effect size: log response ratio  lnRR = ln(stressed / control), reported as % change.

Outputs (results/tables/)
  effect_sizes.csv         40 contrasts: 5 markers x 2 regimes x 2 phases x 2 sexes
  control_age_change.csv   control day 20 vs control day 10
  defence_damage.csv       antioxidant index (mean lnRR of SOD and CAT) vs LPO lnRR
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zbico.stats import log_response_ratio, hedges_g  # noqa: E402

TABLES = ROOT / "results" / "tables"
MARKERS = ["SOD", "CAT", "LPO", "glucose", "triglycerides"]


def holm(p):
    """Holm step-down adjusted p-values."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (len(p) - rank) * p[idx])
        adj[idx] = min(1.0, running)
    return adj


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(ROOT / "data" / "cell_summary.csv").set_index(["marker", "sex", "phase", "regime"])

    rows = []
    for marker in MARKERS:
        for sex in ["male", "female"]:
            for phase in ["stress", "recovery"]:
                c = d.loc[(marker, sex, phase, "control")]
                for regime in ["cold", "heat"]:
                    t = d.loc[(marker, sex, phase, regime)]
                    lnrr, se, lo, hi = log_response_ratio(t["mean"], t["se"], c["mean"], c["se"])
                    g, g_se = hedges_g(t["mean"], t["sd"], t["n"], c["mean"], c["sd"], c["n"])
                    rows.append({"marker": marker, "sex": sex, "phase": phase, "regime": regime,
                                 "mean_stressed": t["mean"], "mean_control": c["mean"],
                                 "lnRR": lnrr, "lnRR_se": se, "lnRR_lo": lo, "lnRR_hi": hi,
                                 "pct_change": 100 * (np.exp(lnrr) - 1),
                                 "pct_lo": 100 * (np.exp(lo) - 1), "pct_hi": 100 * (np.exp(hi) - 1),
                                 "hedges_g": g, "hedges_g_se": g_se,
                                 "p": 2 * stats.norm.sf(abs(lnrr / se))})
    es = pd.DataFrame(rows)
    es["p_holm"] = holm(es["p"])
    es["significant_holm"] = es["p_holm"] < 0.05
    es.round(5).to_csv(TABLES / "effect_sizes.csv", index=False)

    # how much do the controls themselves change between the two sampling days?
    age = []
    for marker in MARKERS:
        for sex in ["male", "female"]:
            c10 = d.loc[(marker, sex, "stress", "control")]
            c20 = d.loc[(marker, sex, "recovery", "control")]
            lnrr, se, lo, hi = log_response_ratio(c20["mean"], c20["se"], c10["mean"], c10["se"])
            age.append({"marker": marker, "sex": sex, "control_day10": c10["mean"], "control_day20": c20["mean"],
                        "lnRR": lnrr, "pct_change": 100 * (np.exp(lnrr) - 1),
                        "pct_lo": 100 * (np.exp(lo) - 1), "pct_hi": 100 * (np.exp(hi) - 1)})
    age = pd.DataFrame(age)
    age.round(5).to_csv(TABLES / "control_age_change.csv", index=False)

    # defence (SOD, CAT) against damage (LPO)
    w = es.pivot_table(index=["regime", "sex", "phase"], columns="marker", values="lnRR")
    v = es.assign(var=es.lnRR_se ** 2).pivot_table(index=["regime", "sex", "phase"], columns="marker", values="var")
    dd = pd.DataFrame({"defence_lnRR": (w["SOD"] + w["CAT"]) / 2,
                       "defence_se": np.sqrt(v["SOD"] + v["CAT"]) / 2,
                       "damage_lnRR": w["LPO"], "damage_se": np.sqrt(v["LPO"])}).reset_index()
    dd["defence_pct"] = 100 * (np.exp(dd.defence_lnRR) - 1)
    dd["damage_pct"] = 100 * (np.exp(dd.damage_lnRR) - 1)
    dd.round(5).to_csv(TABLES / "defence_damage.csv", index=False)

    show = es[["marker", "sex", "phase", "regime", "pct_change", "pct_lo", "pct_hi", "p_holm"]].copy()
    print("Percent change vs same-day control (95% CI), Holm-adjusted p:")
    print(show.round(3).to_string(index=False))
    print("\nControl beetles, day 20 vs day 10:")
    print(age[["marker", "sex", "pct_change", "pct_lo", "pct_hi"]].round(1).to_string(index=False))
    print("\nDefence (SOD+CAT) vs damage (LPO), % vs same-day control:")
    print(dd[["regime", "sex", "phase", "defence_pct", "damage_pct"]].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
