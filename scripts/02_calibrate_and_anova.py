#!/usr/bin/env python3
"""Step 2. Calibrate the digitised SEs against the published ANOVA, then recompute the ANOVA.

Why: the error bars in the published figures are only 1-6 px long, so reading them
off the image is the weakest part of the digitisation. The paper's Table 1, however,
gives exact F-ratios. For a factorial ANOVA,  F = MS_effect / MS_error,  and MS_effect
depends only on the group MEANS (which digitise well). So

        MS_error (implied) = MS_effect (from digitised means) / F (published)

gives the pooled within-group variance without touching the error bars. Comparing that
with the variance implied by the digitised whiskers tells us how biased the whisker
readings are. It turns out a single additive offset (in pixels) reconciles all five
markers, so that offset is fitted and applied.

Outputs
  data/cell_summary.csv                         means with raw and calibrated SE / SD
  results/tables/implied_error_variance.csv     MS_error implied by each published F
  results/tables/anova_recomputed.csv           full 3-way factorial ANOVA, all 7 effects
  results/tables/anova_vs_published.csv         recomputed F next to published F
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zbico.stats import factorial_anova_from_summary  # noqa: E402

FACTORS = ["regime", "phase", "sex"]
MAIN = ["regime", "phase", "sex"]
TABLES = ROOT / "results" / "tables"


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(ROOT / "data" / "digitised_means.csv")
    pub = pd.read_csv(ROOT / "data" / "published_anova.csv")

    # 1. error variance implied by each published F (uses digitised means only)
    d["sd_raw"] = d["se"] * np.sqrt(d["n"])
    rows = []
    for marker, cells in d.groupby("marker", sort=False):
        aov = factorial_anova_from_summary(cells, FACTORS, sd="sd_raw")
        for _, r in pub[pub.marker == marker].iterrows():
            rows.append({"marker": marker, "effect": r.effect, "F_published": r.F_published,
                         "MS_effect_from_means": aov.loc[r.effect, "MS"],
                         "MS_error_implied": aov.loc[r.effect, "MS"] / r.F_published,
                         "MS_error_raw_whiskers": aov.loc["error", "MS"]})
    implied = pd.DataFrame(rows)
    implied.to_csv(TABLES / "implied_error_variance.csv", index=False)
    # the three main effects have large F and stable MS; the interaction is kept as a hold-out check
    target = implied[implied.effect.isin(MAIN)].groupby("marker", sort=False)["MS_error_implied"].median()

    # 2. one additive whisker offset (pixels) for all panels
    def pooled_var(cells, offset):
        return float(((((cells["se_px"] + offset) * cells["unit_per_px"]) ** 2) * cells["n"]).mean())

    def loss(offset):
        return sum(np.log(pooled_var(d[d.marker == m], offset) / target[m]) ** 2 for m in target.index)

    fit = minimize_scalar(loss, bounds=(0, 3), method="bounded")
    offset = float(fit.x)
    d["se_raw"] = d["se"]
    d["se"] = (d["se_px"] + offset) * d["unit_per_px"]
    d["sd"] = d["se"] * np.sqrt(d["n"])
    cols = ["marker", "sex", "regime", "phase", "day", "mean", "se", "sd", "n", "unit", "se_raw", "se_px", "unit_per_px", "source"]
    d[cols].round(6).to_csv(ROOT / "data" / "cell_summary.csv", index=False)

    # 3. recompute the ANOVA with calibrated SDs and compare with Table 1
    full, comp = [], []
    for marker, cells in d.groupby("marker", sort=False):
        aov = factorial_anova_from_summary(cells, FACTORS).reset_index(names="effect")
        aov.insert(0, "marker", marker)
        full.append(aov)
        raw = factorial_anova_from_summary(cells, FACTORS, sd="sd_raw")
        for _, r in pub[pub.marker == marker].iterrows():
            f_new = float(aov.loc[aov.effect == r.effect, "F"].iloc[0])
            comp.append({"marker": marker, "effect": r.effect, "F_published": r.F_published,
                         "F_raw_whiskers": raw.loc[r.effect, "F"], "F_calibrated": f_new,
                         "ratio_calibrated_to_published": f_new / r.F_published,
                         "used_for_calibration": r.effect in MAIN})
    full = pd.concat(full)
    comp = pd.DataFrame(comp)
    full.to_csv(TABLES / "anova_recomputed.csv", index=False)
    comp.to_csv(TABLES / "anova_vs_published.csv", index=False)

    print(f"Fitted whisker offset: {offset:.2f} px")
    print("\nPooled error variance, raw whiskers vs implied by published F (ratio implied/raw):")
    for m in target.index:
        raw = pooled_var(d[d.marker == m], 0.0)
        print(f"  {m:14s} raw {raw:.4g}   implied {target[m]:.4g}   ratio {target[m] / raw:.2f}"
              f"   after calibration {pooled_var(d[d.marker == m], offset) / target[m]:.2f}")
    print("\nRecomputed vs published F:")
    print(comp.round(2).to_string(index=False))
    held = comp[~comp.used_for_calibration]
    print(f"\nHold-out (3-way interaction) ratio calibrated/published: "
          f"median {held.ratio_calibrated_to_published.median():.2f}, "
          f"range {held.ratio_calibrated_to_published.min():.2f}-{held.ratio_calibrated_to_published.max():.2f}")
    print("\nEffects not reported in the paper (two-way interactions):")
    two = full[full.effect.isin(["regime x phase", "regime x sex", "phase x sex"])]
    print(two[["marker", "effect", "df", "F", "p"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
