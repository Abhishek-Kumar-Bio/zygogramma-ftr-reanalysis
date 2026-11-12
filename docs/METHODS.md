# Methods

## 1. Digitising the figures
`scripts/01_digitise_figures.py` reads each bar chart as an image.
- The y-axis tick marks are located and a straight line is fitted through all of them, which
  converts pixel rows to data values.
- A bar is a run of saturated colour rising from the zero line; its height is measured with
  sub-pixel coverage at the top edge.
- The SE is the length of the upper whisker above the bar top.

Checks
- Mortality: of the 64 sets of whole-number percentages compatible with the digitised bars,
  exactly one (27, 20, 15, 10, 47, 40 %) reproduces the published chi-square of 53.88.
- Means: the effect mean squares computed from the digitised means, divided by the published
  F-ratios, give almost the same error variance for all four effects reported for each marker
  (`results/tables/implied_error_variance.csv`). That can only happen if the digitised means
  are close to the ones behind the published table.

## 2. Calibrating the standard errors
The error bars are 1-6 px long, so reading them is the weak point. But for a factorial ANOVA

    F = MS_effect / MS_error      ->      MS_error = MS_effect / F

and MS_effect needs only the means. So the published F-ratios give the pooled within-group
variance directly. The whisker readings turned out to be too short by a constant amount: one
additive offset of 0.83 px, fitted to the 15 main-effect F-ratios, brings the pooled variance of
all five markers to within 7 % of the implied value. The five three-way-interaction F-ratios
were not used in the fit; after calibration they are reproduced to within 0.87-1.31x.

What this assumes: the bias in reading a whisker is the same number of pixels in every panel.
What it does not assume: equal variances across groups. Each group keeps its own SE.

## 3. ANOVA from summary statistics
`zbico.stats.factorial_anova_from_summary` computes a between-subjects factorial ANOVA from
cell means, SDs and n. For a balanced design the sums of squares follow from the marginal
means, and the error mean square is the average within-cell variance. The tests compare it with
`scipy.stats.f_oneway` and with `statsmodels` on simulated raw data.

Notes on the comparison with the published Table 1
- The paper calls the analysis a repeated-measures ANOVA. Each beetle was homogenised for the
  assays, so the day-10 and day-20 samples are different individuals, and exposure duration is
  a between-subjects factor. The published F-ratios are consistent with an ordinary factorial
  ANOVA with one error term, which is what is recomputed here.
- Table 1 lists denominator degrees of freedom of 119. With 120 observations in 12 groups the
  error term has 108. The F-ratios are not affected.
- Table 1 reports the three main effects and the three-way interaction. The two-way
  interactions are in `results/tables/anova_recomputed.csv`; regime x phase is large for every
  marker, which is the statistical form of "the effect of stress depends on when you look".

## 4. Effect sizes
Each stressed group is compared with the control sampled on the same day:

    lnRR = ln(mean_stressed / mean_control),   var = SE_s^2/mean_s^2 + SE_c^2/mean_c^2

reported as a percentage change with a 95 % interval. p-values are Holm-adjusted over the 40
contrasts. Hedges' g is also in the table.

Why same-day controls: the control beetles themselves differ between day 10 and day 20
(`results/tables/control_age_change.csv`). Comparing a stressed group on day 20 with itself on
day 10 therefore mixes the after-effect of stress with whatever changed in all beetles over
those ten days (age, and the fact that the two days were assayed in separate runs).

## 5. Mortality
The published chi-square is reproduced by running the test on the six percentages, i.e. as if
each group had 100 beetles. The group sizes were 500, 500 and 150 per sex, so the test is rerun
on counts. 15 % of 150 is 22.5, so the control-male count is ambiguous (22 or 23); both are
reported. Risk ratios use a log-scale Wald interval; the sex difference is pooled over the three
regimes with the Mantel-Haenszel estimator.

## Limits of this re-analysis
- It uses 60 group means, not the 600 individual measurements. Individual-level questions
  (outliers, distributions, correlations between markers within a beetle) cannot be addressed.
- Intervals use a normal approximation with n = 10 per group.
- SEs depend on the calibration in section 2.
- One species, one population, one laboratory.
