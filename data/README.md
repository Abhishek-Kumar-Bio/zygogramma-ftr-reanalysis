# Data

All values come from one published paper:

> Patel, A.K., Yadav, P. & Kumar, B. (2025). Physiological responses of Parthenium beetles,
> *Zygogramma bicolorata* Pallister (Coleoptera: Chrysomelidae), to fluctuating temperature
> stress. *The Canadian Entomologist* 157: e26. https://doi.org/10.4039/tce.2025.10012

The paper reports group means as bar charts and does not provide raw data, so the means and
standard errors here were measured from the figure images (`scripts/01_digitise_figures.py`).

| File | What it is | How it was obtained |
|---|---|---|
| `digitised_means.csv` | mean and SE for 5 markers x 3 regimes x 2 phases x 2 sexes (60 groups, n = 10 each) | digitised from Figs 2 and 3 |
| `mortality.csv` | mortality (%) for 3 regimes x 2 sexes, with group sizes | digitised from Fig. 4; group sizes from the Methods |
| `published_anova.csv` | the 20 F-ratios in Table 1 | typed in from the paper |
| `cell_summary.csv` | `digitised_means.csv` plus calibrated SE and SD | written by `scripts/02_calibrate_and_anova.py` |
| `literature/extracted_means.csv` | table for the cross-species extension | Z. bicolorata rows written by step 6; other rows to be added by hand |

Column notes for `digitised_means.csv`
- `regime`: `cold` (10 °C for 1 h/day), `control` (27 °C throughout), `heat` (42 °C for 1 h/day).
- `phase`: `stress` = sampled after the 10 days of daily exposure; `recovery` = sampled after a
  further 10 days at 27 °C.
- `unit`: as printed on the figure axis. (The Methods give glucose and triglycerides in mg/dL;
  the axes say mg/ml. Ratios against the control are unaffected.)
- `se_px`: length of the error bar in pixels. `unit_per_px`: size of one pixel in data units.

Accuracy: means are good to about 1 % of the axis range. Raw SEs are coarse (most error bars are
1-6 px long); see `docs/METHODS.md` for how they are calibrated.
