#!/usr/bin/env python3
"""Step 6 (extension, to be filled in by hand). Put Z. bicolorata next to other insects.

data/literature/extracted_means.csv holds one row per comparison: a stressed group and its
control, taken from a published table or digitised from a figure (WebPlotDigitizer works
well: https://automeris.io). The Z. bicolorata rows are written automatically from this
repository's own data; every other row has to be extracted and checked by a person.

The script computes the log response ratio for each row and, once at least three studies
are present, draws a forest plot per marker.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zbico.stats import log_response_ratio  # noqa: E402

LIT = ROOT / "data" / "literature" / "extracted_means.csv"
COLUMNS = ["study", "doi", "species", "order", "stressor", "protocol", "temperature_C", "exposure",
           "marker", "sex_or_stage", "mean_control", "se_control", "n_control",
           "mean_stress", "se_stress", "n_stress", "source", "extracted_by", "notes"]


def own_rows():
    d = pd.read_csv(ROOT / "data" / "cell_summary.csv")
    d = d[(d.phase == "stress") & d.marker.isin(["SOD", "CAT", "LPO"])].set_index(["marker", "sex", "regime"])
    rows = []
    for marker in ["SOD", "CAT", "LPO"]:
        for sex in ["male", "female"]:
            c = d.loc[(marker, sex, "control")]
            for regime, temp in (("heat", 42), ("cold", 10)):
                t = d.loc[(marker, sex, regime)]
                rows.append({"study": "Patel et al. 2025", "doi": "10.4039/tce.2025.10012",
                             "species": "Zygogramma bicolorata", "order": "Coleoptera", "stressor": regime,
                             "protocol": "repeated", "temperature_C": temp, "exposure": "1 h/day x 10 days",
                             "marker": marker, "sex_or_stage": f"adult {sex}",
                             "mean_control": c["mean"], "se_control": c["se"], "n_control": 10,
                             "mean_stress": t["mean"], "se_stress": t["se"], "n_stress": 10,
                             "source": c["source"], "extracted_by": "digitised by script (step 1)",
                             "notes": "SE calibrated against published F (step 2)"})
    return pd.DataFrame(rows, columns=COLUMNS)


def main():
    own = own_rows()
    if LIT.exists():
        lit = pd.read_csv(LIT)
        lit = lit[lit.study != "Patel et al. 2025"]
        table = pd.concat([own, lit], ignore_index=True)
    else:
        table = own
    table.to_csv(LIT, index=False)

    res = table.copy()
    out = [log_response_ratio(r.mean_stress, r.se_stress, r.mean_control, r.se_control) for r in res.itertuples()]
    res["lnRR"], res["lnRR_se"], res["lnRR_lo"], res["lnRR_hi"] = zip(*out)
    res["pct_change"] = 100 * (np.exp(res.lnRR) - 1)
    res.round(5).to_csv(ROOT / "results" / "tables" / "cross_species_effects.csv", index=False)

    n_studies = res.study.nunique()
    print(f"{len(res)} comparisons from {n_studies} stud{'y' if n_studies == 1 else 'ies'}.")
    if n_studies < 3:
        print("Add at least two more studies to data/literature/extracted_means.csv to get the forest plot.")
        print("Candidate papers are listed in docs/EXTENSION_cross_species.md.")
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colour = {"heat": "#eb6834", "cold": "#2a78d6"}
    markers = [m for m in ["SOD", "CAT", "LPO"] if m in set(res.marker)]
    fig, axes = plt.subplots(1, len(markers), figsize=(4.2 * len(markers), 0.32 * len(res) / len(markers) + 2), sharey=False)
    for ax, marker in zip(np.atleast_1d(axes), markers):
        sub = res[res.marker == marker].sort_values(["stressor", "lnRR"]).reset_index(drop=True)
        ax.axvline(0, color="#52514e", lw=0.8)
        for i, r in sub.iterrows():
            ax.plot([r.lnRR_lo, r.lnRR_hi], [i, i], color=colour.get(r.stressor, "#6b6a66"), lw=2)
            ax.plot(r.lnRR, i, "o", color=colour.get(r.stressor, "#6b6a66"), ms=6)
        ax.set_yticks(range(len(sub)), [f"{r.species} ({r.sex_or_stage}, {r.stressor})" for r in sub.itertuples()], fontsize=7)
        ax.set_title(marker, loc="left")
        ax.set_xlabel("ln(stressed / control)")
    fig.tight_layout()
    fig.savefig(ROOT / "results" / "figures" / "cross_species_forest.png", dpi=200)
    print("wrote results/figures/cross_species_forest.png")


if __name__ == "__main__":
    main()
