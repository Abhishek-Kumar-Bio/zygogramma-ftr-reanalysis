#!/usr/bin/env python3
"""Step 5. Figures (written to results/figures/). Run steps 1-4 first."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3de"
COLOUR = {"cold": "#2a78d6", "heat": "#eb6834", "control": "#6b6a66"}     # colour = thermal regime
SHAPE = {"male": "o", "female": "D"}                                      # shape  = sex
MARKERS = ["SOD", "CAT", "LPO", "glucose", "triglycerides"]
LABEL = {"SOD": "SOD activity", "CAT": "Catalase activity", "LPO": "Lipid peroxidation",
         "glucose": "Glucose", "triglycerides": "Triglycerides"}
PCT_TICKS = [-75, -50, -25, 0, 25, 50, 100]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False,
})


def pct_axis(ax, axis="x", ticks=PCT_TICKS):
    pos = [np.log(1 + t / 100) for t in ticks]
    lab = [f"{t:+d}%" if t else "0" for t in ticks]
    if axis == "x":
        ax.set_xticks(pos, lab)
    else:
        ax.set_yticks(pos, lab)


def legend_handles(regimes=("cold", "heat")):
    name = {"cold": "Repeated cold (10 °C, 1 h/day)", "heat": "Repeated heat (42 °C, 1 h/day)", "control": "Control (27 °C)"}
    h = [Line2D([], [], marker="s", ls="", color=COLOUR[r], ms=7, label=name[r]) for r in regimes]
    h += [Line2D([], [], marker=SHAPE[s], ls="", color=INK2, ms=6, label=s.capitalize()) for s in ("male", "female")]
    return h


def fig_effect_sizes():
    es = pd.read_csv(TABLES / "effect_sizes.csv")
    fig, axes = plt.subplots(5, 2, figsize=(8.2, 8.4), sharex=True)
    order = [("cold", "male"), ("cold", "female"), ("heat", "male"), ("heat", "female")]
    for i, marker in enumerate(MARKERS):
        for j, phase in enumerate(["stress", "recovery"]):
            ax = axes[i, j]
            ax.axvline(0, color=INK2, lw=0.8)
            ax.grid(axis="x", color=GRID, lw=0.6)
            ax.set_axisbelow(True)
            for k, (regime, sex) in enumerate(order):
                r = es[(es.marker == marker) & (es.phase == phase) & (es.regime == regime) & (es.sex == sex)].iloc[0]
                y = 3 - k
                ax.plot([r.lnRR_lo, r.lnRR_hi], [y, y], color=COLOUR[regime], lw=2, solid_capstyle="round")
                ax.plot(r.lnRR, y, SHAPE[sex], ms=7, mec=COLOUR[regime], mew=1.6,
                        mfc=COLOUR[regime] if r.significant_holm else SURFACE, zorder=3)
            ax.set_ylim(-0.6, 3.6)
            ax.set_yticks([3, 2, 1, 0], ["Cold, male", "Cold, female", "Heat, male", "Heat, female"] if j == 0 else [""] * 4)
            ax.tick_params(length=0)
            if j == 0:
                ax.set_title(LABEL[marker], loc="left", pad=4)
            if i == 0:
                ax.text(1.0, 1.18, ["During stress (day 10)", "After 10 days of recovery (day 20)"][j],
                        transform=ax.transAxes, ha="right", va="bottom", color=INK2, fontsize=9.5, fontweight="bold")
    for ax in axes[-1]:
        pct_axis(ax, "x", [-75, -50, -25, 0, 50])
        ax.set_xlim(np.log(0.19), np.log(2.0))
    fig.supxlabel("Change relative to control beetles sampled on the same day (log scale)", y=0.075, fontsize=9, color=INK2)
    handles = legend_handles() + [
        Line2D([], [], marker="o", ls="", mfc=INK2, mec=INK2, ms=6, label="Filled: differs from control (Holm p < 0.05)"),
        Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=INK2, mew=1.4, ms=6, label="Open: not distinguishable from control")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8.5, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Thermal stress in Zygogramma bicolorata: effect sizes with 95% CI", x=0.02, y=0.99, ha="left", fontsize=12, fontweight="bold")
    fig.text(0.02, 0.953, "Re-analysis of group means digitised from Patel, Yadav & Kumar (2025), Can. Entomol. 157:e26. n = 10 beetles per group.",
             color=INK2, fontsize=8.5)
    fig.subplots_adjust(left=0.115, right=0.985, top=0.885, bottom=0.125, hspace=0.42, wspace=0.07)
    fig.savefig(FIGS / "effect_sizes.png", dpi=200)
    plt.close(fig)


def fig_defence_damage():
    dd = pd.read_csv(TABLES / "defence_damage.csv")
    fig, ax = plt.subplots(figsize=(7.4, 5.9))
    ax.axhline(0, color=INK2, lw=0.8)
    ax.axvline(0, color=INK2, lw=0.8)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    nudge = {("heat", "male"): (8, 8), ("heat", "female"): (8, 6), ("cold", "male"): (12, -19), ("cold", "female"): (-8, 10)}
    for (regime, sex), g in dd.groupby(["regime", "sex"]):
        s = g[g.phase == "stress"].iloc[0]
        r = g[g.phase == "recovery"].iloc[0]
        c = COLOUR[regime]
        for p in (s, r):
            ax.errorbar(p.defence_lnRR, p.damage_lnRR, xerr=p.defence_se, yerr=p.damage_se, fmt="none", ecolor=c, elinewidth=0.9, alpha=0.45)
        ax.annotate("", xy=(r.defence_lnRR, r.damage_lnRR), xytext=(s.defence_lnRR, s.damage_lnRR),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=1.8, shrinkA=5, shrinkB=5, mutation_scale=13))
        ax.plot(s.defence_lnRR, s.damage_lnRR, SHAPE[sex], ms=9, mfc=c, mec=SURFACE, mew=1.5, zorder=3)
        ax.plot(r.defence_lnRR, r.damage_lnRR, SHAPE[sex], ms=9, mfc=SURFACE, mec=c, mew=1.8, zorder=3)
        dx, dy = nudge[(regime, sex)]
        ax.annotate(f"{regime.capitalize()}, {sex}", (s.defence_lnRR, s.damage_lnRR), xytext=(dx, dy),
                    textcoords="offset points", ha="left" if dx > 0 else "right", fontsize=8.5, color=INK)
    pct_axis(ax, "x", [-60, -40, -20, 0, 20, 40])
    pct_axis(ax, "y", [0, 20, 40, 60, 80])
    ax.set_xlim(np.log(0.30), np.log(1.5))
    ax.set_ylim(np.log(0.93), np.log(1.95))
    ax.set_xlabel("Antioxidant defence: mean change in SOD and catalase activity vs same-day control")
    ax.set_ylabel("Oxidative damage: change in lipid peroxidation vs same-day control")
    kw = dict(color=INK2, fontsize=8.5, style="italic")
    ax.text(0.985, 0.975, "defence raised,\ndamage still high", transform=ax.transAxes, ha="right", va="top", **kw)
    ax.text(0.015, 0.975, "defence lowered,\ndamage raised", transform=ax.transAxes, ha="left", va="top", **kw)
    handles = legend_handles() + [
        Line2D([], [], marker="o", ls="", mfc=INK2, mec=INK2, ms=7, label="Filled: during stress (day 10)"),
        Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=INK2, mew=1.5, ms=7, label="Open: after recovery (day 20)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8.5, bbox_to_anchor=(0.5, 0.035))
    ax.set_title("Heat raises antioxidant defence during stress; ten days later it is below control in every group", loc="left", fontsize=9.5)
    fig.text(0.012, 0.012, "Arrows run from the stress phase to the recovery phase. Bars: ±1 SE. Data digitised from Patel et al. (2025).", color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(FIGS / "defence_damage_plane.png", dpi=200)
    plt.close(fig)


def fig_control_age():
    age = pd.read_csv(TABLES / "control_age_change.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.axvline(0, color=INK2, lw=0.8)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for i, marker in enumerate(MARKERS):
        for sex, off in (("male", 0.16), ("female", -0.16)):
            r = age[(age.marker == marker) & (age.sex == sex)].iloc[0]
            y = len(MARKERS) - 1 - i + off
            lo, hi = np.log(1 + r.pct_lo / 100), np.log(1 + r.pct_hi / 100)
            ax.plot([lo, hi], [y, y], color=COLOUR["control"], lw=2, solid_capstyle="round")
            ax.plot(r.lnRR, y, SHAPE[sex], ms=7, mfc=COLOUR["control"], mec=SURFACE, mew=1.2, zorder=3)
            ax.annotate(f"{r['pct_change']:+.0f}%", (hi, y), xytext=(9, 0), textcoords="offset points", va="center", fontsize=8, color=INK2)
    ax.set_yticks(range(len(MARKERS)), [LABEL[m] for m in MARKERS][::-1])
    ax.tick_params(length=0)
    pct_axis(ax, "x", [-50, -25, 0, 25, 50, 100, 150])
    ax.set_xlim(np.log(0.5), np.log(2.9))
    ax.set_xlabel("Control beetles: day 20 relative to day 10")
    ax.legend(handles=[Line2D([], [], marker=SHAPE[s], ls="", color=COLOUR["control"], ms=6, label=s.capitalize()) for s in ("male", "female")],
              loc="lower right", fontsize=8.5)
    ax.set_title("Unstressed control beetles change between the two sampling days", loc="left")
    fig.tight_layout()
    fig.savefig(FIGS / "control_age_change.png", dpi=200)
    plt.close(fig)


def fig_mortality():
    m = pd.read_csv(TABLES / "mortality_counts.csv")
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for i, regime in enumerate(["control", "cold", "heat"]):
        for sex, off in (("male", -0.13), ("female", 0.13)):
            r = m[(m.regime == regime) & (m.sex == sex)].iloc[0]
            ax.plot([i + off] * 2, [r.mortality_lo, r.mortality_hi], color=COLOUR[regime], lw=2, solid_capstyle="round")
            ax.plot(i + off, 100 * r.dead / r.n_start, SHAPE[sex], ms=8, mfc=COLOUR[regime], mec=SURFACE, mew=1.2, zorder=3)
            ax.annotate(f"{100 * r.dead / r.n_start:.0f}%", (i + off, r.mortality_hi), xytext=(0, 4), textcoords="offset points",
                        ha="center", fontsize=8, color=INK2)
    ax.set_xticks(range(3), ["Control\n(n = 150 per sex)", "Repeated cold\n(n = 500 per sex)", "Repeated heat\n(n = 500 per sex)"])
    ax.set_ylim(0, 58)
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylabel("Mortality over 20 days (%)")
    ax.tick_params(length=0)
    ax.legend(handles=[Line2D([], [], marker=SHAPE[s], ls="", color=INK2, ms=6, label=s.capitalize()) for s in ("male", "female")],
              loc="upper left", fontsize=8.5)
    ax.set_title("Mortality with 95% Wilson intervals, computed on beetle counts", loc="left")
    fig.tight_layout()
    fig.savefig(FIGS / "mortality.png", dpi=200)
    plt.close(fig)


def fig_anova_check():
    c = pd.read_csv(TABLES / "anova_vs_published.csv")
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    lim = (1, 2500)
    ax.plot(lim, lim, color=INK2, lw=0.8)
    ax.plot(c.F_published, c.F_raw_whiskers, "o", ms=6, mfc=SURFACE, mec="#8a8985", mew=1.3, label="SE read directly off the error bars")
    ax.plot(c.F_published, c.F_calibrated, "o", ms=6, mfc="#2a78d6", mec=SURFACE, mew=0.8, label="SE after the one-parameter calibration")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel("F published in Table 1")
    ax.set_ylabel("F recomputed from the digitised data")
    ax.text(900, 1150, "1 : 1", color=INK2, fontsize=8.5, rotation=45, ha="center", va="center")
    ax.legend(loc="lower right", fontsize=8.5)
    ax.set_title("The digitised data reproduce the published ANOVA", loc="left")
    fig.tight_layout()
    fig.savefig(FIGS / "anova_check.png", dpi=200)
    plt.close(fig)


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    for f in (fig_effect_sizes, fig_defence_damage, fig_control_age, fig_mortality, fig_anova_check):
        f()
        print("wrote", f.__name__)


if __name__ == "__main__":
    main()
