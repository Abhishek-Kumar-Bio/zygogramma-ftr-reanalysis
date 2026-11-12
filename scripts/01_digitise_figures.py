#!/usr/bin/env python3
"""Step 1. Digitise group means and standard errors from the published bar charts.

Source: Patel, Yadav & Kumar (2025) Can. Entomol. 157:e26, Figs 2-4.
The paper reports group means only as bar charts, so the numbers are read off the
figure images by pixel measurement. The images themselves are NOT redistributed in
this repository: extract them from your own copy of the PDF with

    pdfimages -png -f 10 -l 12 paper.pdf figures_src/im
    # then rename: im-000.png -> fig2.png, im-001.png -> fig3.png, im-002.png -> fig4.png

How it works, per panel
  1. y-axis  = the grey vertical line; tick marks = short dashes left of it.
     A straight line is fitted through all ticks -> pixel row <-> data value.
  2. bars    = runs of saturated colour rising from the zero line.
     Height is measured with sub-pixel coverage at the top edge.
  3. SE      = length of the upper whisker above the bar top, minus half the cap.

Resolution limit: one pixel is 0.5-1 % of the axis range, and most whiskers are
only 1-6 px long. Means are therefore good to about 1 % of the axis range, but an
individual SE can be off by +-1 px. `se_px` and `unit_per_px` are written to the
output so every downstream step can see how coarse each SE is.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "figures_src"
WHITE = np.array([255.0, 255.0, 255.0])

# panel -> (figure file, crop box x0,y0,x1,y1, axis maximum, number of ticks, marker, sex, unit)
PANELS = {
    "SOD_M": ("fig2.png", (0, 0, 550, 365), 1.2, 4, "SOD", "male", "U/mg protein"),
    "SOD_F": ("fig2.png", (555, 0, 1119, 365), 1.2, 4, "SOD", "female", "U/mg protein"),
    "CAT_M": ("fig2.png", (0, 365, 550, 745), 200, 6, "CAT", "male", "U/mg protein"),
    "CAT_F": ("fig2.png", (555, 365, 1119, 745), 200, 6, "CAT", "female", "U/mg protein"),
    "LPO_M": ("fig2.png", (0, 745, 550, 1096), 0.16, 5, "LPO", "male", "nmol TBARS/mg protein"),
    "LPO_F": ("fig2.png", (555, 745, 1119, 1096), 0.16, 5, "LPO", "female", "nmol TBARS/mg protein"),
    "GLU_M": ("fig3.png", (0, 0, 550, 385), 0.30, 6, "glucose", "male", "mg/ml"),
    "GLU_F": ("fig3.png", (555, 0, 1117, 385), 0.30, 6, "glucose", "female", "mg/ml"),
    "TAG_M": ("fig3.png", (0, 385, 550, 767), 1.0, 6, "triglycerides", "male", "mg/ml"),
    "TAG_F": ("fig3.png", (555, 385, 1117, 767), 1.0, 6, "triglycerides", "female", "mg/ml"),
}
# left-to-right order of the six bars in every panel
BAR_ORDER = [("cold", "stress", 10), ("control", "stress", 10), ("heat", "stress", 10),
             ("cold", "recovery", 20), ("control", "recovery", 20), ("heat", "recovery", 20)]


def load(path):
    return np.asarray(Image.open(path).convert("RGB")).astype(float)


def find_axis(a):
    """Column index of the y-axis line: the column with most grey, unsaturated pixels."""
    grey = a.mean(axis=2)
    sat = a.max(axis=2) - a.min(axis=2)
    count = ((grey < 200) & (sat < 40)).sum(axis=0)
    return int(np.argmax(count[: int(a.shape[1] * 0.45)]))


def find_ticks(a, xa, n_ticks):
    """Rows of the y-axis tick marks (top tick first), from a line fitted through all of them."""
    block = a[:, xa - 5:xa - 1].mean(axis=2)
    darkness = np.clip((235 - block) / 235, 0, 1).sum(axis=1)
    rows = np.where(darkness > 0.35)[0]
    groups = []
    for r in rows:
        if groups and r - groups[-1][-1] <= 1:
            groups[-1].append(r)
        else:
            groups.append([r])
    cents = [float(np.sum(np.array(g) * darkness[g]) / np.sum(darkness[g])) for g in groups if len(g) <= 4]
    best = None
    for i in range(len(cents)):                      # candidate top tick
        for j in range(i + 1, len(cents)):           # candidate zero tick
            step = (cents[j] - cents[i]) / (n_ticks - 1)
            if step < 15:
                continue
            seq = [cents[i] + k * step for k in range(n_ticks)]
            if all(min(abs(s - c) for c in cents) < 1.6 for s in seq):
                cand = [min(cents, key=lambda c: abs(c - s)) for s in seq]
                if best is None or cand[-1] > best[-1]:
                    best = cand
    if best is None:
        raise RuntimeError("tick marks not found")
    k = np.arange(n_ticks)
    slope, intercept = np.polyfit(k, best, 1)
    return [float(intercept + slope * i) for i in k]


def digitise_panel(img, box, ymax, n_ticks, n_bars=6):
    x0, y0, x1, y1 = box
    a = img[y0:y1, x0:x1]
    xa = find_axis(a)
    ticks = find_ticks(a, xa, n_ticks)
    y_top, y_zero = ticks[0], ticks[-1]
    px_per_unit = (y_zero - y_top) / ymax
    sat = a.max(axis=2) - a.min(axis=2)
    base = int(round(y_zero))
    coloured = sat > 60
    dark = (a.max(axis=2) < 130) & (sat < 60)

    # columns that carry a bar: a coloured/dark run rising from the zero line
    height = np.zeros(a.shape[1], int)
    for x in range(xa + 3, a.shape[1]):
        y, n = base - 2, 0
        while y > 0 and (coloured[y, x] or dark[y, x]):
            y -= 1
            n += 1
        height[x] = n
    runs = []
    for x in range(len(height)):
        if height[x] > 4 and coloured[base - 3, x]:
            if runs and x - runs[-1][1] <= 1:
                runs[-1][1] = x
            else:
                runs.append([x, x])
    runs = [r for r in runs if r[1] - r[0] >= 15]
    if len(runs) != n_bars:
        raise RuntimeError(f"expected {n_bars} bars, found {len(runs)}")

    out = []
    for bx0, bx1 in runs:
        xc = (bx0 + bx1) // 2
        side = [x for x in range(bx0 + 3, bx1 - 2) if abs(x - xc) > 5]      # avoid the whisker
        bar_colour = np.median(a[base - 6:base - 2, side].reshape(-1, 3), axis=0)
        contrast = np.linalg.norm(bar_colour - WHITE)
        heights = []
        for x in side:                                                     # sub-pixel bar height
            y, h = base - 1, 0.0
            while y > 0:
                cover = np.clip(np.linalg.norm(a[y, x] - WHITE) / contrast, 0, 1)
                if cover < 0.15:
                    break
                h += 1.0 if cover >= 0.85 else cover
                y -= 1
            heights.append(h)
        h_px = float(np.median(heights))
        top = base - h_px
        # upper whisker: longest dark run just above the bar top, near the bar centre
        best = (0, None, None)
        for x in range(xc - 4, xc + 5):
            col_dark = (a[:, x].max(axis=1) < 140) & (sat[:, x] < 70)
            t = int(round(top))
            y = t - 1
            while y > 0 and not col_dark[y] and t - y < 3:
                y -= 1
            n = 0
            while y > 0 and col_dark[y]:
                y -= 1
                n += 1
            if n > best[0]:
                best = (n, x, y + 1)
        n_up, xw, tip = best
        if xw is None:
            se_px = float("nan")
        else:
            cap = 0
            y = tip
            while y < tip + 4 and (a[y, xw - 3:xw + 4].max(axis=1) < 150).sum() >= 5:
                cap += 1
                y += 1
            se_px = (top - tip) - cap / 2.0
        out.append({"mean": (y_zero - top) / px_per_unit, "se": se_px / px_per_unit,
                    "se_px": se_px, "unit_per_px": 1 / px_per_unit})
    return out


def digitise_mortality(img):
    """Fig 4: mortality (red) and survival (green) bars for 3 regimes x 2 sexes. No error bars.

    The bars have a dark outline, so a raw height reads ~0.5 % low. Mortality and survival
    must add up to 100 %, so mortality is taken as red / (red + green), which cancels that offset.
    """
    a = img
    xa = find_axis(a)
    ticks = find_ticks(a, xa, 6)
    y_zero = ticks[-1]
    base = int(round(y_zero))
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    masks = {"red": (r > 140) & (g < 90) & (b < 90), "green": (g > 110) & (r < 90) & (b < 110)}
    heights = {}
    for name, mask in masks.items():
        cols = np.where(mask[base - 8:base - 4].all(axis=0))[0]
        runs = []
        for x in cols:
            if runs and x - runs[-1][1] <= 1:
                runs[-1][1] = x
            else:
                runs.append([x, x])
        runs = [rr for rr in runs if rr[1] - rr[0] >= 10]
        if len(runs) != 6:
            raise RuntimeError(f"expected 6 {name} bars, found {len(runs)}")
        hs = []
        for bx0, bx1 in runs:
            tops = []
            for x in range(bx0 + 3, bx1 - 2):
                y = base - 5                       # start inside the bar, above the axis line
                while y > 0 and mask[y, x]:
                    y -= 1
                tops.append(y + 1)
            hs.append(y_zero - float(np.median(tops)))
        heights[name] = np.array(hs)
    return list(100 * heights["red"] / (heights["red"] + heights["green"]))


def main():
    rows = []
    cache = {}
    for panel, (fname, box, ymax, n_ticks, marker, sex, unit) in PANELS.items():
        img = cache.setdefault(fname, load(SRC / fname))
        bars = digitise_panel(img, box, ymax, n_ticks)
        for (regime, phase, day), bar in zip(BAR_ORDER, bars):
            rows.append({"marker": marker, "sex": sex, "regime": regime, "phase": phase, "day": day,
                         "mean": round(bar["mean"], 5), "se": round(bar["se"], 5), "n": 10, "unit": unit,
                         "se_px": round(bar["se_px"], 1), "unit_per_px": round(bar["unit_per_px"], 6),
                         "source": f"Patel et al. 2025, {'Fig. 2' if fname == 'fig2.png' else 'Fig. 3'}"})
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "data" / "digitised_means.csv", index=False)
    print(df.to_string(index=False))

    pct = digitise_mortality(load(SRC / "fig4.png"))
    groups = [("cold", "male", 500), ("cold", "female", 500), ("control", "male", 150),
              ("control", "female", 150), ("heat", "male", 500), ("heat", "female", 500)]
    mort = pd.DataFrame([{"regime": rg, "sex": sx, "n_start": n, "mortality_pct_digitised": round(p, 2),
                          "source": "Patel et al. 2025, Fig. 4"}
                         for (rg, sx, n), p in zip(groups, pct)])
    mort.to_csv(ROOT / "data" / "mortality.csv", index=False)
    print(mort.to_string(index=False))


if __name__ == "__main__":
    main()
