# Extension: how does *Z. bicolorata* compare with other insects?

Status: **not done yet.** The table and script are ready; the data extraction is manual work.

Question: is the size of the antioxidant response in *Z. bicolorata* (SOD +22 to +28 %, LPO
+54 to +74 % under repeated heat) typical for insects, or unusual?

How
1. For each paper below, find the figure or table giving SOD, CAT and MDA/LPO in a control
   group and a heat- or cold-stressed group.
2. Read off mean, SE (or SD) and n. For figures use WebPlotDigitizer (https://automeris.io).
3. Add one row per comparison to `data/literature/extracted_means.csv`. Fill in `source`
   (figure/table number) and `extracted_by`, and note anything uncertain in `notes`.
4. Run `python scripts/06_cross_species.py`. With three or more studies it draws a forest plot.

Candidate papers (all are cited in Patel et al. 2025; check what each one actually reports
before extracting)
- Chen et al. 2018, *Ophraella communa*, Frontiers in Physiology 9: 808 (open access; a leaf
  beetle and biocontrol agent, the closest comparison)
- Das, Malik & Barik 2018, *Galerucella placida*, Int. J. Trop. Insect Sci. 38: 400-409
- Jia et al. 2011, *Bactrocera dorsalis*, Florida Entomologist 94: 956-963
- Cui et al. 2011, *Chilo suppressalis*, J. Thermal Biology 36: 292-297
- Wang et al. 2022, *Sitophilus oryzae* and *Lasioderma serricorne*, J. Stored Prod. Res. 96: 101949
- Lalouette et al. 2011, *Alphitobius diaperinus* under fluctuating regimes, Comp. Biochem.
  Physiol. A 158: 229-234
- Lubawy et al. 2023, cockroach under fluctuating regimes, J. Exp. Biol. 226: jeb246287
- Farahani et al. 2020, *Ectomyelois ceratoniae*, PLOS One 15: e0228104 (open access)
- Yadav et al. 2024, *Z. bicolorata* under single short and long exposures, Physiol. Entomol.
  50: 139-151 (same species, different protocol)

Things to decide and write down before starting
- Which contrast to take when a paper has several temperatures or durations (for example: the
  highest sub-lethal temperature).
- How to handle different units and assay methods (lnRR is unit-free, which is the reason for
  using it).
- Whether whole-body and tissue-specific measurements can be pooled.
