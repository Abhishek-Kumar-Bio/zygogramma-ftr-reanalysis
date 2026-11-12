# Learning log

Answer each question in your own words before making the repository public. If you cannot
answer one, that is the part to study next.

## The biology
- [ ] What do SOD and catalase each do, and why are they measured together?
- [ ] What does the TBARS assay measure, and what are its known weaknesses?
- [ ] Why would antioxidant activity fall under cold but rise under heat?
- [ ] Enzyme activity went up under heat. Does that mean the genes were transcribed more?
      List three other ways activity per mg protein can change.

## The data
- [ ] Why are the SEs the weakest part of the digitised table? How many pixels long is a
      typical error bar?
- [ ] Explain in three sentences how the published F-ratios were used to calibrate the SEs.
- [ ] What does the calibration assume? What would break it?
- [ ] What can a re-analysis of 60 group means not tell you that 600 raw values could?

## The statistics
- [ ] What is a log response ratio, and why use it instead of a raw difference?
- [ ] Why compare with the same-day control rather than with the same group ten days earlier?
- [ ] What does the Holm adjustment do, and why are 40 unadjusted tests a problem?
- [ ] Why does a chi-square on percentages give a different answer from one on counts?
- [ ] What does the Mantel-Haenszel odds ratio of 1.40 mean in plain words?
- [ ] The regime x phase interaction is large for every marker. What does that say about
      interpreting the main effect of regime?

## The code
- [ ] Walk through `factorial_anova_from_summary` line by line. Where does each SS come from?
- [ ] Why do the tests compare against statsmodels on simulated data?
- [ ] Change one digitised mean by 10 % and rerun. Which results move, and by how much?
