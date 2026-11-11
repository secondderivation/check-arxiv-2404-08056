# Interpreting DESI's evidence for evolving dark energy

[2404.08056v3](https://arxiv.org/abs/2404.08056v3) · Review 0.2.0-rc3

**Reproduced in scope.** The pivot transformation and main interpretation reproduce, with qualifications to the precision and direction of prior effects.

**Limit:** DESY5 does not meet the stated pivot precision; prior effects depend on the chosen prior.

## Checked scope

- The DESI significances, priors and pivot redshifts the paper quotes, against DESI 2024 VI v1, v2 and v3 (C1–C3, C10)
- The pivot values of w(a) behind Fig. 2 and the ± 0.02 / 'a few hundredths' statements, from DESI's chains (C4)
- The algebra of the pivot transformation, and the claim that the induced prior on w\_p is flat (C5)
- The orientation of the w0–w\_a contours, their offset from w\_p = −1, and the pivot-shifted contours the paper could not draw (C6)
- The four zones of Fig. 1, where the posterior lies among them, and the drawn lines (C7)
- The phantom crossing, the maximum of the dark-energy density, and its position against the pivot (C8)
- The prior argument: what DESI's significances and Bayes factors depend on, and what five prior variants fixed before the run do to the Bayes factor and to the credible level of ΛCDM (C9)

Not checked: A new fit or optimisation. All numbers come from reweighting DESI's released chains, so DESI's own maximum-posterior statistic was not recomputed under a changed prior; New nested-sampling evidences. Absolute Bayes factors are DESI's printed ln B (0.65, 2.4, 2.8) plus the change computed here. If DESI normalised its prior over the full box rather than the w0 + w\_a \< 0 wedge, every absolute value would move by 0.255 in magnitude; the changes would not; The never-phantom prior (T5) beyond a rough estimate: only 2 to 31 effective samples survive it, so it needs its own sampler; How wide the observational window behind the pivot is, and how likely an alignment of the contours with w\_p = −1 would be under ΛCDM; Whether the CPL family describes the dark energy beyond the observed range; the 'maximum it will ever have' statement holds only for CPL extrapolated to all times

## Review record

| Claim | Method | Outcome / adequacy | Evidence |
|---|---|---|---|
| C1: C1 in fixed reviewed version — 2.5/3.5/3.9 sigma printed in DESI VI v1-v3; Delta chi2\_MAP -8.7/-15.2/-18.1 with k=2 give 2.486/3.481/3.852 (A and B) | source\_check | held:reproduced / pass | [checks/out/quoted.out](checks/out/quoted.out) |
| C2: C2 in fixed reviewed version — Priors match DESI Table 2; w0+wa 99.9th percentiles -0.98/-1.17/-1.14, below -0.3; no explicit wedge entry in the yaml (Addendum A1) | source\_check | held:reproduced / pass | [checks/out/chain_checks_dr1.out](checks/out/chain_checks_dr1.out) |
| C3: C3 in fixed reviewed version — z\_p 0.26/0.33/0.26 and a\_p 0.79/0.75 exact; chain DETF pivots 0.2615/0.3260/0.2628; fidelity gate passed for all four runs | source\_check | held:reproduced / pass | [checks/out/chain_checks_dr1.out](checks/out/chain_checks_dr1.out) |
| C4: C4 in fixed reviewed version — Fig. 2 values are DESI's marginal means, not best fits (N1); '+-0.02' holds for PP only, 'a few hundredths' for all three (\|E\[w\_p\]+1\| 0.018/0.040/0.055); DY5 carried by N2 at C6 | source\_check | held:discrepancy / pass | [checks/out/chain_checks_dr1.out](checks/out/chain_checks_dr1.out) |
| C5: C5 in fixed reviewed version — Algebra holds (sympy, corr ~1e-16, sigma ratios in \[0.8, 1.25\]); marginal prior on w\_p tilts 11.8-12.3% over \[-1.2,-0.8\], not flat; effect on E\[w\_p\] ~2e-4 (N3, S1 after V6) | source\_check | held:discrepancy / pass | [checks/out/prior_volume.out](checks/out/prior_volume.out) |
| C6: C6 in fixed reviewed version — Axes within 0.18 deg of w\_p=-1; all 95% contours reach the line; d at DESI pivots 0.66/0.97/2.15: DY5 does not lie on the line (N2) | source\_check | held:discrepancy / pass | [checks/out/contours_dr1.out](checks/out/contours_dr1.out) |
| C7: C7 in fixed reviewed version — Zones hold; mass outside Z4 \<= 0.003, no 95% cell outside Z4; Fig. 1 magenta line off by up to 0.068 in w0 (N4) | source\_check | held:discrepancy / pass | [checks/out/contours_dr1.out](checks/out/contours_dr1.out) |
| C8: C8 in fixed reviewed version — Maximum-density derivation holds, P(wa\<0) \>= 0.9989; DY5 pivot 0.26 outside z\_x 68% \[0.307, 0.418\] (N5, same fact as N2); partially\_reproduced | source\_check | held:discrepancy / pass | [checks/out/chain_checks_dr1.out](checks/out/chain_checks_dr1.out) |
| C9: C9 in fixed reviewed version — (a) DESI stated prior dependence of its evidence from v1 (N6, S1 after V6); (b) T1-T4 lower the LCDM level, ln B falls for all three only under T2, rises under T3/T4; post\_hoc T6 lowers both (N7); partially\_reproduced | source\_check | held:discrepancy / pass | [checks/out/priors_dr1.out](checks/out/priors_dr1.out) |
| C10: C10 in fixed reviewed version — wCDM rows within 2 sigma of -1, SDSS swap 0.3/0.4/0.4 sigma, curvature \|Omega\_k\|/sigma \<= 0.21; 'significantly' for DESI's 'mildly' at z=0.51 (N6 wording list) | source\_check | held:discrepancy / pass | [checks/out/quoted.out](checks/out/quoted.out) |

[Protocol and provenance](PROTOCOL.md). This review record is retrospective.

## Review provenance

Prepared with AI assistance. Historical Claude first pass and Codex independent review/critique remain distinct; this assembly reads saved evidence and is neither blind B nor new P2. No scientific code was run in this assembly. Before the replay command below, extract legacy-original.zip into a new directory and read its README.md and requirements.txt for installation, dependency, data and scope instructions.

Tools: Claude (historical A); Codex (historical B/C and current informed assembly). Human involvement: Owner selected the preparation scope and approved the predecessor release published on 10 November 2026; the owner approved publication of this editorial revision on 11 November 2026.

## Replay

~~~sh
mkdir replay-work && python3 -m zipfile -e legacy-original.zip replay-work && cd replay-work && python3 reproduce.py --all
~~~

See original source-payload README and recorded P2; not rerun during assembly.. Unpacks the unchanged original candidate into a new work directory and runs its documented scoped replay. Package validation checked the artifact and inherited the earlier replay after verifying the relevant files were unchanged; it did not rerun the scientific calculation.
