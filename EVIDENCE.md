# Saved evidence — scoped public edition

The pivot transformation and main interpretation reproduce, with qualifications to the precision and direction of prior effects.

DESY5 does not meet the stated pivot precision; prior effects depend on the chosen prior.

## C1

2.5/3.5/3.9 sigma printed in DESI VI v1-v3; Delta chi2_MAP -8.7/-15.2/-18.1 with k=2 give 2.486/3.481/3.852 (A and B)

Outcome: held:reproduced. Evidence adequacy: pass. This is a retrospective record of existing checks.

## C2

Priors match DESI Table 2; w0+wa 99.9th percentiles -0.98/-1.17/-1.14, below -0.3; no explicit wedge entry in the yaml (Addendum A1)

Outcome: held:reproduced. Evidence adequacy: pass. This is a retrospective record of existing checks.

## C3

z_p 0.26/0.33/0.26 and a_p 0.79/0.75 exact; chain DETF pivots 0.2615/0.3260/0.2628; fidelity gate passed for all four runs

Outcome: held:reproduced. Evidence adequacy: pass. This is a retrospective record of existing checks.

## C4

Fig. 2 values are DESI's marginal means, not best fits (N1); '+-0.02' holds for PP only, 'a few hundredths' for all three (|E[w_p]+1| 0.018/0.040/0.055); DY5 carried by N2 at C6

Outcome: held:discrepancy. Evidence adequacy: pass. This is a retrospective record of existing checks.

## C5

Algebra holds (sympy, corr ~1e-16, sigma ratios in [0.8, 1.25]); marginal prior on w_p tilts 11.8-12.3% over [-1.2,-0.8], not flat; effect on E[w_p] ~2e-4 (N3, S1 after V6)

Outcome: held:discrepancy. Evidence adequacy: pass. This is a retrospective record of existing checks.

## C6

Axes within 0.18 deg of w_p=-1; all 95% contours reach the line; d at DESI pivots 0.66/0.97/2.15: DY5 does not lie on the line (N2)

Outcome: held:discrepancy. Evidence adequacy: pass. This is a retrospective record of existing checks.

## C7

Zones hold; mass outside Z4 <= 0.003, no 95% cell outside Z4; Fig. 1 magenta line off by up to 0.068 in w0 (N4)

Outcome: held:discrepancy. Evidence adequacy: pass. This is a retrospective record of existing checks.

## C8

Maximum-density derivation holds, P(wa<0) >= 0.9989; DY5 pivot 0.26 outside z_x 68% [0.307, 0.418] (N5, same fact as N2); partially_reproduced

Outcome: held:discrepancy. Evidence adequacy: pass. This is a retrospective record of existing checks.

## C9

(a) DESI stated prior dependence of its evidence from v1 (N6, S1 after V6); (b) T1-T4 lower the LCDM level, ln B falls for all three only under T2, rises under T3/T4; post_hoc T6 lowers both (N7); partially_reproduced

Outcome: held:discrepancy. Evidence adequacy: pass. This is a retrospective record of existing checks.

## C10

wCDM rows within 2 sigma of -1, SDSS swap 0.3/0.4/0.4 sigma, curvature |Omega_k|/sigma <= 0.21; 'significantly' for DESI's 'mildly' at z=0.51 (N6 wording list)

Outcome: held:discrepancy. Evidence adequacy: pass. This is a retrospective record of existing checks.

The linked outputs, protocol and byte-identical original archive are included. The recorded 10 November P2 traced the applicable original numerical replay. Subsequent editorial validation checks the revised package against that frozen evidence; it does not repeat the scientific calculation.

## Historical unavailable diagnostics

The following retained historical fields contain nonfinite values and are unavailable diagnostics, not finite validated estimates or zero. The original JSON remains unchanged. This is a disclosure of old evidence, not a fresh calculation failure.

- [checks/out/priors_dr2.json](checks/out/priors_dr2.json): `/dr2_w0wa_PP/T5/n_sigma`, `/dr2_w0wa_U3/T5/dlnB_err`, `/dr2_w0wa_U3/T5/n_sigma`, `/dr2_w0wa_CMB/T5/dlnB_err`, `/dr2_w0wa_CMB/T5/n_sigma`.
- [checks/out/priors_dr1.json](checks/out/priors_dr1.json): `/dr1_w0wa_PP/T5/n_sigma`, `/dr1_w0wa_U3/T5/dlnB_err`, `/dr1_w0wa_U3/T5/n_sigma`, `/dr1_w0wa_DY5/T5/dlnB_err`, `/dr1_w0wa_DY5/T5/n_sigma`, `/dr1_w0wa_CMB/T5/dlnB_err`, `/dr1_w0wa_CMB/T5/n_sigma`.
