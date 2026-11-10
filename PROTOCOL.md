# PROTOCOL: arXiv:2404.08056v3 (M. Cortês and A. R. Liddle, "Interpreting DESI's evidence for evolving dark energy", JCAP 12 (2024) 007)

Committed 2026-10-18 before any computation on data (G2). Changes go in a dated addendum at the end, also before the run they affect. Anything not planned here is labelled `post_hoc` in the report.

Inputs: `paper.json` (V0 freeze and hashes), `claims/C1–C10-r1.json`, `conventions.json`, `paper_values.json` (the paper's printed numbers and the DESI 2024 VI numbers it quotes, parsed from the two TeX files by `checks/v1_transcribe.py`), `datasets.json` (pinned chain files and hashes).

**What kind of paper this is.** A three-page interpretation note. It runs no fit of its own: its numbers are DESI 2024 VI's (DESI VI below), read off DESI's Figure 6 and text, and its new content is a reading of the DESI posterior in pivot coordinates (the "PhantomX coincidence") and an argument that DESI's evidence depends on the prior. The authors say they could not do the pivot transformation because DESI's chains were not public (footnote, l.119). They are public now. This review therefore (i) checks every number the paper quotes against DESI VI, (ii) re-derives the pivot algebra, (iii) does the transformation the paper estimated by eye, on DESI's own DR1 chains, and (iv) runs the impact test the paper's argument implies: what changing the prior actually does to the evidence.

Roles: this session is **A** (Claude). Everything it produces is **V2 only** (G19). V4 needs a blind **B** run by Codex in a new conversation with only its packet open (G23): the paper, DESI VI, the chain links, the claims and the tolerance table of §3, nothing else. A's expectations are in `A-EXPECTATIONS.md`, committed with this file and sealed from B until B freezes its results. V5 needs a critic **C**, also a new Codex conversation.

Status of what A has seen before this commit: the paper's TeX (all of it, v3) and a word diff of the v1, v2 and v3 PDFs (v1 had a section "Other issues" on apparent inconsistencies in DESI VI's w_a errors, removed in v3; v3 adds the maximum-density statement); renderings of Figs 1 and 2; DESI VI v3 TeX: Section 5.2 (w0waCDM, including its pivot values, Δχ²_MAP, Bayes factors, prior remarks), Table 2 (priors), Table 3 (parameter table) and Appendix A (DESI+SDSS swap). No chain file has been opened and nothing has been computed. From training A knows that DESI DR2 (2025) reported 2.8σ to 4.2σ for the same model, that DESI DR1 chains were released publicly in 2025, and that several later papers discussed prior and parametrisation effects on the DESI w0wa result; none of this is used for a claim outcome.

## 1. Implementation (fixed now)

**Code.** Own Python (numpy, scipy, getdist 1.7.7 for loading Cobaya chains and 2D densities, sympy for the algebra), written from the paper's equations and `conventions.json`. No author code exists (G13).

**Chains.** DESI's public DR1 Cobaya chains for flat w0waCDM with DESI BAO + CMB and each SN sample (PantheonPlus, Union3, DESY5), DESI BAO + CMB without SN, and the matching wCDM runs, as listed in `datasets.json` with SHA-256 per file. Burn-in, weights and convergence as `conventions.json` `chains`. Each chain directory's `.updated.yaml` is read to confirm the likelihoods, model and priors before use (C0.1).

**Pivot, zones, crossing, prior variants, Bayes-factor shift, LCDM level.** Exactly as `conventions.json` (`pivot`, `zones`, `crossing`, `prior_variants`, `bayes_factor_shift`, `lcdm_level`). The five prior variants T1–T5 are fixed now and none is added after seeing results.

**Algebra.** `checks/algebra.py` (sympy): decorrelation condition, invariance of the w_a marginal, Jacobian of (w0, wa) → (w_p, wa), the Gaussian identity between σ(w_p) and the conditional (constant-w) error, the zone partition, d ln ρ_de/d ln a and the location of the maximum of ρ_de.

## 2. Controls (run first)

- **C0.1 Chain identity.** SHA-256 of every file equal to `datasets.json`; each run's `.updated.yaml` names DESI DR1 BAO, the CMB likelihoods of DESI VI and the stated SN likelihood, the w0waCDM (or wCDM) model, and the priors w0 ∈ [−3, 1], wa ∈ [−3, 2] plus w0 + wa < 0. Gelman–Rubin R − 1 over the chain files reported (< 0.05 expected; above that, the run is used with a note).
- **C0.2 Fidelity gate (G15).** From the chains, for the three +SN combinations and DESI+CMB: weighted mean of w0 within ± 0.005 and of wa within ± 0.02 of DESI VI Table 3; their standard deviations within 10 % of the printed symmetric errors (or the mean of the two asymmetric ones); z_p within ± 0.02, w_p within ± 0.005 and σ(w_p) within 10 % of DESI VI's pivot list (l.917–920). wCDM: mean w within ± 0.005, σ within 10 %. A combination that fails is reported, and its chain-based numbers are not used for claim outcomes (DESI's printed values are used instead, with a note).
- **C0.3 Reweighting numerics.** t ≡ 1 gives Δ ln B = 0 exactly; E_π0[t] by quadrature and by Monte Carlo agree within 0.2 % for T1–T5; a box cut that removes no posterior sample (wa > −2.9 for the +SN runs) gives Δ ln B = −ln(area fraction), compared with the analytic value.

## 3. Checks, tolerances and outcomes

A miss is first checked against our transcription and the conventions (rounding, which DESI version, mean vs MAP) before it is called a discrepancy (G5).

| Claim | Check | Tolerance for `held:reproduced` |
|---|---|---|
| C1 | DESI VI v1 and v3 text; σ from Δχ²_MAP with k = 2 | Exact match of 2.5, 3.5, 3.9 in DESI VI; our conversion of −8.7, −15.2, −18.1 equals the printed σ after half-up rounding to one decimal; DESI states the numbers come from a likelihood (MAP) ratio |
| C2 | DESI VI Table 2; chain `.updated.yaml`; max of w0 + wa over posterior samples | Exact match; "nowhere near": the 99.9th weighted percentile of w0 + wa below −0.3 in all three +SN runs |
| C3 | DESI VI; 1/(1 + z) half-up; DETF pivot from chains | Exact match of 0.26, 0.33 with DESI VI; a_p conversions exact; chain z_p within ± 0.02 of 0.26, 0.33, 0.26 (this is also gate C0.2) |
| C4 | Fig. 2 parameters vs DESI VI; w_p from them and from the chains (mean, sd, at the MAP sample) | Fig. 2 parameters equal DESI VI numbers (which ones: marginal means or MAP, reported). "± 0.02": \|E[w_p] + 1\| ≤ 0.025 for each SN sample. "a few hundredths": \|E[w_p] + 1\| ≤ 0.06 for each. The blobs equal w0 + wa(1 − a_p) of the drawn lines within 0.005 |
| C5 | sympy; chains | Algebra: each property holds exactly (or, for the Gaussian statement, as an identity of the Gaussian). Numbers: \|corr(w_p, wa)\| < 0.01 at the DETF pivot; sd(wa) unchanged; σ(w_p)/σ(w_const) within [0.8, 1.25] for each SN sample, with σ(w_const) from DESI VI's wCDM rows (and the wCDM chains, reported). "Induced prior on w_p is still flat": the marginal prior density of w_p implied by π0 varies by < 5 % over w_p ∈ [−1.2, −0.8] |
| C6 | Chains in (w0, wa) and (w_p, wa) | "Major axes nearly coincide with the w_p = −1 lines": angle between the principal axis and the line ≤ 5°, and the offset d = (E[w_p] + 1)/σ(w_p) with \|d\| ≤ 1, at the DESI pivot of each SN sample. \|d\| in (1, 2]: partially, S0 note; \|d\| > 2: the wording does not hold for that sample (S1 if the conclusion relies on it). "Accurately vertical": \|corr(w_p, wa)\| < 0.05. "All three 95 % contours overlap w_p = −1": the getdist 2D density in (w_p, wa) on the line w_p = −1 reaches the 95 % contour level for each sample |
| C7 | sympy; Fig. 1 rendering at 300 dpi; chains | Zones as stated; Fig. 1's red, magenta, blue and orange lines within 0.02 in w0 of their loci (w0 = −1; w0 + wa = −1; w_p = −1 at z 0.26 and 0.33) at three heights each; posterior mass outside Z4 < 5 % and no cell of the getdist 95 % region outside Z4, for each SN sample |
| C8 | sympy; chains | Derivation holds; P(wa < 0) ≥ 0.95 for each SN sample (so the crossing is a maximum of ρ_de); DESI's z_p lies inside the 68 % equal-tail interval of the posterior of z_x for each |
| C9 | (a) DESI VI text; (b) reweighting T1–T5 | (a) holds as written if DESI's significances depend on the shape of the prior or DESI is silent on prior dependence; it is narrowed (S0–S1, G20) if the significances are likelihood ratios that do not depend on the prior inside the box and DESI states the prior dependence of its contours and evidence. (b) "Downweighting ... would reduce the discrepancy": for T1–T4, held if both the Bayes factor (ln B_T < ln B_DESI) and the Gaussian LCDM level (n σ_T < n σ_0) decrease for all three SN samples; partially if only one of the two measures decreases (or only for some samples); not held if neither. "Diminish or perhaps disappear": disappear means ln B_T ≤ 0 or n σ_T < 2, reported per sample. The likelihood-ratio statistic is reported as unchanged whenever the MAP lies inside the prior's support with t > 0 |
| C10 | DESI VI v1 and v3 | Constant w: \|w + 1\| < 2σ (σ on the side towards −1) in every wCDM row of DESI VI Table 3, else partially with the row named. SDSS swap: each reduction (2.5 → 2.2 etc.) within [0.25, 0.75]σ. z = 0.51: DESI VI names it as the single most anomalous point; wording difference ("significantly" vs DESI's "mildly", ~2σ) is S0. Curvature: \|Ω_k\|/σ < 0.5 in the three kw0wa rows and w0, wa within 0.5σ of the flat values |

**Severity, fixed now (G18, 05 §3).** Quoted-number mismatches (C1, C2, C3, C10) are S0 unless a conclusion of this paper rests on the number, then S1. C4/C6: a PhantomX statement that fails for one SN sample is S0 if the paper's other wording ("a few hundredths") still holds, S1 if it fails for the sample DESI calls the strongest evidence (DESY5) and the paper's conclusion rests on it. C9 is the paper's main conclusion: if the impact run shows that down-weighting deep-phantom or far-from-LCDM models increases the Bayes factor for w0waCDM, the "would diminish" sentence is contradicted for that measure; this is an S2 candidate only if it reverses the abstract's conclusion, otherwise S1 (narrowing, G20). The paper is interpretive: "the prior is subjective" is an opinion and is not graded.

## 4. Data re-tests (only after the C0.2 fidelity gate passes)

Data traces, not verification; they never change a claim outcome (G17).
- **M1 DESI DR2** (update role): the DESI DR2 Results II public chains, DESI DR2 BAO + CMB + PantheonPlus, Union3, DESY5 (and DES-Dovekie if a public run exists): z_p, w_p and σ(w_p), d and the axis angle (C6), zone masses (C7), z_x (C8), and T1–T5 Δ ln B and n σ (C9), with the same code. If DESI DR2 prints pivot values or Bayes factors, they are compared as reference values.

## 5. Literature (G21), after the runs

INSPIRE citing list of 2404.08056 (saved on the owner's Mac): papers that redo the pivot analysis, test prior dependence of the DESI w0wa result (profile likelihoods, alternative priors, Bayesian evidence), the Shlivko & Steinhardt exchange the note responds to, and DESI DR2's own treatment of priors. Also the v1 → v2 → v3 changes of this paper and of DESI VI (v1 PDFs and sources kept). Results to `literature.md` with queries and dates.

## 6. Stop rules

Round cap 2 per claim. Claim cap 10, 10 used. Every number in the report is printed by a script in `checks/` into `checks/out/` (G14). Effort is measured from the session log at the end of each stage (G28). Compute cap for the A pass: 2 h of wall clock for all chain analyses and re-tests; what has not finished by then is `not_run`.

## Addendum A1 (2026-10-19, before any chain number is computed)

Written after the chain files arrived and were hashed and their `.yaml` files read, but before any chain column was loaded or any statistic computed.

1. **Files and pins.** The DR1 and DR2 runs are pinned file by file (URL and SHA-256) in `datasets/registry.json` as `desi-dr1-cosmo-chains` and `desi-dr2-cosmo-chains`; `checks/chains.json` maps each run key to its file root. Every local file equals DESI's official checksum list (re-hashed after the transfer from the owner's Mac).
2. **Burn-in.** DESI's DR1 runs with CMB are published as importance-sampled files (`chain.post.importance.*`) whose `.input.yaml` says `skip: 0.3, thin: 4`, i.e. DESI already removed 30 % burn-in. For these files the burn-in is 0; for the raw chains (DR1 DESI-only and every DR2 run) it is 0.3 as in `conventions.json`. A variant with 0.3 removed everywhere (`BURN_ALL=0.3`) is reported for C0.2 only; it does not decide anything.
3. **The w0 + wa < 0 condition.** No run's `.updated.yaml` lists it as a prior; the sampled box is w ∈ [−3, 1], wa ∈ [−3, 2]. C0.1 is read as "box priors as in DESI VI Table 2", and C2's "nowhere near" test (99.9th percentile of w0 + wa) is unchanged. `conventions.json` `priors_reference` keeps the wedge because DESI VI Table 2 states it; the prior volume E_π0 is computed with it as already committed.
4. **Likelihoods.** DR1 CMB is Planck 2018 plik TTTEEE + low-ℓ TT, EE with Planck/ACT DR6 lensing added by importance sampling (as DESI VI says); DR2 uses Planck NPIPE CamSpec. C0.1 checks these names.
5. **M1 re-test scope.** DESI's server has no DES-Dovekie run for DR2, so M1 has no Dovekie part. DR2 runs are used with the same code and the same gate quantities, compared with DESI DR2 Results II (2503.14738) where it prints them (`desi-dr2-cosmology`, reference value). DESI's own `.margestats` files in the DR2 directories are a cross-check of our loading, reported, not used for outcomes.
