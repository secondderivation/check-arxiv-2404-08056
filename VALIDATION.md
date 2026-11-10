# Release validation contract

Written 19 October 2026, at P1. These are engineering rules for replaying the package. They are not part of the preregistration and do not replace the scientific tolerances in `PROTOCOL.md`.

## Three separate outcomes

1. **Scientific reproduction.** This uses the frozen claim criteria in `PROTOCOL.md` with Addendum A1. The final rulings after review are the claim outcomes in the README. A replay PASS changes none of them.
2. **Numerical replay.** This re-runs what can be re-run from the package and compares the result with the saved outputs, under the rules below.
3. **Distribution integrity.**
   - Every distributed file is checked against `MANIFEST.txt` (SHA-256).
   - DESI's chain files are downloaded from DESI's server and must match their pinned hashes (`datasets.json`, `fetch_data.py`). They are 56 DR1 files (about 0.65 GB) and 90 DR2 files (about 0.6 GB).
   - This verifies consistency, not scientific correctness.

## What is replayed, and what is not

Each check runs with `--replay`. A check run that way sends what it writes besides stdout (`chain_summary_*.json`, `priors_*.json`, `v6_t6.json`, the contour figures) to `_runs/<time>/replay-out/` instead of `checks/out/`. The work copy therefore stays identical to the manifest. Each run's stdout is compared with the published output in the package root (`compare_results.py`). `test_validation.py` tests this end to end.

**Deterministic** (`algebra`, `prior_volume`). sympy identities and a seeded quadrature. The text and verdict words must match exactly, with every number within 1e-9.

**Chain statistics** (`chain_checks`, `contours`, `priors` for DR1; `v6_t6`; the same three for DR2). These are weighted moments, getdist densities and reweighted expectations of DESI's released chains.
- The words and verdicts must match.
- Every printed value must be within 1e-3 absolute or 1e-4 relative. The comparison uses the printed decimal values and includes equality at the tolerance boundary: |new - old| <= max(1e-3, 1e-4 * max(|old|, |new|)). These tolerances can allow more than one last printed digit, depending on the field's precision.
- Correlations that are zero by construction print as round-off (around 1e-16). Any value below 1e-12 counts as zero.
- On our machine a rerun was byte-identical.

The prior-variant run `v6_t6` and the three DR2 runs are marked slow and run only with `--all`. The `v6_t6` variant was defined at V6, after the first results (see the README, N7).

**Not replayed.** Their saved outputs are the record.
- `quoted`, `desi_versions`, `fig1_lines` and `m1_dr2_ref` read the TeX and figure sources of this paper and of DESI's papers, which the package does not redistribute. Each script header names the source files, and the arXiv e-prints are 2404.08056v3, 2404.03002v1–v3 and 2503.14738v3.
- `c01_identity` compares the chain files with the project's registry. In the package, `fetch_data.py` does the same with the pinned hashes in `datasets.json`.
- `chain_checks_dr1_burn030` and `contours_dr2_at_dr1_pivots` are variants run with the environment settings `BURN_ALL=0.3` and `ZP_FROM=dr1`. To rerun them by hand, run `chain_checks.py dr1` or `contours.py dr2` with that setting from `checks/`.

The suite is complete only with `python3 reproduce.py --all`. No CI step publishes, creates a release or records owner approval.
