"""C9 impact run: importance-reweight DESI's w0waCDM chains with the prior variants T1-T5
(conventions.json prior_variants, bayes_factor_shift, lcdm_level) and control C0.3.

Run: CHAINS=<dir> python3 -I checks/priors.py <dr1|dr2> > checks/out/priors_<release>.out
"""
import json, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RUNS, OUT, load, cols, wstats, taper, prior_mean_quad, gauss_level

rel = sys.argv[1]
V = json.load(open(os.path.join(os.path.dirname(__file__), "..", "paper_values.json")))["desi_vi_v3"]
lnB_desi = dict(zip(("PP", "U3", "DY5"), [float(x) for x in re.findall(r"\d+(?:\.\d+)?", V["desi.lnB"]["value"])])) if rel == "dr1" else {}
EPI = {T: prior_mean_quad(T)[0] for T in ("T0", "T1", "T2", "T3", "T4", "T5", "CUT-2.9")}
out = {}
for sn in ("PP", "U3", "DY5", "CMB"):
    key = f"{rel}_w0wa_{sn}"
    if key not in RUNS: continue
    s = load(key); w = s.weights; w0, wa = cols(s, "w0", "wa")
    # per chain-file split for the spread of E_post[t]
    chain_id = np.concatenate([np.full(len(c.weights), j) for j, c in enumerate(s.getSeparateChains())]) if hasattr(s, "getSeparateChains") else np.zeros(len(w))
    d2, lvl0, n0 = gauss_level(w0, wa, w)
    m_a, s_a = wstats(wa, w)
    print(f"\n=== {key} ===  baseline (DESI prior): Gaussian LCDM level {lvl0:.6f} = {n0:.2f} sigma (2D); wa = {m_a:.3f} +- {s_a:.3f}, wa=0 at {m_a / s_a:+.2f} sd"
          + (f"; DESI ln B = {lnB_desi[sn]}" if sn in lnB_desi else ""))
    out[key] = {"baseline": dict(n_sigma=n0, wa=m_a, wa_sd=s_a)}
    for T in ("T0", "CUT-2.9", "T1", "T2", "T3", "T4", "T5"):
        t = taper(T, w0, wa)
        Ep = (w * t).sum() / w.sum()
        per = [((w * t)[chain_id == j].sum() / w[chain_id == j].sum()) for j in np.unique(chain_id)]
        sp = np.std(np.log(per), ddof=1) / np.sqrt(len(per)) if len(per) > 1 and min(per) > 0 else float("nan")  # nan: a chain file has no sample with t > 0
        neff = (w * t).sum() ** 2 / ((w * t) ** 2).sum() if Ep > 0 else 0
        if Ep > 0:
            dlnB = np.log(Ep) - np.log(EPI[T])
            wt = w * t
            ma, sa = wstats(wa, wt)
            if neff >= 50:
                _, lvl, n = gauss_level(w0, wa, wt)
            else:  # too few effective samples for a covariance; reported as not computed
                n = float("nan")
            flag = "" if neff >= 50 else "  [n_eff < 50: E_post[t] rests on a handful of samples, a rough estimate only]"
            line = (f"  {T:<8} E_post[t] = {Ep:.5f} (n_eff {neff:.0f}), E_pi0[t] = {EPI[T]:.5f}: Delta ln B = {dlnB:+.3f} (+- {sp:.3f} over chain files)"
                    + (f" -> ln B = {lnB_desi[sn] + dlnB:+.2f}" if sn in lnB_desi else "")
                    + f"; LCDM level {n:.2f} sigma (2D, Gaussian; change {n - n0:+.2f}); wa = {ma:.3f} +- {sa:.3f} ({ma / sa:+.2f} sd)" + flag)
            out[key][T] = dict(E_post=Ep, E_prior=EPI[T], dlnB=dlnB, dlnB_err=sp, n_eff=neff, n_sigma=n, wa=ma, wa_sd=sa)
        else:
            ub = np.log(w.max() / w.sum()) - np.log(EPI[T])
            line = f"  {T:<8} no sample with t > 0: Delta ln B < {ub:+.2f} (one-sample bound)"
            out[key][T] = dict(E_post=0.0, E_prior=EPI[T], dlnB_upper=ub)
        if T == "CUT-2.9":
            line += f"   [control: analytic -ln(area fraction) = {-np.log(EPI[T]):+.4f} if no posterior mass below wa = -2.9; mass there {1 - Ep:.2e}]"
        print(line)
json.dump(out, open(os.path.join(OUT, f"priors_{rel}.json"), "w"), indent=1, default=float)
