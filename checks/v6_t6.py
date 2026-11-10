"""V6: A's own check of critic C's post_hoc taper T6 (C's results/T6-DESIGN.md), with A's loader and
quadrature, and I(w0, wa) integrated numerically over a (C used the closed form).
T6(s) = exp[-I / (2 s^2)], I = int_0^1 max(0, -1 - w(a))^2 da, s = 0.10 (0.075, 0.15 as sensitivity).

Run: CHAINS=<dir> python3 -I checks/v6_t6.py > checks/out/v6_t6.out
"""
import json, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT, load, cols, gauss_level, in_box, BOX

A = (np.arange(4000) + 0.5) / 4000  # midpoint rule in a


def phantom_integral(w0, wa, chunk=20000):
    out = np.empty(len(w0))
    for i in range(0, len(w0), chunk):
        w = w0[i:i + chunk, None] + wa[i:i + chunk, None] * (1 - A[None, :])
        out[i:i + chunk] = (np.maximum(0, -1 - w) ** 2).mean(axis=1)
    return out


def t6(w0, wa, s):
    return np.exp(-phantom_integral(w0, wa) / (2 * s**2))


n = 1001  # E_pi0 by midpoint quadrature on DESI's box with the w0 + wa < 0 wedge (as prior_mean_quad)
g0 = np.linspace(*BOX["w0"], n + 1); g0 = (g0[1:] + g0[:-1]) / 2
ga = np.linspace(*BOX["wa"], n + 1); ga = (ga[1:] + ga[:-1]) / 2
W0, WA = np.meshgrid(g0, ga, indexing="ij"); m = in_box(W0, WA)
I_grid = phantom_integral(W0[m], WA[m])
V = json.load(open(os.path.join(os.path.dirname(__file__), "..", "paper_values.json")))["desi_vi_v3"]
lnB = dict(zip(("PP", "U3", "DY5"), [float(x) for x in re.findall(r"\d+(?:\.\d+)?", V["desi.lnB"]["value"])]))
out = {}
for s in (0.10, 0.075, 0.15):
    Epi = np.exp(-I_grid / (2 * s**2)).mean()
    print(f"\nT6 s = {s}: E_pi0[t] = {Epi:.6f}")
    for sn in ("PP", "U3", "DY5"):
        S = load(f"dr1_w0wa_{sn}"); w = S.weights; w0, wa = cols(S, "w0", "wa")
        t = t6(w0, wa, s); wt = w * t
        Ep = wt.sum() / w.sum(); neff = wt.sum() ** 2 / (wt ** 2).sum()
        _, _, n0 = gauss_level(w0, wa, w); _, _, nT = gauss_level(w0, wa, wt)
        d = np.log(Ep) - np.log(Epi)
        out[f"{sn}_s{s}"] = dict(E_prior=Epi, E_post=Ep, dlnB=d, lnB=lnB[sn] + d, n_eff=neff, n_sigma=nT, n_sigma_T0=n0)
        print(f"  {sn}: E_post {Ep:.5f}, Delta ln B {d:+.3f} -> ln B {lnB[sn] + d:+.3f}; LCDM level {nT:.3f} sigma (T0 {n0:.3f}); n_eff {neff:.0f}")
json.dump(out, open(os.path.join(OUT, "v6_t6.json"), "w"), indent=1, default=float)
