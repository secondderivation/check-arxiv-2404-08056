"""C9 and C5, prior-only parts (no data): E_pi0[t] for the prior variants T1-T5 by quadrature and
Monte Carlo (control C0.3), and the marginal prior on w_p implied by DESI's prior (C5 'the induced
prior on w_p is still flat').

Run: python3 -I checks/prior_volume.py > checks/out/prior_volume.out
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import prior_mean_quad, prior_mean_mc, in_box, BOX

print("== C0.3 E_pi0[t] on the DESI box (w0 in [-3,1], wa in [-3,2], w0+wa<0) ==")
q0, area = prior_mean_quad("T0")
print(f"box area by quadrature: {area:.4f} (analytic 15.5)")
for T in ("T1", "T2", "T3", "T4", "T5", "CUT-2.9"):
    q, _ = prior_mean_quad(T); m = prior_mean_mc(T)
    rel = abs(q - m) / q
    print(f"  {T:<8} quadrature {q:.6f}  Monte Carlo {m:.6f}  rel. diff {rel:.2e} < 2e-3: {'PASS' if rel < 2e-3 else 'FAIL'}   ln E = {np.log(q):+.4f}")
print("  analytic CUT-2.9 fraction: 1 - (0.1 x 4)/15.5 =", f"{1 - 0.1 * 4 / 15.5:.6f}", "(the strip wa in [-3, -2.9] has the full w0 width 4, since w0 + wa < 0 holds there)")

print("\n== C5 marginal prior density of w_p implied by pi0 (the joint prior in (w_p, wa) is flat: Jacobian 1; the tilt below comes from the condition w0 + wa < 0, i.e. wa < -w_p/a_p) ==")
for zp in (0.26, 0.33):
    om = 1 - 1 / (1 + zp)
    wa = np.linspace(-3, 2, 200001); dwa = wa[1] - wa[0]
    for wp in (-1.2, -1.1, -1.0, -0.9, -0.8):
        w0 = wp - wa * om
        dens = in_box(w0, wa).sum() * dwa / 15.5
        print(f"  z_p={zp}: p(w_p={wp:+.1f}) = {dens:.5f}", end="")
        if wp == -1.2: ref = dens
        print(f"   relative to w_p=-1.2: {dens / ref:.4f}")
    vals = [in_box(w - wa * om, wa).sum() * dwa / 15.5 for w in np.linspace(-1.2, -0.8, 41)]
    var = (max(vals) - min(vals)) / np.mean(vals)
    print(f"  z_p={zp}: variation over w_p in [-1.2, -0.8]: {100 * var:.2f}% < 5%: {'PASS' if var < 0.05 else 'FAIL'}")
