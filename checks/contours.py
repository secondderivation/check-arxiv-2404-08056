"""C6 and C7 contour tests with getdist's default 2D densities (the method of DESI VI's figures).

C6: in (w_p, wa) at DESI's pivot z_p of each SN sample, does the 95% region reach w_p = -1?
C7: does any grid cell of the 95% region in (w0, wa) lie outside zone Z4? Also a figure of the
pivot-shifted contours (the replot the paper could not make).

Run: CHAINS=<dir> python3 -I checks/contours.py <dr1|dr2> > checks/out/contours_<release>.out
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RUNS, OUT, load, cols, zone

rel = sys.argv[1]
ZP = {"PP": 0.26, "U3": 0.33, "DY5": 0.26, "CMB": 0.57}  # DESI VI's DR1 pivots
if rel != "dr1" and os.environ.get("ZP_FROM") != "dr1":  # M1: each release at its own DETF pivot from chain_checks.py, rounded like DESI's
    _S = json.load(open(os.path.join(os.path.dirname(__file__), "out", f"chain_summary_{rel}.json")))
    ZP = {sn: round(_S[f"{rel}_w0wa_{sn}"]["z_p"], 2) for sn in ZP if f"{rel}_w0wa_{sn}" in _S}
ok = lambda b: "PASS" if b else "FAIL"
samples = {}
for sn in ("PP", "U3", "DY5"):
    key = f"{rel}_w0wa_{sn}"
    if key not in RUNS: continue
    s = load(key); w0, wa = cols(s, "w0", "wa")
    om = 1 - 1 / (1 + ZP[sn])
    s.addDerived(w0 + wa * om, name="wp", label=r"w_{\rm p}")
    s.addDerived(w0, name="w0x", label="w_0")
    s.updateBaseStatistics()
    for (x, y, tag) in (("wp", "wa", "pivot"), ("w0x", "wa", "w0wa")):
        d = s.get2DDensityGridData(x, y, num_plot_contours=2)
        lv68, lv95 = sorted(d.contours, reverse=True)[0], sorted(d.contours, reverse=True)[1]
        X, Y = np.meshgrid(d.x, d.y)  # d.P has shape (len(y), len(x))
        in95 = d.P >= lv95
        if tag == "pivot":
            j = int(np.argmin(np.abs(d.x + 1)))
            onl = d.P[:, j] if abs(d.x[j] + 1) < 2 * (d.x[1] - d.x[0]) else np.array([0.0])
            reach95, reach68 = onl.max() >= lv95, onl.max() >= lv68
            xmin95 = X[in95].min(); xmax95 = X[in95].max()
            print(f"{key} (w_p at z_p={ZP[sn]}): 95% region spans w_p {xmin95:.3f} .. {xmax95:.3f}; the line w_p = -1 (grid col {d.x[j]:.4f}) reaches the 95% level: {ok(reach95)}, the 68% level: {reach68}")
            print(f"   max density on the line / 95% level = {onl.max() / lv95:.3f}")
        else:
            zn = zone(X[in95], Y[in95])
            outside = (zn != 4).sum()
            print(f"{key} (w0, wa): 95% region has {in95.sum()} grid cells, {outside} outside zone Z4: {ok(outside == 0)};"
                  f" max of w0 + wa over the 95% region {(X + Y)[in95].max():.3f}; min of w0 {X[in95].min():.3f}")
    samples[sn] = s

if samples:
    import matplotlib; matplotlib.use("Agg")
    from getdist import plots
    g = plots.get_single_plotter(width_inch=5)
    g.plot_2d(list(samples.values()), "wp", "wa", filled=True, colors=["#0173b2", "#ff9c34", "#33ad88"][: len(samples)])
    g.add_x_marker(-1); g.add_y_marker(0)
    g.add_legend([f"DESI {rel.upper()} + CMB + {k} (z_p = {ZP[k]})" for k in samples], legend_loc="lower left")
    g.export(os.path.join(OUT, f"pivot_contours_{rel}.png"))
    print(f"figure: checks/out/pivot_contours_{rel}.png")
