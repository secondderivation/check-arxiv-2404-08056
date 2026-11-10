"""C0.1-C0.2 and the chain parts of C2-C8 on DESI's public chains.

Run: CHAINS=<dir> python3 -I checks/chain_checks.py <release: dr1|dr2> > checks/out/chain_checks_<release>.out
Run keys and roots: checks/chains.json. Conventions: conventions.json (pivot, zones, crossing, chains).
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RUNS, OUT, CHAINS, load, cols, wstats, wcov, wquant, pivot, zone, crossing_z, level_to_sigma

rel = sys.argv[1]
V = json.load(open(os.path.join(os.path.dirname(__file__), "..", "paper_values.json")))["desi_vi_v3"]
ok = lambda b: "PASS" if b else "FAIL"
SN = ("PP", "U3", "DY5")
summary = {}

def fnum(s):
    import re
    return [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d+)?", s)]

def fields(s):
    """'Om 0.3085 +- 0.0068; w0 -0.827 +- 0.063' -> {'Om': [0.3085, 0.0068], 'w0': [...]} (label = first token)."""
    out = {}
    for part in s.split(";"):
        lab, rest = part.strip().split(" ", 1)
        out[lab] = fnum(rest)
    return out

def sym(v):
    """Mean and symmetric error from [mean, err] or [mean, +up, -down]."""
    return v[0], (v[1] if len(v) == 2 else (v[1] + abs(v[2])) / 2)

for key in [k for k in RUNS if k.startswith(rel + "_")]:
    r = RUNS[key]
    s = load(key)
    w = s.weights
    print(f"\n=== {key}: {r['root']} ===")
    print(f"  samples after burn-in ({100 * float(os.environ.get('BURN_ALL', r.get('burn', 0.3))):.0f}% removed here, Addendum A1): {len(w)}, effective {w.sum() ** 2 / (w ** 2).sum():.0f}")
    try:
        print(f"  Gelman-Rubin R-1 (getdist, over chain files): {s.getGelmanRubin():.4f}")
    except Exception as e:
        print(f"  Gelman-Rubin: not available ({type(e).__name__}: {e})")
    names = [p.name for p in s.paramNames.names]
    if "wa" not in names:
        if "w" in names:
            m, sd = wstats(np.asarray(s.getParams().w), w)
            print(f"  wCDM: w = {m:.4f} +- {sd:.4f}")
            summary[key] = dict(w_mean=m, w_sd=sd)
        continue
    w0, wa = cols(s, "w0", "wa")
    m0, s0 = wstats(w0, w); ma, sa = wstats(wa, w)
    om = np.asarray(s.getParams().omegam) if "omegam" in names else None
    pv = pivot(w0, wa, w)
    wp = w0 + wa * (1 - pv["a_p"])
    rho = wcov(wp, wa, w) / np.sqrt(wcov(wp, wp, w) * wcov(wa, wa, w))
    print(f"  w0 = {m0:.4f} +- {s0:.4f}; wa = {ma:.4f} +- {sa:.4f} (68% equal-tail {wquant(wa, w, 0.16):.3f} .. {wquant(wa, w, 0.84):.3f})" + (f"; Om = {wstats(om, w)[0]:.4f} +- {wstats(om, w)[1]:.4f}" if om is not None else ""))
    print(f"  DETF pivot: a_p = {pv['a_p']:.4f}, z_p = {pv['z_p']:.4f}, w_p = {pv['wp_mean']:.4f} +- {pv['wp_sd']:.4f} (formula sd {pv['wp_sd_formula']:.4f}); corr(w_p, wa) = {rho:+.2e}; sd(wa) unchanged by construction")
    # MAP-like sample
    lp = np.asarray(s.loglikes)  # getdist loglikes = minuslogpost
    i = int(np.argmin(lp))
    print(f"  highest-posterior sample: w0 = {w0[i]:.4f}, wa = {wa[i]:.4f}, w_p(DETF) = {wp[i]:.4f}, w at the paper's a_p 0.79/0.75 = {w0[i] + wa[i] * 0.21:.4f}/{w0[i] + wa[i] * 0.25:.4f}")
    # orientation (C6)
    C = np.array([[wcov(w0, w0, w), wcov(w0, wa, w)], [wcov(w0, wa, w), wcov(wa, wa, w)]])
    ev, evec = np.linalg.eigh(C); v = evec[:, -1]
    ang_axis = np.degrees(np.arctan2(v[1], v[0])) % 180
    d = (pv["wp_mean"] + 1) / pv["wp_sd"]
    print(f"  principal axis angle {ang_axis:.2f} deg; axis ratio {np.sqrt(ev[-1] / ev[0]):.1f}; offset d = (E[w_p]+1)/sd(w_p) = {d:+.3f}")
    for zp in (0.26, 0.33):
        om_ = 1 - 1 / (1 + zp)
        ang_line = np.degrees(np.arctan2(1, -om_)) % 180
        wpz = w0 + wa * om_
        mz, sz = wstats(wpz, w)
        print(f"    line w(z={zp}) = -1: angle {ang_line:.2f} deg, |axis - line| = {abs(ang_axis - ang_line):.2f} deg; E[w(z={zp})] = {mz:.4f} +- {sz:.4f}, offset {(mz + 1) / sz:+.3f} sd")
    # zones (C7)
    zn = zone(w0, wa)
    mass = {k: w[zn == k].sum() / w.sum() for k in (1, 2, 3, 4)}
    print(f"  zone masses: Z1 {mass[1]:.5f}  Z2 {mass[2]:.5f}  Z3 {mass[3]:.5f}  Z4 {mass[4]:.5f}; outside Z4 {1 - mass[4]:.5f}")
    early = w0 + wa
    print(f"  w0 + wa: mean {wstats(early, w)[0]:.3f}; 99.9th weighted percentile {wquant(early, w, 0.999):.3f}; max over samples {early.max():.3f}")
    # crossing (C8)
    zx = crossing_z(w0, wa)
    good = np.isfinite(zx)
    pneg = w[wa < 0].sum() / w.sum()
    lo, med, hi = (wquant(zx[good], w[good], q) for q in (0.16, 0.5, 0.84))
    print(f"  P(wa < 0) = {pneg:.5f}; crossing in the past for {w[good].sum() / w.sum():.4f} of the mass; z_x median {med:.3f}, 68% {lo:.3f} .. {hi:.3f}")
    summary[key] = dict(w0=m0, w0_sd=s0, wa=ma, wa_sd=sa, a_p=pv["a_p"], z_p=pv["z_p"], wp=pv["wp_mean"], wp_sd=pv["wp_sd"], corr=rho,
                        d=d, axis_deg=ang_axis, zones=mass, early_p999=wquant(early, w, 0.999), p_wa_neg=pneg, zx=(lo, med, hi),
                        map=(float(w0[i]), float(wa[i]), float(wp[i])))

print("\n=== summary and tolerances ===")
ref = {"PP": ("desi.table_w0wa_PP", "desi.pivot_PP", "desi.wcdm_PP"), "U3": ("desi.table_w0wa_U3", "desi.pivot_U3", "desi.wcdm_U3"),
       "DY5": ("desi.table_w0wa_DY5", "desi.pivot_DY5", "desi.wcdm_DY5"), "CMB": ("desi.table_w0wa_CMB", "desi.pivot_CMB", "desi.wcdm_CMB")}
for sn, (kt, kp, kw) in ref.items():
    k = f"{rel}_w0wa_{sn}"
    if k not in summary: continue
    S = summary[k]
    if rel == "dr1":
        t = fields(V[kt]["value"]); p = fnum(V[kp]["value"])
        (w0r, w0s), (war, was) = sym(t["w0"]), sym(t["wa"])
        g = [abs(S["w0"] - w0r) <= 0.005, abs(S["wa"] - war) <= 0.02, abs(S["w0_sd"] / w0s - 1) <= 0.1, abs(S["wa_sd"] / was - 1) <= 0.1,
             abs(S["z_p"] - p[0]) <= 0.02, abs(S["wp"] - p[1]) <= 0.005, abs(S["wp_sd"] / p[2] - 1) <= 0.1]
        print(f"{sn}: gate C0.2 w0 {S['w0']:.4f} vs {w0r} ({ok(g[0])}), wa {S['wa']:.3f} vs {war} ({ok(g[1])}), sd w0 {S['w0_sd']:.4f} vs {w0s:.3f} ({ok(g[2])}), sd wa {S['wa_sd']:.3f} vs {was:.3f} ({ok(g[3])}),"
              f" z_p {S['z_p']:.3f} vs {p[0]} ({ok(g[4])}), w_p {S['wp']:.4f} vs {p[1]} ({ok(g[5])}), sd w_p {S['wp_sd']:.4f} vs {p[2]} ({ok(g[6])}) -> gate {ok(all(g))}")
        kwk = f"{rel}_w_{sn}"
        if kwk in summary:
            wr = fnum(V[kw]["value"]); ws = wr[1] if len(wr) == 2 else (wr[1] + abs(wr[2])) / 2
            print(f"     wCDM chain w {summary[kwk]['w_mean']:.4f} +- {summary[kwk]['w_sd']:.4f} vs DESI {wr[0]} +- {ws:.3f}: mean {ok(abs(summary[kwk]['w_mean'] - wr[0]) <= 0.005)}, sd {ok(abs(summary[kwk]['w_sd'] / ws - 1) <= 0.1)}")
    if sn == "CMB": continue
    kwk = f"{rel}_w_{sn}"
    ratio = S["wp_sd"] / summary[kwk]["w_sd"] if kwk in summary else float("nan")
    print(f"  {sn}: C3 z_p {S['z_p']:.3f} | C4 |E[w_p]+1| = {abs(S['wp'] + 1):.4f} (<=0.025 {ok(abs(S['wp'] + 1) <= 0.025)}, <=0.06 {ok(abs(S['wp'] + 1) <= 0.06)}) | C5 |corr| {abs(S['corr']):.1e} ({ok(abs(S['corr']) < 0.01)}), sd(w_p)/sd(w_const chain) {ratio:.3f} ({ok(0.8 <= ratio <= 1.25)})"
          f" | C6 |d| = {abs(S['d']):.2f} ({'<=1' if abs(S['d']) <= 1 else '(1,2]' if abs(S['d']) <= 2 else '>2'}) | C7 outside Z4 {1 - S['zones'][4]:.4f} ({ok(1 - S['zones'][4] < 0.05)}), w0+wa p99.9 {S['early_p999']:.3f} (< -0.3 {ok(S['early_p999'] < -0.3)})"
          f" | C8 P(wa<0) {S['p_wa_neg']:.4f} ({ok(S['p_wa_neg'] >= 0.95)}), z_p in z_x 68% [{S['zx'][0]:.3f}, {S['zx'][2]:.3f}]: {ok(S['zx'][0] <= S['z_p'] <= S['zx'][2])}")
json.dump(summary, open(os.path.join(OUT, f"chain_summary_{rel}" + (f"_burn{round(100 * float(os.environ['BURN_ALL'])):03d}" if os.environ.get("BURN_ALL") else "") + ".json"), "w"), indent=1, default=float)

if rel == "dr2":
    # Addendum A1 item 5: DESI's own margestats in each DR2 run directory, a check of our loading (reported only)
    import re
    print("\n=== DR2: our means and sds vs DESI's .margestats (68% column) ===")
    for sn in ("PP", "U3", "DY5", "CMB"):
        k = f"{rel}_w0wa_{sn}"
        if k not in summary: continue
        ms = os.path.join(CHAINS, RUNS[k]["root"] + ".margestats")
        rows = {l.split()[0]: l[40:80].strip() for l in open(ms).read().splitlines()[1:] if l.strip()}
        for par, mine, sd in (("w", summary[k]["w0"], summary[k]["w0_sd"]), ("wa", summary[k]["wa"], summary[k]["wa_sd"])):
            v = fnum(rows[par].replace("\\pm", " ").replace("^{", " ").replace("}_{", " ").replace("}", " "))
            m, e = sym(v)
            print(f"  {sn} {par}: ours {mine:.4f} +- {sd:.4f} | DESI margestats {rows[par]} -> mean diff {mine - m:+.4f}, sd ratio {sd / e:.3f}")
