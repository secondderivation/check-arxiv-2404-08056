"""C1, C2, C3 (quoted part) and C10: the paper's quoted DESI numbers against DESI VI, and the
sigma conversions (G5, G10). DESI VI v1/v2 strings are checked when their TeX is given.

Run: python3 -I checks/quoted.py [<DESI VI v1 main tex> [<v2 main tex>]] > checks/out/quoted.out
"""
import json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import sigma_from_dchi2, half_up

V = json.load(open(os.path.join(os.path.dirname(__file__), "..", "paper_values.json")))
P, D = V["paper"], V["desi_vi_v3"]
ok = lambda b: "PASS" if b else "FAIL"

def nums(s):
    return [x for x in re.findall(r"[-+]?\d+(?:\.\d+)?", s)]

def fields(s):
    """'Om 0.3085 +- 0.0068; w0 -0.827 +- 0.063' -> {'Om': [...], 'w0': [...]} (label = first token)."""
    out = {}
    for part in s.split(";"):
        lab, rest = part.strip().split(" ", 1)
        out[lab] = nums(rest)
    return out

print("== C1 significances ==")
print("paper l.55:", P["intro.sigmas"]["value"], "| l.128:", P["sec2.sigmas_again"]["value"], "| DESI VI v3 l.907:", D["desi.sigmas"]["value"])
print("match:", ok(P["intro.sigmas"]["value"] == P["sec2.sigmas_again"]["value"] == D["desi.sigmas"]["value"]))
dchi = [float(x) for x in nums(D["desi.dchi2"]["value"])]
conv = [sigma_from_dchi2(-d, 2) for d in dchi]
printed = nums(D["desi.sigmas"]["value"])
for d, c, p in zip(dchi, conv, printed):
    print(f"  Delta chi2_MAP {d:+.1f}, k=2 -> {c:.4f} sigma -> half-up {half_up(c, 1)} vs printed {p}: {ok(half_up(c, 1) == p)}")
    print(f"     (k=1 would give {sigma_from_dchi2(-d, 1):.3f}; sqrt(-dchi2) = {(-d) ** 0.5:.3f})")
print("DESI VI describes them as Delta chi2_MAP between the w0waCDM MAP and the MAP at (w0, wa) = (-1, 0): a likelihood (posterior-maximum) ratio;",
      "the paper's 'based on a likelihood-ratio comparison' (l.128):", ok("likelihood" in P["sec2.likelihood_ratio"]["value"]))
print("DESI VI's stated MAP scatter on the significances:", D["desi.sigma_unc"]["value"], "sigma")

print("\n== C2 priors (text part; chain part in c0_controls.py) ==")
for k1, k2 in (("sec2.prior_w0", "desi.prior_w0"), ("sec2.prior_wa", "desi.prior_wa"), ("sec2.prior_sum", "desi.prior_sum")):
    a, b = P[k1]["value"], D[k2]["value"]
    print(f"  paper {a!r:>14}  DESI {b!r:>14}: {ok(a.replace(' ', '') in b.replace(' ', '') or a == b)}")

print("\n== C3 pivot redshifts (quoted) ==")
desi_zp = [nums(D[k]["value"])[0] for k in ("desi.pivot_PP", "desi.pivot_U3", "desi.pivot_DY5")]
print("DESI VI z_p (PP, U3, DY5):", desi_zp, "| paper: 0.26 (PP, DY5), 0.33 (U3):", ok(desi_zp == ["0.26", "0.33", "0.26"]))
for z, apr in (("0.26", "0.79"), ("0.33", "0.75")):
    a = 1 / (1 + float(z))
    print(f"  a_p = 1/(1+{z}) = {a:.5f} -> half-up {half_up(a, 2)} vs printed {apr}: {ok(half_up(a, 2) == apr)}")
print("DESI+CMB (no SN) pivot, not discussed by the paper:", D["desi.pivot_CMB"]["value"])

print("\n== C4 (quoted part): Fig. 2 parameters vs DESI VI ==")
lines = [tuple(float(v) for v in nums(t)) for t in P["fig2.lines"]["value"].split(";")]
eq = {k: (nums(D[f"desi.eq_w0_{k}"]["value"])[0], nums(D[f"desi.eq_wa_{k}"]["value"])[0]) for k in ("PP", "U3", "DY5")}  # DESI VI text equations l.880-895
tab = {k: fields(D[f"desi.table_w0wa_{k}"]["value"]) for k in ("PP", "U3", "DY5")}
for (k, (e0, ea)), (l0, la) in zip(eq.items(), lines):
    t0, ta = tab[k]["w0"][0], tab[k]["wa"][0]
    print(f"  {k}: Fig. 2 ({l0}, {la}) | DESI eq. ({e0}, {ea}) | DESI Table 3 ({t0}, {ta}): eq {ok((l0, la) == (float(e0), float(ea)))}, table {ok((l0, la) == (float(t0), float(ta)))}")
print(f"  DESI VI labels them: {D['desi.means_label']['value']!r} (l.877); the paper's Fig. 2 caption and l.100 call them 'best-fit'.")
print("  DESI's MAP for DESI+CMB+DESY5 (Appendix A l.1296):", D["desi.map_DY5"]["value"])
blobs = [tuple(float(v) for v in nums(t)) for t in P["fig2.blobs"]["value"].split(";")]
for (k, (l0, la)), (ab, wb), dk in zip(zip(eq, lines), blobs, ("desi.pivot_PP", "desi.pivot_U3", "desi.pivot_DY5")):
    wline = l0 + la * (1 - ab)
    wexact = l0 + la * (1 - 1 / (1 + float(nums(D[dk]["value"])[0])))
    print(f"  {k}: blob w = {wb}; DESI w_p = {nums(D[dk]['value'])[1]}; line at a={ab}: {wline:.4f}; line at exact 1/(1+z_p): {wexact:.4f};"
          f" |blob - line(a)| = {abs(wb - wline):.4f} <= 0.005: {ok(abs(wb - wline) <= 0.005)}")
print("  |w_p + 1| for the blobs:", [f"{abs(wb + 1):.3f}" for _, wb in blobs], "(tolerance 0.025 for '+-0.02', 0.06 for 'a few hundredths')")
print("  in units of DESI's sigma(w_p):", [f"{abs(float(nums(D[k]['value'])[1]) + 1) / float(nums(D[k]['value'])[2]):.2f}" for k in ("desi.pivot_PP", "desi.pivot_U3", "desi.pivot_DY5")])

print("\n== C5 (quoted part): sigma(w_p) vs constant-w sigma, DESI VI printed ==")
pairs = (("PP", "desi.pivot_PP", "desi.wcdm_PP"), ("U3", "desi.pivot_U3", "desi.wcdm_U3"), ("DY5", "desi.pivot_DY5", "desi.wcdm_DY5"), ("CMB", "desi.pivot_CMB", "desi.wcdm_CMB"))
for k, kp, kw in pairs:
    sp_ = float(nums(D[kp]["value"])[2]); wn = nums(D[kw]["value"])
    sw = (float(wn[1]) + abs(float(wn[2]))) / 2 if len(wn) == 3 else float(wn[1])
    print(f"  {k}: sigma(w_p) {sp_} / sigma(w) {sw:.4f} = {sp_ / sw:.3f} in [0.8, 1.25]: {ok(0.8 <= sp_ / sw <= 1.25)};"
          f"  w_p {nums(D[kp]['value'])[1]} vs constant w {wn[0]} (Gaussian identity: equal means; difference {float(nums(D[kp]['value'])[1]) - float(wn[0]):+.3f})")

print("\n== C10 ==")
print("(i) constant w consistent with LCDM in all cases (|w+1| < 2 sigma, sigma on the side towards -1):")
for k in ("desi.wcdm_DESI", "desi.wcdm_DESI_BBN_theta", "desi.wcdm_CMB", "desi.wcdm_PP", "desi.wcdm_U3", "desi.wcdm_DY5"):
    n = [float(x) for x in nums(D[k]["value"])]
    m = n[0]; s = (n[1] if m < -1 else abs(n[2])) if len(n) == 3 else n[1]
    print(f"  {D[k]['meaning']:<30} w = {D[k]['value']:<22} |w+1|/sigma = {abs(m + 1) / s:.2f}: {ok(abs(m + 1) / s < 2)}")
print("(ii) SDSS swap, DESI VI Appendix A (l.1284):", D["desi.sdss_sigmas"]["value"], "vs", D["desi.sigmas"]["value"])
red = [float(a) - float(b) for a, b in zip(nums(D["desi.sigmas"]["value"]), nums(D["desi.sdss_sigmas"]["value"]))]
print("   reductions:", [f"{r:.1f}" for r in red], "each in [0.25, 0.75]:", ok(all(0.25 <= r <= 0.75 for r in red)))
print("(iii) DESI VI on z=0.51 (l.576):", D["desi.lrg_offset"]["printed"][:160], "...")
print("(iv) curvature: DESI VI Table 3 kw0wa rows, read in c10_curvature below")
flat = {k: (fields(D[f"desi.table_w0wa_{k}"]["value"])["w0"][0], fields(D[f"desi.table_w0wa_{k}"]["value"])["wa"][0]) for k in ("PP", "U3", "DY5")}
for k in ("PP", "U3", "DY5"):
    f = fields(D[f"desi.kw0wa_{k}"]["value"].replace("10^3 Ok", "Ok"))
    m, s_ = float(f["Ok"][0]), float(f["Ok"][1])
    sa = (float(f["wa"][1]) + abs(float(f["wa"][2]))) / 2
    print(f"  {k}: 10^3 Omega_k = {m} +- {s_}: |Ok|/sigma = {abs(m) / s_:.2f} < 0.5: {ok(abs(m) / s_ < 0.5)};  w0 {f['w0'][0]} vs flat {flat[k][0]} (shift {abs(float(f['w0'][0]) - float(flat[k][0])) / float(f['w0'][1]):.2f} sigma), wa {f['wa'][0]} vs flat {flat[k][1]} (shift {abs(float(f['wa'][0]) - float(flat[k][1])) / sa:.2f} sigma)")

# DESI VI earlier versions: same strings present?
for path in sys.argv[1:]:
    txt = open(path, encoding="utf8", errors="replace").read()
    print(f"\n== DESI VI version check: {os.path.basename(os.path.dirname(path))} ==")
    for k, v in D.items():
        print(f"  {k:<28} {'present' if v['printed'] in txt else 'ABSENT'}")
