"""M1 reference values: numbers printed in DESI DR2 Results II (arXiv:2503.14738v3 main.tex), transcribed
with line and SHA-256 asserts, compared with our DR2 chain results (checks/out/chain_summary_dr2.json).
A data trace only (G17).

Run: python3 -I checks/m1_dr2_ref.py <path to 2503.14738v3 main.tex> > checks/out/m1_dr2_ref.out
"""
import hashlib, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import sigma_from_dchi2, half_up

tex = sys.argv[1]
raw = open(tex, "rb").read()
print("main.tex sha256", hashlib.sha256(raw).hexdigest())
L = raw.decode("utf-8").splitlines()
def line(n, must):
    assert must in L[n - 1], (n, must)
    return L[n - 1]
l926 = line(926, r"For DESI+CMB, $\zpiv=0.53$ and $\wpiv=w(\zpiv)=-1.024\pm0.043$")
assert r"we find $\zpiv=0.31$ and $\wpiv=-0.954\pm0.024$" in l926
l894 = line(894, r"The $\dchisq$ values are $-10.7$, $-17.4$, and $-21.0$")
assert r"$2.8\sigma$, $3.8\sigma$, and $4.2\sigma$" in l894
l907 = line(907, r"DESI+CMB & $-12.5$ & $3.1\sigma$")
l1119 = line(1119, "This constitutes a weak coincidence")
l857 = line(857, "the values of Bayesian quantities such as the model evidence will always depend on the particular choice used")

S = json.load(open(os.path.join(os.path.dirname(__file__), "out", "chain_summary_dr2.json")))
print("\n== pivots (DESI DR2 l.926 vs our DR2 chains) ==")
for sn, zp, wp, e in (("CMB", 0.53, -1.024, 0.043), ("DY5", 0.31, -0.954, 0.024)):
    s = S[f"dr2_w0wa_{sn}"]
    print(f"{sn}: DESI z_p {zp}, w_p {wp} +- {e} | ours z_p {s['z_p']:.4f} ({half_up(s['z_p'], 2)}), w_p {s['wp']:.4f} +- {s['wp_sd']:.4f}"
          f" -> z_p {'PASS' if half_up(s['z_p'], 2) == str(zp) else 'FAIL'}, w_p {'PASS' if abs(s['wp'] - wp) <= 0.005 else 'FAIL'}, sd {'PASS' if abs(s['wp_sd'] / e - 1) <= 0.1 else 'FAIL'}")
print("DESI DR2 says only that the PP and U3 pivots are 'slightly larger' with larger errors; ours:",
      ", ".join(f"{sn} z_p {S[f'dr2_w0wa_{sn}']['z_p']:.3f}, w_p {S[f'dr2_w0wa_{sn}']['wp']:.4f} +- {S[f'dr2_w0wa_{sn}']['wp_sd']:.4f}" for sn in ("PP", "U3")))
print("\n== significances (l.894, l.907): sigma from Delta chi2_MAP with k = 2 ==")
for lab, d, p in (("PP", -10.7, "2.8"), ("U3", -17.4, "3.8"), ("DY5", -21.0, "4.2"), ("CMB", -12.5, "3.1")):
    sg = sigma_from_dchi2(-d, 2)
    print(f"{lab}: Delta chi2 {d} -> {sg:.4f} sigma -> {half_up(sg, 1)} vs printed {p}: {'PASS' if half_up(sg, 1) == p else 'FAIL'}")
print("\n== DESI DR2's own words on the pivot reading (l.1119) and prior dependence (l.857), quoted for literature.md ==")
print("l.1119:", re.search(r"We note that the degeneracy direction.*?observation\.", l1119).group(0))
print("l.857: ...", re.search(r"While other justifiable choices.*?approach\.", l857).group(0))
