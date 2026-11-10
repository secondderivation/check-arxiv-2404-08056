"""C0.1 chain identity: SHA-256 of every local chain file against datasets/registry.json, and the
likelihoods, sampled dark-energy parameters and their priors read from each run's .updated.yaml.

Run: CHAINS=<dir> python3 -I checks/c01_identity.py <dr1|dr2> > checks/out/c01_identity_<rel>.out
"""
import hashlib, json, os, sys
import yaml
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RUNS, CHAINS

rel = sys.argv[1]
REG = json.load(open(os.path.join(os.path.dirname(__file__), "..", "..", "..", "datasets", "registry.json")))
entry = next(e for e in REG["entries"] if e["id"] == f"desi-{rel}-cosmo-chains")
bad = 0
for f in entry["distribution"]["files"]:
    p = os.path.join(CHAINS, f"desi-{rel}", f["path"])
    h = hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.exists(p) else "missing"
    if h != f["sha256"]:
        bad += 1
        print("MISMATCH", f["path"], h)
print(f"files pinned: {len(entry['distribution']['files'])}, mismatches: {bad} -> {'PASS' if bad == 0 else 'FAIL'}")

for key, r in RUNS.items():
    if not key.startswith(rel + "_"): continue
    y = yaml.safe_load(open(os.path.join(CHAINS, r["root"] + ".updated.yaml")))
    like = sorted(k.split(".")[-1] for k in (y.get("likelihood") or {}))
    post = y.get("post") or {}
    added = sorted(k.split(".")[-1] for k in ((post.get("add") or {}).get("likelihood") or {}))
    removed = sorted(k.split(".")[-1] for k in ((post.get("remove") or {}).get("likelihood") or {}))
    params = y.get("params", {})
    de = {}
    for n in ("w", "w0", "wa"):
        p = params.get(n)
        if isinstance(p, dict) and "prior" in p:
            de[n] = (p["prior"].get("min"), p["prior"].get("max"))
        elif p is not None:
            de[n] = p if not isinstance(p, dict) else p.get("value", "derived")
    de_model = ((y.get("theory") or {}).get("camb") or {}).get("extra_args", {}).get("dark_energy_model")
    wedge = [k for k in (y.get("prior") or {})]
    print(f"\n{key}\n  likelihoods: {', '.join(like)}" + (f"\n  post add: {', '.join(added)}" if added else "") + (f"\n  post remove: {', '.join(removed)}" if removed else ""))
    print(f"  dark energy params: {de}; CAMB dark_energy_model: {de_model}; extra prior blocks: {wedge or 'none'}")
