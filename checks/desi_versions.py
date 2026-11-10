"""Which DESI 2024 VI version printed each number the paper quotes: searches every desi_vi_v3 'printed'
string of paper_values.json in the v1, v2 and v3 main TeX files (exact substring). For the literature
and version notes (protocol section 5); never changes an outcome.

Run: python3 -I checks/desi_versions.py <v1 FinalReview_adds_kp7_paper1.tex> <v2 FinalReview_adds_kp7_paper1.tex> <v3 JournalReview_adds_v2_kp7_paper1.tex> > checks/out/desi_versions.out
"""
import hashlib, json, os, sys
V = json.load(open(os.path.join(os.path.dirname(__file__), "..", "paper_values.json")))["desi_vi_v3"]
texts = []
for p in sys.argv[1:4]:
    raw = open(p, "rb").read()
    print(os.path.basename(os.path.dirname(p)), hashlib.sha256(raw).hexdigest())
    texts.append(raw.decode("utf-8"))
miss = 0
for k, v in V.items():
    pr = v.get("printed")
    if not pr: continue
    found = ["yes" if pr in t else "no" for t in texts]
    if found != ["yes"] * 3: miss += 1
    print(f"{k:<28} v1 {found[0]:<3} v2 {found[1]:<3} v3 {found[2]:<3}  {pr[:70]}")
print(f"\n{miss} of {sum(1 for v in V.values() if v.get('printed'))} printed strings are not in all three versions")
