"""Shared code for the 2404.08056 checks: chain loading, pivot, zones, crossing, prior tapers.

Written from conventions.json; no author code exists (G13). Chain roots are listed in
checks/chains.json (key -> directory and file root under CHAINS: the environment variable, else
data/chains/ beside checks/ in the published package, else the review's local copy of the files
pinned in datasets.json).
"""
import json, os, sys
import numpy as np
from scipy.special import erfinv
from scipy.stats import chi2 as chi2dist

HERE = os.path.dirname(os.path.abspath(__file__))
CHAINS = os.environ.get("CHAINS") or (os.path.join(HERE, "..", "data", "chains") if os.path.isdir(os.path.join(HERE, "..", "data", "chains"))
                                     else "/home/claude/work/chains")
# Machine-readable outputs (JSON, figures): checks/out/, or with --replay (reproduce.py in the published package)
# a separate folder beside the work copy, so a replay never rewrites the published files it is compared against.
REPLAY = "--replay" in sys.argv
OUT = os.path.join(HERE, "..", "..", "replay-out") if REPLAY else os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
RUNS = json.load(open(os.path.join(HERE, "chains.json"))) if os.path.exists(os.path.join(HERE, "chains.json")) else {}
BURN = float(os.environ.get("BURN", "0.3"))  # raw Cobaya chains
BURN_ALL = os.environ.get("BURN_ALL")  # Addendum A1 variant: one burn-in fraction for every run
BOX = dict(w0=(-3.0, 1.0), wa=(-3.0, 2.0))  # plus w0 + wa < 0


def load(key, burn=None):
    """getdist MCSamples for a run key (conventions.json `chains`)."""
    from getdist import loadMCSamples
    root = os.path.join(CHAINS, RUNS[key]["root"])
    if burn is None:
        # Addendum A1: DESI's post-importance files already had 30 % removed (burn 0 there)
        burn = float(BURN_ALL) if BURN_ALL is not None else RUNS[key].get("burn", BURN)
    import contextlib, sys
    with contextlib.redirect_stdout(sys.stderr):  # getdist prints file names
        return loadMCSamples(root, settings={"ignore_rows": burn})


def cols(s, *names):
    """Columns by name, accepting the CAMB/Cobaya alias w for w0."""
    out = []
    p = s.getParams()
    for n in names:
        alias = {"w0": ["w0", "w"], "wa": ["wa"]}.get(n, [n])
        for al in alias:
            if hasattr(p, al):
                out.append(np.asarray(getattr(p, al)))
                break
        else:
            raise KeyError(n)
    return out


def wstats(x, wt):
    m = np.average(x, weights=wt)
    return m, np.sqrt(np.average((x - m) ** 2, weights=wt))


def wcov(x, y, wt):
    mx, my = np.average(x, weights=wt), np.average(y, weights=wt)
    return np.average((x - mx) * (y - my), weights=wt)


def wquant(x, wt, q):
    i = np.argsort(x)
    c = np.cumsum(wt[i]) / wt.sum()
    return np.interp(q, c, x[i])


def pivot(w0, wa, wt):
    """DETF decorrelation pivot: returns a_p, z_p, mean and sd of w_p (conventions.json `pivot`)."""
    c0a, caa, c00 = wcov(w0, wa, wt), wcov(wa, wa, wt), wcov(w0, w0, wt)
    one_minus_ap = -c0a / caa
    ap = 1 - one_minus_ap
    wp = w0 + wa * one_minus_ap
    m, s = wstats(wp, wt)
    return dict(a_p=ap, z_p=1 / ap - 1, wp_mean=m, wp_sd=s, wp_sd_formula=np.sqrt(c00 - c0a**2 / caa))


def zone(w0, wa):
    early = w0 + wa
    z = np.full(w0.shape, 4)
    z[(w0 < -1) & (early < -1)] = 1
    z[(w0 >= -1) & (early >= -1)] = 2
    z[(w0 < -1) & (early >= -1)] = 3
    return z


def crossing_z(w0, wa):
    with np.errstate(divide="ignore", invalid="ignore"):
        ax = 1 + (1 + w0) / wa
        zx = 1 / ax - 1
    zx[~((ax > 0) & (ax <= 1))] = np.nan
    return zx


def taper(name, w0, wa):
    """Prior variants T1-T5 (conventions.json `prior_variants`); returns t in [0, 1]."""
    w0, wa = np.asarray(w0, float), np.asarray(wa, float)
    wmin = np.minimum(w0, w0 + wa)
    if name in ("T1", "T2"):
        s = 0.5 if name == "T1" else 0.25
        return np.where(wmin >= -1, 1.0, np.exp(-((wmin + 1) ** 2) / (2 * s**2)))
    if name in ("T3", "T4"):
        s0, sa = (0.5, 1.0) if name == "T3" else (0.25, 0.5)
        return np.exp(-((w0 + 1) ** 2) / (2 * s0**2) - wa**2 / (2 * sa**2))
    if name == "T5":
        return ((w0 >= -1) & (w0 + wa >= -1)).astype(float)
    if name == "T0":
        return np.ones_like(w0)
    if name.startswith("CUT"):  # control C0.3: wa > value
        return (wa > float(name[3:])).astype(float)
    raise KeyError(name)


def in_box(w0, wa):
    return (w0 >= BOX["w0"][0]) & (w0 <= BOX["w0"][1]) & (wa >= BOX["wa"][0]) & (wa <= BOX["wa"][1]) & (w0 + wa < 0)


def prior_mean_quad(name, n=2001):
    """E_pi0[t] by midpoint quadrature on the DESI box."""
    g0 = np.linspace(*BOX["w0"], n + 1); g0 = (g0[1:] + g0[:-1]) / 2
    ga = np.linspace(*BOX["wa"], n + 1); ga = (ga[1:] + ga[:-1]) / 2
    W0, WA = np.meshgrid(g0, ga, indexing="ij")
    m = in_box(W0, WA)
    return taper(name, W0[m], WA[m]).mean(), m.mean() * 20.0


def prior_mean_mc(name, n=10_000_000, seed=20261018):
    rng = np.random.default_rng(seed)
    w0 = rng.uniform(*BOX["w0"], n); wa = rng.uniform(*BOX["wa"], n)
    m = in_box(w0, wa)
    return taper(name, w0[m], wa[m]).mean()


def sigma_from_dchi2(d, k=2):
    p = chi2dist.sf(d, k)
    return np.sqrt(2) * erfinv(1 - p)


def level_to_sigma(level):
    return np.sqrt(2) * erfinv(level)


def gauss_level(w0, wa, wt, point=(-1.0, 0.0)):
    mu = np.array([np.average(w0, weights=wt), np.average(wa, weights=wt)])
    C = np.array([[wcov(w0, w0, wt), wcov(w0, wa, wt)], [wcov(w0, wa, wt), wcov(wa, wa, wt)]])
    d = np.array(point) - mu
    d2 = d @ np.linalg.solve(C, d)
    lvl = 1 - np.exp(-d2 / 2)
    return d2, lvl, np.sqrt(2) * erfinv(min(lvl, 1 - 1e-16))


def half_up(x, nd):
    from decimal import Decimal, ROUND_HALF_UP
    q = Decimal(1).scaleb(-nd)
    return str(Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP))
