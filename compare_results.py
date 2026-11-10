"""Release replay checks for case 2404.08056, distinct from the frozen scientific claim criteria (see VALIDATION.md)."""
import re
from decimal import Decimal
from pathlib import Path

NUM = re.compile(r'(?<![\w.])[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?')


def text_compare(old, new, atol=1e-9, rtol=1e-9):
    if NUM.sub('#', old) != NUM.sub('#', new):
        return ['Output structure or verdict changed']
    # Compare the printed decimal values directly: binary floats can put an
    # exact boundary (e.g. -0.781 versus -0.780) just outside the stated tolerance.
    a, b = [Decimal(x) for x in NUM.findall(old)], [Decimal(x) for x in NUM.findall(new)]
    absolute, relative = Decimal(str(atol)), Decimal(str(rtol))
    if len(a) != len(b) or any(abs(x - y) > max(absolute, relative * max(abs(x), abs(y)))
                               for x, y in zip(a, b)):
        return ['Numeric replay tolerance exceeded']
    return []


# Correlations that are zero by construction print as round-off (e.g. 6.1e-16) and differ in sign and size across
# platforms; any |value| below 1e-12 is the same result.
TINY = re.compile(r'(?<![\w.])[-+]?\d+(?:\.\d+)?[eE]-(?:1[2-9]|[2-9]\d)\b')


def _zero_tiny(s):
    return TINY.sub('0', s)


def chain_compare(old, new):
    """Statistics of DESI's released chains (weights, getdist densities): printed values within 1e-3 absolute or
    1e-4 relative, including equality at either tolerance boundary."""
    return text_compare(_zero_tiny(old), _zero_tiny(new), atol=1e-3, rtol=1e-4)


REF = Path(__file__).resolve().parent  # the published package, never the work copy the checks run in


def compare_run(run, work, stdout):
    """Compares against the published output in the package root; `work` is unused, so a check that rewrote its
    work copy cannot turn the comparison into a self-comparison."""
    old = (REF / 'checks' / run['out']).read_text()
    if run['out'].endswith(('algebra.out', 'prior_volume.out')):
        return text_compare(old, stdout)  # symbolic algebra and seeded quadrature: deterministic
    return chain_compare(old, stdout)
