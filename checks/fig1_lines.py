"""C7 (figure part): positions of the lines the authors drew on Fig. 1 (src-v3/DESIpivotplot.pdf).

The figure is DESI VI's Fig. 6 (right) as an embedded raster (820 x 732 px) with the authors' lines
as vector strokes on top. The axes are calibrated from the raster's own tick marks (found below as
dark pixel columns/rows just outside the frame), then each vector line is converted to (w0, wa)
and compared with its stated locus at three heights (PROTOCOL C7: within 0.02 in w0).

Run: python3 -I checks/fig1_lines.py <src-v3/DESIpivotplot.pdf> > checks/out/fig1_lines.out
"""
import sys
import numpy as np
import pdfplumber

pdf = pdfplumber.open(sys.argv[1]); page = pdf.pages[0]
im = page.images[0]; st = im["stream"]
W, H = st.attrs["Width"], st.attrs["Height"]
a = np.frombuffer(st.get_data(), np.uint8)[: W * H * 3].reshape(H, W, 3).astype(int)
dark = a.sum(2) < 250
# frame: rows/columns that are dark over most of the plot
rows = [i for i, r in enumerate(dark.sum(1)) if r > 300]
cols = [i for i, c in enumerate(dark[60:580].sum(0)) if c > 400]
top, bottom, left, right = rows[0], rows[-1], cols[0], cols[-1]
xt = [i for i, v in enumerate(dark[bottom + 1: bottom + 9].sum(0)) if v >= 5 and i > left]
yt = [i for i, v in enumerate(dark[:, left - 8: left].sum(1)) if v >= 5]
print(f"raster {W}x{H}; frame rows {top}-{bottom}, cols {left}-{right}")
print("x tick columns (right of the frame's left edge):", xt, " y tick rows:", yt)
xv = [-1.0, -0.8, -0.6][: len(xt)]  # tick labels read from the figure: -1.0, -0.8, -0.6, -0.4
yv = [0.0, -1.0, -2.0, -3.0][: len(yt)]
fx = np.polyfit(xt[: len(xv)], xv, 1); fy = np.polyfit(yt[: len(yv)], yv, 1)
print(f"w0 per px {fx[0]:.6f} (residual {np.abs(np.polyval(fx, xt[:len(xv)]) - xv).max():.1e}); wa per px {fy[0]:.6f} (residual {np.abs(np.polyval(fy, yt[:len(yv)]) - yv).max():.1e})")
sx = (im["x1"] - im["x0"]) / W; sy = (im["bottom"] - im["top"]) / H

def to_data(x, y):  # page coordinates (top origin) -> (w0, wa)
    px, py = (x - im["x0"]) / sx, (y - im["top"]) / sy
    return np.polyval(fx, px), np.polyval(fy, py)

name = {(0.9321255, 0.1348295, 0.04655506): "red (today, w0 = -1)", (1.0, 0.2527924, 1.0): "magenta (early, w0 + wa = -1)",
        (0.0, 0.4624593, 0.7291366): "blue (w = -1 at z = 0.26)", (1.0, 0.577706, 0.0): "orange (w = -1 at z = 0.33)"}
locus = {"red": lambda wa: -1.0, "magenta": lambda wa: -1 - wa,
         "blue": lambda wa: -1 - wa * (1 - 1 / 1.26), "orange": lambda wa: -1 - wa * (1 - 1 / 1.33)}
worst = 0
for ln in page.lines:
    c = tuple(ln["stroking_color"]); lab = name.get(c)
    if not lab: continue
    (x0, y0), (x1, y1) = ln["pts"][0], ln["pts"][-1]
    p0, p1 = np.array(to_data(x0, y0)), np.array(to_data(x1, y1))
    key = lab.split()[0]
    devs = []
    for t in (0.1, 0.5, 0.9):
        w0, wa = p0 + t * (p1 - p0)
        devs.append(w0 - locus[key](wa))
    worst = max(worst, max(abs(d) for d in devs))
    print(f"{lab:<32} from ({p0[0]:+.3f}, {p0[1]:+.3f}) to ({p1[0]:+.3f}, {p1[1]:+.3f}); w0 minus locus at 10/50/90%: " + ", ".join(f"{d:+.4f}" for d in devs))
print(f"largest |deviation| in w0: {worst:.4f}  (tolerance 0.02): {'PASS' if worst <= 0.02 else 'FAIL'}")
