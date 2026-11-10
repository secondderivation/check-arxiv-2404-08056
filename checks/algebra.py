"""C5, C7, C8 algebra (G6, G9): re-derive the pivot, zone and maximum-density statements.

Run: python3 -I checks/algebra.py > checks/out/algebra.out
"""
import sympy as sp

w0, wa, a, ap, ax = sp.symbols("w0 w_a a a_p a_x", real=True)
S00, S0a, Saa = sp.symbols("S00 S0a Saa", real=True)
m0, ma = sp.symbols("mu0 mu_a", real=True)
ok = lambda b: "PASS" if b else "FAIL"

print("== C5 pivot transformation ==")
w = w0 + wa * (1 - a)
wp = w0 + wa * (1 - ap)
# (i) decorrelation: Cov(wp, wa) = S0a + (1-ap) Saa = 0
cov_p_a = S0a + (1 - ap) * Saa
sol = sp.solve(sp.Eq(cov_p_a, 0), ap)[0]
print("Cov(w_p, w_a) = ", cov_p_a, "-> a_p =", sol, " i.e. 1 - a_p = -S0a/Saa:", ok(sp.simplify(1 - sol + S0a / Saa) == 0))
var_p = sp.simplify((S00 + 2 * (1 - ap) * S0a + (1 - ap) ** 2 * Saa).subs(ap, sol))
print("Var(w_p) at the pivot =", var_p, ":", ok(sp.simplify(var_p - (S00 - S0a**2 / Saa)) == 0))
print("Var(w_p) <= Var(w0) (Saa > 0):", ok(sp.simplify(S00 - var_p - S0a**2 / Saa) == 0), "(difference S0a^2/Saa >= 0)")
# (ii) wa marginal unchanged: the map leaves wa as a coordinate
J = sp.Matrix([wp, wa]).jacobian([w0, wa])
print("Jacobian of (w0, wa) -> (w_p, wa):", J.tolist(), "det =", J.det(), ":", ok(J.det() == 1))
print("wa unchanged, w0 shifted by wa(1-a_p) at fixed wa: horizontal slide in the (w0, wa) plane:", ok(J[1, :].tolist() == [[0, 1]]))
# (iii) Gaussian: conditional mean and variance of w0 given wa = 0 (constant-w model)
cond_mean = m0 + S0a / Saa * (0 - ma)
cond_var = S00 - S0a**2 / Saa
mean_p = (m0 + ma * (1 - ap)).subs(ap, sol)
print("Gaussian: Var(w0 | wa=0) == Var(w_p):", ok(sp.simplify(cond_var - var_p) == 0))
print("Gaussian: E(w0 | wa=0) == E(w_p) (pivot mean equals the constant-w estimate):", ok(sp.simplify(cond_mean - mean_p) == 0))
# (iv) the line w_p = -1 in the (w0, wa) plane
line = sp.solve(sp.Eq(wp, -1), w0)[0]
print("w_p = -1 line: w0 =", line, "; passes through LCDM (-1, 0):", ok(line.subs(wa, 0) == -1),
      "; slope d wa/d w0 =", sp.simplify(1 / sp.diff(line, wa)))
for z in ("0.26", "0.33"):
    A = 1 / (1 + sp.Rational(z))
    print(f"  z_p={z}: a_p={sp.N(A, 6)}, slope d wa/d w0 = {sp.N(-1 / (1 - A), 6)}")

print("\n== C7 zones ==")
print("w(a) is linear in a, so on (0,1] its extremes are w(a->0) = w0 + wa and w(1) = w0:", w.subs(a, 0), ",", w.subs(a, 1))
print("Z1 always phantom: w0 < -1 and w0 + wa < -1   (lower left of both lines)")
print("Z2 never phantom: w0 >= -1 and w0 + wa >= -1  (upper right)")
print("Z3 phantom today only: w0 < -1, w0 + wa >= -1  => wa > -1 - w0 > 0 (top left)")
print("Z4 phantom early only: w0 >= -1, w0 + wa < -1 => wa < -1 - w0 <= 0")
print("Red line (today) w0 = -1 is vertical; magenta line (early) wa = -1 - w0 has slope -1 through (-1, 0)")

print("\n== C8 crossing and maximum of rho_de ==")
lnrho = -3 * (1 + w0 + wa) * sp.log(a) - 3 * wa * (1 - a)
dlnrho = sp.simplify(a * sp.diff(lnrho, a))
print("ln rho_de = ", lnrho, ";  d ln rho/d ln a =", sp.factor(dlnrho), ":", ok(sp.simplify(dlnrho + 3 * (1 + w)) == 0))
axs = sp.solve(sp.Eq(w, -1), a)[0]
print("crossing a_x =", sp.simplify(axs), ":", ok(sp.simplify(axs - (1 + (1 + w0) / wa)) == 0))
d2 = sp.simplify(sp.diff(dlnrho, a).subs(a, axs))
print("d/da (d ln rho / d ln a) at a_x =", d2, " -> negative for wa < 0: the crossing is a maximum of rho_de")
print("For wa < 0, w increases monotonically with a for all a > 0 (dw/da = -wa > 0), so w < -1 before a_x and > -1 after,",
      "i.e. rho_de rises before and falls after: a_x is the global maximum over a in (0, inf) of the CPL extrapolation.")
