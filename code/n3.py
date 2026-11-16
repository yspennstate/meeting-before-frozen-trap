"""Three states (Section 3 of the paper).

With w_z the conductance of the edge opposite z and unnormalised weights p_1, p_2, p_3 (S = p_1 + p_2 + p_3),
    h_z = a_z / (a_z + b_z),   a_z = w_z (S - p_z),   b_z = sum_{x != z} w_x p_x,
and D = F - 1/3 - (2/3) sum pi^3 = 2 e_3 (sum_z h_z - 1) / S^3.  This script checks, in exact arithmetic:
  1. the closed form of h_z and of D against the general solver meeting.py at random integer points;
  2. the identity  sum_z a_z d_x d_y - d_1 d_2 d_3 = S sum_z (prod_{x != z} p_x) w_z (sum_{x != z} w_x^2)
                                                    + 2 (S e_2 - 2 e_3) w_1 w_2 w_3,   d = a + b;
  3. that every pair of conductances carries a monomial of the numerator in those two alone (strictness on the path);
  4. that D -> 0 along w = (eps, eps, 1) for fixed p (sharpness).
Exits with status 1 on any failure.

    python n3.py
"""
import random
import sys
from fractions import Fraction as Fr

import sympy as sp

from meeting import bound, probability

fails = 0


def check(label, ok):
    global fails
    fails += not ok
    print("%-80s %s" % (label, "ok" if ok else "FAILED"))


p1, p2, p3, w1, w2, w3 = sp.symbols("p1 p2 p3 w1 w2 w3", positive=True)
p, w = [p1, p2, p3], [w1, w2, w3]
S = p1 + p2 + p3
e2 = p1 * p2 + p1 * p3 + p2 * p3
e3 = p1 * p2 * p3
a = [w[z] * (S - p[z]) for z in range(3)]
b = [sum(w[x] * p[x] for x in range(3) if x != z) for z in range(3)]
d = [a[z] + b[z] for z in range(3)]
num = sp.expand(sum(a[z] * d[(z + 1) % 3] * d[(z + 2) % 3] for z in range(3)) - d[0] * d[1] * d[2])
closed = sp.expand(S * sum(sp.prod([p[x] for x in range(3) if x != z]) * w[z] * sum(w[x] ** 2 for x in range(3) if x != z)
                           for z in range(3)) + 2 * (S * e2 - 2 * e3) * w1 * w2 * w3)
check("numerator equals the closed form", sp.expand(num - closed) == 0)
P = sp.Poly(num, *p, *w)
coeffs = [c for _, c in P.terms()]
check("numerator: %d terms, all coefficients positive (min %s)" % (len(coeffs), min(coeffs)), min(coeffs) > 0)
check("S e_2 - 2 e_3 has positive coefficients", min(sp.Poly(sp.expand(S * e2 - 2 * e3), *p).coeffs()) > 0)
for i, j in ((0, 1), (0, 2), (1, 2)):
    k = 3 - i - j
    alone = [m for m, _ in P.terms() if m[3 + k] == 0]
    check("a monomial in w_%d, w_%d alone (absent edge w_%d)" % (i + 1, j + 1, k + 1), len(alone) > 0)

rnd = random.Random(3)
worst = 0
for t in range(30):
    pv = [rnd.randint(1, 40) for _ in range(3)]
    wv = [rnd.randint(0 if t % 3 == 0 else 1, 40) for _ in range(3)]
    if sum(1 for x in wv if x > 0) < 2:
        wv[0] = wv[1] = 1
    c = [[0, wv[2], wv[1]], [wv[2], 0, wv[0]], [wv[1], wv[0], 0]]
    exact = probability(pv, c) - bound(pv)
    subs = dict(zip(p + w, pv + wv))
    hz = [Fr(int(a[z].subs(subs)), int(d[z].subs(subs))) for z in range(3)]
    sv = sum(pv)
    formula = 2 * Fr(pv[0] * pv[1] * pv[2]) * (sum(hz) - 1) / Fr(sv) ** 3
    worst += exact != formula
    worst += not exact > 0
check("closed form of D equals the general solver, and D > 0, at 30 integer points", worst == 0)

pv = [2, 3, 7]
prev = None
mono = True
for k in range(1, 7):
    eps = Fr(1, 10 ** k)
    c = [[0, 1, eps], [1, 0, eps], [eps, eps, 0]]
    dval = probability(pv, c) - bound(pv)
    print("   eps = 1e-%d   D = %.3e" % (k, float(dval)))
    mono &= dval > 0 and (prev is None or dval < prev)
    prev = dval
check("D decreases to 0 along w = (eps, eps, 1), p = (2, 3, 7)", mono and prev < Fr(1, 10 ** 5))
sys.exit(1 if fails else 0)
