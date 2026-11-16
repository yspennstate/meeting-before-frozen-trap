"""Four states (Section 4 of the paper): the certificate built with SymPy polynomial arithmetic.

For a trap z, the pair equation for the unordered pair {x, y} (x, y != z) is multiplied by p_x p_y, giving a 3 x 3
system A_z h_z = b_z with polynomial entries in the weights p and the conductances c of the support graph.  With
Q = prod_z det A_z and N = 3 F S^3 - S^3 - 2 sum p^3 (so N = 3 D S^3),
    R = N Q = 3 sum_z p_z [ (sum p^2) Q + 2 sum_{x<y} p_x p_y (adj(A_z) b_z)_{xy} prod_{w != z} det A_w ]
              - (S^3 + 2 sum p^3) Q
is a polynomial.  The script checks that every det A_z and R has only positive coefficients, prints the counts quoted
in Table 1, and checks N = R / Q exactly against the general solver meeting.py at random integer points.
Exits with status 1 on any failure.

    python n4_cramer.py star path cycle paw diamond k4      (k4 takes several minutes and about 1 GB)
"""
import itertools
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

from meeting import probability

GRAPHS = {"star": [(0, 1), (0, 2), (0, 3)], "path": [(0, 1), (1, 2), (2, 3)],
          "cycle": [(0, 1), (1, 2), (2, 3), (3, 0)], "paw": [(0, 1), (1, 2), (2, 0), (2, 3)],
          "diamond": [(0, 1), (1, 2), (2, 0), (1, 3), (2, 3)],
          "k4": [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]}


def build(E, n=4):
    p = sp.symbols("p0:%d" % n, positive=True)
    c = {}
    for a, b in E:
        c[(a, b)] = c[(b, a)] = sp.Symbol("c%d%d" % (min(a, b), max(a, b)), positive=True)
    gens = (*p, *sorted(set(c.values()), key=str))
    P = lambda e: sp.Poly(e, *gens, domain="ZZ")
    dets, tops = [], []
    for z in range(n):
        others = [x for x in range(n) if x != z]
        pairs = [frozenset(q) for q in itertools.combinations(others, 2)]
        idx = {q: i for i, q in enumerate(pairs)}
        m = len(pairs)
        A = [[P(0) for _ in range(m)] for _ in range(m)]
        b = [P(0) for _ in range(m)]
        for q in pairs:
            i = idx[q]
            x, y = tuple(q)
            for mover, other in ((x, y), (y, x)):
                for v in range(n):
                    if (mover, v) in c:
                        r = P(c[(mover, v)] * p[other])  # rate c / p_mover, times p_x p_y
                        A[i][i] += r
                        if v == other:
                            b[i] += r
                        elif v != z:
                            A[i][idx[frozenset((v, other))]] -= r
        M = sp.Matrix(m, m, lambda i, j: A[i][j].as_expr())
        det = P(M.det(method="berkowitz"))
        adj = M.adjugate(method="berkowitz")
        sol = [sum((P(adj[i, j]) * b[j] for j in range(m)), P(0)) for i in range(m)]
        top = P(0)
        for q in pairs:
            x, y = tuple(q)
            top += P(2 * p[x] * p[y]) * sol[idx[q]]
        dets.append(det)
        tops.append(top)
    s2 = P(sum(x ** 2 for x in p))
    Q = P(1)
    for dz in dets:
        Q *= dz
    R = P(0)
    for z in range(n):
        rest = P(1)
        for k, dz in enumerate(dets):
            if k != z:
                rest *= dz
        R += P(3 * p[z]) * (s2 * Q + tops[z] * rest)
    S = sum(p)
    R -= P(S ** 3 + 2 * sum(x ** 3 for x in p)) * Q
    return p, gens, dets, R


def peval(poly, vals):
    tot = 0
    for mon, co in poly.terms():
        t = int(co)
        for v, e in zip(vals, mon):
            if e:
                t *= v ** e
        tot += t
    return tot


if __name__ == "__main__":
    fails = 0
    for name in sys.argv[1:] or list(GRAPHS):
        E = GRAPHS[name]
        t0 = time.time()
        p, gens, dets, R = build(E)
        dc = [d.coeffs() for d in dets]
        rc = R.coeffs()
        ok = all(min(x) > 0 for x in dc) and min(rc) > 0
        print("%-8s det A_z terms %-18s min det coeff %s | R terms %8d, degree %d, min coeff %s, negative %d | %.0f s"
              % (name, sorted(len(x) for x in dc), min(min(x) for x in dc), len(rc), R.total_degree(), min(rc),
                 sum(1 for x in rc if x < 0), time.time() - t0), flush=True)
        names = [str(g) for g in gens]
        rnd = random.Random(2026)
        agree = 0
        for _ in range(10):
            pv = [rnd.randint(1, 50) for _ in range(4)]
            cv = {tuple(sorted(e)): rnd.randint(1, 50) for e in E}
            vals = pv + [cv[tuple(int(ch) for ch in g[1:])] for g in names[4:]]
            cm = [[0] * 4 for _ in range(4)]
            for (i, j), v in cv.items():
                cm[i][j] = cm[j][i] = v
            S = sum(pv)
            N = 3 * probability(pv, cm) * S ** 3 - S ** 3 - 2 * sum(x ** 3 for x in pv)
            Q = 1
            for dz in dets:
                Q *= peval(dz, vals)
            agree += N == Fr(peval(R, vals), Q)
        print("%-8s exact identity N = R / Q at %d of 10 integer points" % (name, agree), flush=True)
        fails += (not ok) + (agree != 10)
    print("ALL PASS" if fails == 0 else "FAILURES: %d" % fails)
    sys.exit(1 if fails else 0)
