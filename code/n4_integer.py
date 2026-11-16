"""Four states: a second, independent construction of the certificate of Section 4.

Integer polynomials are dictionaries from packed exponents to integer coefficients (base-32 Kronecker packing; every
monomial that occurs has total degree at most 27 < 32, so packed exponents never carry).  Determinants are expanded
over permutations and adjugates by signed minors; nothing is shared with n4_cramer.py.  In this form
    R = N Q = 6 [ sum_z p_z (v_z . adj(A_z) b_z) prod_{w != z} det A_w - e_3(p) Q ],   v_z(ij) = p_i p_j,
which equals the expression in n4_cramer.py after expanding S^3.  With --compare the two constructions are checked to
agree coefficient by coefficient.  Exits with status 1 on any failure.

    python n4_integer.py [--compare] C4 K4 P4 star paw diamond
"""
import argparse
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations, permutations

from meeting import probability

GRAPHS = {"C4": [(0, 1), (1, 2), (2, 3), (0, 3)], "K4": list(combinations(range(4), 2)),
          "P4": [(0, 1), (1, 2), (2, 3)], "star": [(0, 1), (0, 2), (0, 3)],
          "paw": [(0, 1), (0, 2), (1, 2), (2, 3)], "diamond": [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]}
TO_CRAMER = {"C4": "cycle", "K4": "k4", "P4": "path", "star": "star", "paw": "paw", "diamond": "diamond"}
BASE = 32
ONE = {0: 1}


def var(i):
    return {BASE ** i: 1}


def add(*items):
    out = {}
    for poly, scale in items:
        for e, c in poly.items():
            out[e] = out.get(e, 0) + scale * c
    return {e: c for e, c in out.items() if c}


def mul(a, b):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            out[ea + eb] = out.get(ea + eb, 0) + ca * cb
    return {e: c for e, c in out.items() if c}


def prod(items):
    out = ONE
    for x in items:
        out = mul(out, x)
    return out


def det(a):
    n = len(a)
    terms = []
    for perm in permutations(range(n)):
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        terms.append((prod([a[i][perm[i]] for i in range(n)]), sign))
    return add(*terms)


def evaluate(poly, values):
    tot = 0
    for code, coeff in poly.items():
        t = coeff
        for x in values:
            t *= x ** (code % BASE)
            code //= BASE
        tot += t
    return tot


def certificate(name):
    edges = GRAPHS[name]
    p = [var(i) for i in range(4)]
    c = [[{} for _ in range(4)] for _ in range(4)]
    for t, (i, j) in enumerate(edges):
        c[i][j] = c[j][i] = var(4 + t)
    D, H = [], []
    for z in range(4):
        pairs = list(combinations([i for i in range(4) if i != z], 2))
        index = {s: t for t, s in enumerate(pairs)}
        a = [[{} for _ in pairs] for _ in pairs]
        b = []
        for row, (i, j) in enumerate(pairs):
            a[row][row] = add((mul(p[j], add(*[(c[i][k], 1) for k in range(4) if k != i])), 1),
                              (mul(p[i], add(*[(c[j][k], 1) for k in range(4) if k != j])), 1))
            b.append(mul(add((p[i], 1), (p[j], 1)), c[i][j]))
            for k in range(4):
                if k not in (i, j, z):
                    a[row][index[tuple(sorted((k, j)))]] = add((mul(p[j], c[i][k]), -1))
                    a[row][index[tuple(sorted((i, k)))]] = add((mul(p[i], c[j][k]), -1))
        numer = []
        for row in range(3):
            pieces = []
            for col in range(3):
                minor = [[a[i][j] for j in range(3) if j != row] for i in range(3) if i != col]
                pieces.append((mul(det(minor), b[col]), (-1) ** (row + col)))
            numer.append(add(*pieces))
        D.append(det(a))
        H.append(add(*[(prod([p[i], p[j], num]), 1) for (i, j), num in zip(pairs, numer)]))
    Q = prod(D)
    e3 = add(*[(prod([p[i] for i in t]), 1) for t in combinations(range(4), 3)])
    R = add((mul(e3, Q), -6))
    for z in range(4):
        R = add((R, 1), (prod([p[z], H[z]] + [D[w] for w in range(4) if w != z]), 6))
    return edges, D, Q, R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("graphs", nargs="*")
    ap.add_argument("--compare", action="store_true")
    args = ap.parse_args()
    fails = 0
    for name in args.graphs or list(GRAPHS):
        t0 = time.time()
        edges, D, Q, R = certificate(name)
        ok = all(min(d.values()) > 0 for d in D) and min(R.values()) > 0
        agree = 0
        for weights, ev in (([1] * 4, [1] * len(edges)), ([1, 2, 3, 5], list(range(1, len(edges) + 1)))):
            cm = [[0] * 4 for _ in range(4)]
            for (i, j), v in zip(edges, ev):
                cm[i][j] = cm[j][i] = v
            S = sum(weights)
            N = 3 * probability(weights, cm) * S ** 3 - S ** 3 - 2 * sum(x ** 3 for x in weights)
            agree += N == Fr(evaluate(R, weights + ev), evaluate(Q, weights + ev))
        line = "%-8s det A_z terms %-18s R terms %8d, min coeff %d, negative %d, identity %d/2" % (
            name, sorted(len(d) for d in D), len(R), min(R.values()), sum(v < 0 for v in R.values()), agree)
        fails += (not ok) + (agree != 2)
        if args.compare:
            from n4_cramer import GRAPHS as CG, build
            ps, gens, sdets, sR = build(CG[TO_CRAMER[name]])
            pos = {str(x): i for i, x in enumerate(ps)}
            for k, (i, j) in enumerate(edges):
                pos["c%d%d" % (min(i, j), max(i, j))] = 4 + k
            enc = lambda poly: {sum(e * BASE ** pos[str(g)] for g, e in zip(gens, mon)): int(co) for mon, co in poly.terms()}
            # vertex labels of the two graph lists may differ by a relabelling; compare the invariant counts then,
            # and the full coefficient table when the labelled edge sets coincide
            same = {tuple(sorted(e)) for e in CG[TO_CRAMER[name]]} == {tuple(sorted(e)) for e in edges}
            if same:
                match = enc(sR) == R and sorted(map(len, map(enc, sdets))) == sorted(map(len, D))
                line += " | coefficient-by-coefficient match with n4_cramer: %s" % match
            else:
                match = len(sR.terms()) == len(R) and sorted(len(d.terms()) for d in sdets) == sorted(map(len, D))
                line += " | relabelled support; term counts match n4_cramer: %s" % match
            fails += not match
        print(line + " | %.0f s" % (time.time() - t0), flush=True)
    print("ALL PASS" if fails == 0 else "FAILURES: %d" % fails)
    sys.exit(1 if fails else 0)
