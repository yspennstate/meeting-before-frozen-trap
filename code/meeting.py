"""Exact meeting-before-trapping probabilities for small reversible chains.

A reversible chain on {0, ..., n-1} is given by positive weights p (the stationary law is p / sum(p)) and symmetric
conductances c[i][j] >= 0; the walk jumps from i to j at rate c[i][j] / p[i].  X and Y run independently at rate 1,
Z is frozen, and X_0, Y_0, Z_0 are independent with law p / sum(p).  For a frozen Z = z, h_z(i, j) is the probability
that X and Y meet before either of them hits z, started from (i, j).  The first-step equations are solved on ordered
pairs in exact rational arithmetic.

F(p, c) = P(M_good <= M_bad)    (ties at time 0 count as meeting: X_0 = Y_0 = Z_0 is good)
strict=True gives P(M_good < M_bad) instead; discrete=True runs both walkers synchronously in discrete time with the
transition matrix c[i][j] / p[i] (rows must sum to one).

    python meeting.py      prints the examples quoted in the paper and checks them
"""
import sys
from fractions import Fraction as Fr


def solve(a, b):
    n = len(b)
    m = [list(map(Fr, row)) + [Fr(b[i])] for i, row in enumerate(a)]
    for col in range(n):
        piv = next(r for r in range(col, n) if m[r][col] != 0)
        m[col], m[piv] = m[piv], m[col]
        d = m[col][col]
        m[col] = [v / d for v in m[col]]
        for r in range(n):
            if r != col and m[r][col] != 0:
                f = m[r][col]
                m[r] = [u - f * v for u, v in zip(m[r], m[col])]
    return [m[i][n] for i in range(n)]


def probability(p, c, strict=False, discrete=False):
    n = len(p)
    q = [[Fr(c[i][j], 1) / p[i] if i != j else Fr(0) for j in range(n)] for i in range(n)]
    mu = [Fr(x) / sum(p) for x in p]
    total = Fr(0)
    for z in range(n):
        def boundary(i, j):
            if i == j:
                return Fr(0) if (strict and i == z) else Fr(1)
            if i == z or j == z:
                return Fr(0)
            return None
        states = [(i, j) for i in range(n) for j in range(n) if boundary(i, j) is None]
        ix = {s: k for k, s in enumerate(states)}
        a = [[Fr(0)] * len(states) for _ in states]
        b = [Fr(0)] * len(states)
        for row, (i, j) in enumerate(states):
            if discrete:
                if sum(q[i]) != 1 or sum(q[j]) != 1:
                    sys.exit("discrete mode needs a stochastic matrix")
                a[row][row] = Fr(1)
                moves = [(k, l, q[i][k] * q[j][l]) for k in range(n) for l in range(n)]
            else:
                a[row][row] = sum(q[i]) + sum(q[j])
                moves = [(k, j, q[i][k]) for k in range(n) if k != i] + [(i, l, q[j][l]) for l in range(n) if l != j]
            for k, l, rate in moves:
                if rate == 0:
                    continue
                v = boundary(k, l)
                if v is None:
                    a[row][ix[(k, l)]] -= rate
                else:
                    b[row] += rate * v
        h = solve(a, b) if states else []
        for i in range(n):
            for j in range(n):
                v = boundary(i, j)
                if v is None:
                    v = h[ix[(i, j)]]
                total += mu[z] * mu[i] * mu[j] * v
    return total


def bound(p):
    s = sum(p)
    return Fr(1, 3) + Fr(2, 3) * sum(Fr(x) ** 3 for x in p) / Fr(s) ** 3


def graph(n, edges, w=None):
    c = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(edges):
        c[i][j] = c[j][i] = 1 if w is None else w[k]
    return c


C4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
K4 = [(i, j) for i in range(4) for j in range(i + 1, 4)]

if __name__ == "__main__":
    fails = 0

    def check(label, got, want):
        global fails
        ok = got == want
        fails += not ok
        print("%-58s %-12s %s" % (label, got, "ok" if ok else "MISMATCH, expected %s" % want))

    check("uniform 4-cycle, F", probability([1] * 4, graph(4, C4)), Fr(3, 7))
    check("uniform K4, F", probability([1] * 4, graph(4, K4)), Fr(7, 16))
    check("bound 1/3 + (2/3) sum pi^3 at uniform pi, n = 4", bound([1] * 4), Fr(3, 8))
    check("uniform 4-cycle, strict event M_good < M_bad", probability([1] * 4, graph(4, C4), strict=True), Fr(41, 112))
    check("uniform K4, strict event M_good < M_bad", probability([1] * 4, graph(4, K4), strict=True), Fr(3, 8))
    half = [[Fr(1, 2) if (i, j) in C4 or (j, i) in C4 else 0 for j in range(4)] for i in range(4)]
    check("4-cycle, synchronous discrete time, F", probability([1] * 4, half, discrete=True), Fr(5, 16))
    two = probability([1, 3], [[0, 5], [5, 0]])
    check("two states, p = (1, 3): F equals the bound", two, bound([1, 3]))
    sys.exit(1 if fails else 0)
