# SPDX-License-Identifier: LGPL-3.0-or-later
#
# This file is part of MOO - Modelica / Model Optimizer
# Copyright (C) 2025 University of Applied Sciences and Arts
# Bielefeld, Faculty of Engineering and Mathematics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from mpmath import *
from sympy import symbols, diff, expand, Poly, solve, N, re, im
import time


def genCppCode(clist, c0list, blist, Dlist, w0list, wlist, dps=250, mindps=40):
    mp.dps = mindps

    with open("radauConstantsC.txt", "w") as f:
        # c
        outStr = "{"
        for c in clist:
            for i in range(len(c)):
                outStr += str(c[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # c0
    with open("radauConstantsC0.txt", "w") as f:
        outStr = "{"
        for c0 in c0list:
            for i in range(len(c0)):
                outStr += str(c0[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # b
    with open("radauConstantsB.txt", "w") as f:
        outStr = "{"
        for b in blist:
            for i in range(len(b)):
                outStr += str(b[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # D1 matrix
    with open("radauConstantsD.txt", "w") as f:
        outStr = "{"
        for D1 in Dlist:
            for i in range(len(D1)):
                for j in range(len(D1[0])):
                    outStr += str(D1[i][j]) + ","
                outStr += "\n"
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # barycentric weights (incl. 0)
    with open("radauConstantsW0.txt", "w") as f:
        outStr = "{"
        for w0 in w0list:
            for i in range(len(w0)):
                outStr += str(w0[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # barycentric weights (excl. 0)
    with open("radauConstantsW.txt", "w") as f:
        outStr = "{"
        for w in wlist:
            for i in range(len(w)):
                outStr += str(w[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    return outStr


x = symbols("x")


def fixRoots(roots, dps=150):
    allRoots = []
    for r in roots:
        rDps = N(r, dps)
        if rDps.is_real:
            allRoots.append(rDps)
        elif abs(im(rDps)) <= 1e-150:
            allRoots.append(re(rDps))
        else:
            rDps = N(r, dps)
            allRoots.append(re(rDps))

    allRoots.sort()
    return allRoots


def lagrange(c, j, dps=150):
    p = mpf(1)
    for i in range(len(c)):
        if i != j:
            p *= (x - c[i]) / (c[j] - c[i])
    return p


def lagrangeToCoeffs(p, dps=150):
    return [mpf(coeff) for coeff in Poly(expand(p)).all_coeffs()[::-1]]


def integrate(poly, roots, i, dps=150):
    poly = lagrangeToCoeffs(poly, dps=dps)
    S = mpf(0)
    for j in range(len(roots)):
        S += (
            N(1, dps)
            / N(j + 1, dps)
            * poly[j]
            * (pow(N(1, dps), j + 1) - pow(N(0, dps), j + 1))
        )
    return S


def weight(c, i, dps=150):
    mp.dps = dps
    p = mpf(1)
    for j in range(len(c)):
        if i != j:
            p *= c[i] - c[j]
    return mpf(1) / p


def tridiagonal_eigenvalues(n, dps=150):
    mp.dps = dps

    alpha = mp.mpf(1)

    # tridiag entries
    a_0 = -mp.mpf(1) / mp.mpf(3)
    a_j = lambda j: -mp.mpf(1) / (
        (mp.mpf(2) * j + mp.mpf(1)) * (mp.mpf(2) * j + mp.mpf(3))
    )

    b_1 = mp.sqrt(mp.mpf(8) / (mp.mpf(9) * (mp.mpf(3) + alpha)))
    b_j = lambda j: mp.sqrt(
        (mp.mpf(4) * j * j * (j + mp.mpf(1)) * (j + mp.mpf(1)))
        / (
            (mp.mpf(2) * j)
            * (mp.mpf(2) * j + mp.mpf(1))
            * (mp.mpf(2) * j + mp.mpf(1))
            * (mp.mpf(2) * j + mp.mpf(2))
        )
    )

    # set matrix
    matrix = mp.zeros(n)
    for i in range(n):
        if i == 0:
            matrix[i, i] = a_0  # First diagonal element
            if n > 1:
                matrix[i, i + 1] = b_1
                matrix[i + 1, i] = b_1
        else:
            matrix[i, i] = a_j(i)  # Diagonal element
            if i < n - 1:  # Off-diagonal elements
                b_val = b_j(i + 1)
                matrix[i, i + 1] = b_val
                matrix[i + 1, i] = b_val

    # QR goes brummmmmmmm
    eigenvalues = mp.eigsy(matrix)[0]
    return eigenvalues


def generate(s, dps=150):
    mp.dps = dps

    # spectrum calculations
    roots = []

    # fLGR nodes
    if s > 1:
        roots = (sorted(tridiagonal_eigenvalues(s - 1, dps=dps)))
    roots.append(mp.mpf(1))

    # radau nodes
    c = [mpf(N(0.5, dps)) * (mpf(N(1, dps) + mpf(N(r, dps)))) for r in roots]

    # radau nodes with 0
    c0 = [mpf(0)] + c

    # quadrature weights
    b = []

    # barycentric weigths
    weights0 = [weight(c0, j) for j in range(len(c0))]
    weights = [weight(c, j) for j in range(len(c))]
    D1 = [[None for _ in range(len(c0))] for _ in range(len(c0))]  # diff matrix at c0

    # Barycentric Formulas from http://richard.baltensp.home.hefr.ch/Publications/3.pdf

    # eval d/dx L
    for i in range(s + 1):
        lagr = lagrange(c0, i)
        if i > 0:
            if s == 1:
                b = [mpf(1)]
            else:
                b.append(integrate(lagr, c0, s))
        for j in range(s + 1):
            if i != j:
                D1[i][j] = weights0[j] / (weights0[i] * (c0[i] - c0[j]))
        D1[i][i] = -sum(D1[i][j] for j in range(s + 1) if j != i)
    return c, c0, b, D1, weights0, weights

# ------------------------------------------------------------------------
# Access pattern for the flattened Radau constant arrays in C
#
# Assume:
#   - `scheme` is the number of collocation points
#   - c, c0, b, w, w0, D are all flattened 1D arrays in row-major order
#
# Offsets:
#   offset_linear(scheme)    = (scheme - 1) * scheme / 2                    for arrays of size 'scheme'
#   offset_linear0(scheme)   = scheme * (scheme + 1) / 2                    for arrays of size 'scheme + 1'
#   offset_quadratic(scheme) = scheme * (scheme + 1) * (2 * scheme + 1) / 6 for (scheme + 1) x (scheme + 1) D-matrix
#
# Access pattern:
#   c[scheme][i]     → c[offset_linear(scheme) + i]        (i = 0 .. scheme - 1)
#   c0[scheme][i]    → c0[offset_linear0(scheme) + i]      (i = 0 .. scheme)
#   b[scheme][i]     → b[offset_linear(scheme) + i]        (i = 0 .. scheme - 1)
#   w[scheme][i]     → w[offset_linear(scheme) + i]        (i = 0 .. scheme - 1)
#   w0[scheme][i]    → w0[offset_linear0(scheme) + i]      (i = 0 .. scheme)
#   D[scheme][i][j]  → D[offset_quadratic(scheme) + i * (scheme + 1) + j]  (i, j = 0 .. scheme)
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Why c0, w0, D start with zero entries:
#   We initialize c0list, w0list, and Dlist with zero-padding at index 0:
#     c0list = [[0.0]], w0list = [[0.0]], Dlist = [[[0.0]]]
#   This ensures that:
#     - offset_linear0(1) = 1 and accesses c0[1] at correct flat offset
#     - offset_quadratic(1) = 1 gives correct start for D[1]
#   In other words, by padding these arrays, we align scheme indices with their offsets
#   and make indexing cleaner and consistent starting from scheme = 1.
#   Otherwise we would need an additional -1 when indexing.
# ------------------------------------------------------------------------
clist, c0list, blist, Dlist, w0list, wlist = [], [[mpf(0.0)]], [], [[[mpf(0.0)]]], [[mpf(0.0)]], []

for m in range(1, 101):
    startTime = time.time()
    c, c0, b, D, w0, w = generate(m, 200)
    clist.append(c)
    c0list.append(c0)
    blist.append(b)
    Dlist.append(D)
    w0list.append(w0)
    wlist.append(w)
    print(f"{m} {time.time() - startTime}")

genCppCode(clist, c0list, blist, Dlist, w0list, wlist)
