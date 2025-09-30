// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2025 University of Applied Sciences and Arts
// Bielefeld, Faculty of Engineering and Mathematics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef MOO_FLGR_H
#define MOO_FLGR_H

#include <base/util.h>
#include <base/linalg.h>
#include <base/export.h>

class MOO_EXPORT fLGR {
public:
    fLGR() = delete;

    static constexpr size_t get_max_scheme() { return max_scheme; }

    static inline const f64* get_c(int scheme) { return &c[offset_linear(scheme)]; }
    static inline const f64* get_b(int scheme) { return &b[offset_linear(scheme)]; }
    static inline const f64* get_w(int scheme) { return &w[offset_linear(scheme)]; }

    static inline f64 get_c(int scheme, int index) { return c[offset_linear(scheme) + index]; }
    static inline f64 get_b(int scheme, int index) { return b[offset_linear(scheme) + index]; }
    static inline f64 get_w(int scheme, int index) { return w[offset_linear(scheme) + index]; }

    static inline const f64* get_c0(int scheme) { return &c0[offset_linear0(scheme)]; }
    static inline const f64* get_w0(int scheme) { return &w0[offset_linear0(scheme)]; }

    static inline f64 get_c0(int scheme, int index) { return c0[offset_linear0(scheme) + index]; }
    static inline f64 get_w0(int scheme, int index) { return w0[offset_linear0(scheme) + index]; }

    static inline const f64* get_D(int scheme)                   { return &D[offset_quadratic(scheme)]; }
    static inline const f64* get_D(int scheme, int row)          { return &D[offset_quadratic(scheme) + row * (scheme + 1)]; }
    static inline const f64  get_D(int scheme, int row, int col) { return  D[offset_quadratic(scheme) + row * (scheme + 1) + col]; }

    // integral 0 to 1 of some values
    static f64 integrate(int scheme, const f64* values);

    // multiply with the complete diff matrix: out := D * in
    static void diff_matrix_multiply(int scheme, const f64* in, f64* out);

    // multiply given differentiation matrix scheme with x_prev, (x_i0, ui0, xi1, ui1, ..., xnm, unm) u only for offset
    static void diff_matrix_multiply_block_strided(int scheme, int x_size, int xu_size, int fg_size, const f64* x_prev, const f64* x_new, f64* out);

    // evaluate the interpolating polynomial at some point T
    static f64 interpolate(int scheme, bool contains_zero, const f64* values, int stride, f64 interval_start, f64 interval_end, f64 T);

private:
    static constexpr size_t max_scheme = 100;

    // TODO: test if a lookup table is faster
    static inline size_t offset_linear(int scheme)     { return static_cast<size_t>((scheme - 1) * scheme / 2); }
    static inline size_t offset_quadratic(int scheme)  { return static_cast<size_t>(scheme * (scheme + 1) * (2 * scheme + 1) / 6); }
    static inline size_t offset_linear0(int scheme)    { return static_cast<size_t>(scheme * (scheme + 1) / 2); }

    // nodes, i.e. [c1, c2, ..., cm]
    static constexpr f64 c[] =
#include "../data/radauConstantsC.data"

    // nodes including 0, i.e. [0 = c0, c1, c2, ..., cm = 1]
    static constexpr f64 c0[] =
#include "../data/radauConstantsC0.data"

    // quadrature weights {}, {1}, ...
    static constexpr f64 b[] =
#include "../data/radauConstantsB.data"

    // differentiation matrices
    static constexpr f64 D[] =
#include "../data/radauConstantsD.data"

    // barycentric weights (all nodes [0 = c0, c1, ..., cm = 1])
    static constexpr f64 w[] =
#include "../data/radauConstantsW.data"

    // barycentric weights (only inner nodes [c1, ..., cm])
    static constexpr f64 w0[] =
#include "../data/radauConstantsW0.data"

};

#endif  // MOO_FLGR_H
