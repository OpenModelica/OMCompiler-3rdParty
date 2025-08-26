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

#ifndef MOO_LINALG_H
#define MOO_LINALG_H

#include <base/util.h>
#include <base/export.h>

// This Linear Algebra namespace implements several useful subroutines, which are (in part) straight BLAS wrappers, but also custom implementations

namespace Linalg {

enum class Norm {
    NORM_1,
    NORM_2,
    NORM_INF
};

inline std::string norm_to_string(Norm norm) {
    switch (norm) {
        case Norm::NORM_1:   return "l_1";
        case Norm::NORM_2:   return "l_2";
        case Norm::NORM_INF: return "l_inf";
        default:             return "Unknown";
    }
}

MOO_EXPORT f64 dot(int size, const f64* x, const f64* y);
MOO_EXPORT void matrix_vector(int size, char format, const f64* matrix, const f64* vector, f64* out);
MOO_EXPORT void square(int size, f64* x);
MOO_EXPORT void dsaxpy(int size, const f64* x, const f64* y, const f64* D, f64 beta, bool invD, f64* out);
MOO_EXPORT void dgmv(int size, const f64* x, const f64* y, const f64* D, f64 beta, bool invD, f64* out);
MOO_EXPORT void diagmat_vec(const f64* D, bool invD, const f64* x, int size, f64* out);
MOO_EXPORT void diagmat_vec_inplace(const f64* D, bool invD, f64* x, int size);

} // namespace Linalg

#endif  // MOO_LINALG_H
