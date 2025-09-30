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

#ifndef MOO_INTEGRATOR_UTIL_H
#define MOO_INTEGRATOR_UTIL_H

#include <base/util.h>
#include <base/export.h>

namespace Simulation {

enum class JacobianFormat {
    DENSE,
    COO,
    CSC
};

class MOO_EXPORT Jacobian {
public:
    static Jacobian dense();
    static Jacobian sparse(JacobianFormat sparse_fmt, int* i_row, int* j_col, int nnz);

    JacobianFormat jfmt;
    int* i_row;
    int* j_col;
    int nnz;

private:
    Jacobian(JacobianFormat jfmt, int* i_row, int* j_col, int nnz);
};

using ODEFunction = std::function<void(const f64* x, const f64* u, const f64* p, f64 t, f64* f, void* user_data)>;
using JacobianFunction = std::function<void(const f64* x, const f64* u, const f64* p, f64 t, f64* dfdx, void* user_data)>;

} // namespace Simulation

#endif // MOO_INTEGRATOR_UTIL_H
