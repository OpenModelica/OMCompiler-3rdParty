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

#include <cassert>

#include <simulation/integrator/integrator_util.h>

namespace Simulation {

Jacobian::Jacobian(JacobianFormat jfmt,
                   int* i_row,
                   int* j_col,
                   int nnz)
    : jfmt(jfmt),
      i_row(i_row),
      j_col(j_col),
      nnz(nnz) {}

Jacobian Jacobian::dense() {
    return Jacobian(JacobianFormat::DENSE, nullptr, nullptr, 0);
}

Jacobian Jacobian::sparse(JacobianFormat sparse_fmt, int* i_row, int* j_col, int nnz) {
    assert(sparse_fmt != JacobianFormat::DENSE);
    return Jacobian(sparse_fmt, i_row, j_col, nnz);
}

} // namespace Simulation
