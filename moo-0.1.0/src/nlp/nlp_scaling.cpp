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

#include <base/linalg.h>
#include <base/util.h>

#include "nlp_scaling.h"

// TODO: use BLAS here, and everywhere possible; these are literally straight vector ops

namespace NLP {

NominalScaling::NominalScaling(FixedVector<f64>&& x_nominal, 
                               FixedVector<f64>&& g_nominal,
                               f64 f_nominal)
    : x_scaling(std::move(x_nominal)),
      g_scaling(std::move(g_nominal)),
      f_scaling(1.0 / f_nominal)
    {
        for (size_t x = 0; x < x_scaling.size(); x++) {
            x_scaling[x] = 1.0 / x_scaling[x];
            assert(std::isfinite(x_scaling[x]) && "Variable vector scaling (x_scaling) contains inf or nan.");
        }

        for (size_t g = 0; g < g_scaling.size(); g++) {
            g_scaling[g] = 1.0 / g_scaling[g];
            assert(std::isfinite(g_scaling[g]) && "Constraint vector scaling (g_scaling) contains inf or nan.");
        }

        assert(std::isfinite(f_scaling) && "Objective function scaling (f_scaling) is inf or nan.");
    }

void NominalScaling::create_grad_scaling(int number_vars) {
    grad_scaling = FixedVector<f64>(number_vars);

    for (int i = 0; i < number_vars; i++) {
        grad_scaling[i] = f_scaling / x_scaling[i];
    }
}

void NominalScaling::create_jac_scaling(int* i_row_jac, int* j_col_jac, int jac_nnz) {
    jac_scaling = FixedVector<f64>(jac_nnz);

    for (int nz = 0; nz < jac_nnz; nz++) {
        jac_scaling[nz] = g_scaling[i_row_jac[nz]] / x_scaling[j_col_jac[nz]];
    }
}

void NominalScaling::create_hes_scaling(int* i_row_hes, int* j_col_hes, int hes_nnz) {
    hes_scaling = FixedVector<f64>(hes_nnz);

    for (int nz = 0; nz < hes_nnz; nz++) {
        hes_scaling[nz] = 1.0 / (x_scaling[i_row_hes[nz]] * x_scaling[j_col_hes[nz]]);
    }
}

// x_unscaled := x_nom * x_scaled
void NominalScaling::unscale_x(const f64* x_scaled, f64* x_unscaled, int number_vars) {
    Linalg::diagmat_vec(x_scaling.raw(), true, x_scaled, number_vars, x_unscaled);
}

// g_unscaled := g_nom * g_scaled
void NominalScaling::unscale_g(const f64* g_scaled, f64* g_unscaled, int number_constraints) {
    Linalg::diagmat_vec(g_scaling.raw(), true, g_scaled, number_constraints, g_unscaled);
}

// f_scaled := f_nom * f_unscaled 
void NominalScaling::unscale_f(const f64* f_scaled, f64* f_unscaled) {
    *f_unscaled = (*f_scaled) / f_scaling;
}

// x_scaled := x_nom^{-1} * x_unscaled
void NominalScaling::inplace_scale_x(f64* x_unscaled) {
    Linalg::diagmat_vec_inplace(x_scaling.raw(), false, x_unscaled, x_scaling.size());
}

// g_scaled := g_nom^{-1} * g_unscaled
void NominalScaling::inplace_scale_g(f64* g_unscaled) {
    Linalg::diagmat_vec_inplace(g_scaling.raw(), false, g_unscaled, g_scaling.size());
}

// f_scaled := f_nom^{-1} * f_unscaled 
void NominalScaling::scale_f(const f64* f_unscaled, f64* f_scaled) {
    *f_scaled = f_scaling * (*f_unscaled);
}

// g_scaled := g_nom^{-1} * g_unscaled 
void NominalScaling::scale_g(const f64* g_unscaled, f64* g_scaled, int number_constraints) {
    Linalg::diagmat_vec(g_scaling.raw(), false, g_unscaled, number_constraints, g_scaled);
}

// grad_scaled := f_nom^{-1} * grad_unscaled * x_nom = grad_scaling * grad_unscaled
void NominalScaling::scale_grad_f(const f64* grad_unscaled, f64* grad_scaled, int number_vars) {
    if (grad_scaling.empty()) {
        create_grad_scaling(number_vars);
    }

    Linalg::diagmat_vec(grad_scaling.raw(), false, grad_unscaled, number_vars, grad_scaled);
}

// jac_scaled := g_nom^{-1} * jac_unscaled * x_nom = jac_scaling * jac_unscaled
void NominalScaling::scale_jac(const f64* jac_unscaled, f64* jac_scaled,
                               int* i_row_jac, int* j_col_jac, int jac_nnz) {
    if (jac_scaling.empty()) {
        create_jac_scaling(i_row_jac, j_col_jac, jac_nnz);
    }

    Linalg::diagmat_vec(jac_scaling.raw(), false, jac_unscaled, jac_nnz, jac_scaled);
}

// hes_scaled := x_nom * H(sigma * f + lambda^T * g) * x_nom = hes_scaling * hes_unscaled
// while the f and g scaling is absorbed into lambda and sigma 
void NominalScaling::scale_hes(const f64* hes_unscaled, f64* hes_scaled,
                               int* i_row_hes, int* j_col_hes, int hes_nnz) {
    if (hes_scaling.empty()) {
        create_hes_scaling(i_row_hes, j_col_hes, hes_nnz);
    }

    Linalg::diagmat_vec(hes_scaling.raw(), false, hes_unscaled, hes_nnz, hes_scaled);
}

} // namespace NLP
