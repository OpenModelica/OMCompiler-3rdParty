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

#ifndef MOO_NLP_SCALING_H
#define MOO_NLP_SCALING_H

#include <base/fixed_vector.h>

namespace NLP {

/**
 * @brief interface for NLP variable and function scaling.
 *
 * set this in the NLP instance and these overloaded scaling routines
 * will be automatically called in the solver.
 */
class Scaling {
public:
    virtual ~Scaling() = default;

    virtual void inplace_scale_x(f64* x_unscaled) = 0;
    virtual void inplace_scale_g(f64* g_unscaled) = 0;

    virtual void unscale_x(const f64* x_scaled, f64* x_unscaled, int number_vars) = 0;
    virtual void unscale_g(const f64* g_scaled, f64* g_unscaled, int number_constraints) = 0;
    virtual void unscale_f(const f64* f_scaled, f64* f_unscaled) = 0;

    virtual void scale_f(const f64* f_unscaled, f64* f_scaled) = 0;
    virtual void scale_g(const f64* g_unscaled, f64* g_scaled, int number_constraints) = 0;
    virtual void scale_grad_f(const f64* grad_unscaled, f64* grad_scaled, int number_vars) = 0;
    virtual void scale_jac(const f64* jac_unscaled, f64* jac_scaled, int* i_row_jac, int* j_col_jac, int jac_nnz) = 0;
    virtual void scale_hes(const f64* hes_unscaled, f64* hes_scaled, int* i_row_hes, int* j_col_hes, int hes_nnz) = 0;
};

/**
 * @brief trivial no-op scaling: simply copies or does nothing.
 * this is useful when the NLP is already well-scaled.
 */
class NoScaling : public Scaling {
public:
    void inplace_scale_x(f64* x_unscaled) override {}

    void inplace_scale_g(f64* g_unscaled) override {}

    void unscale_x(const f64* x_scaled, f64* x_unscaled, int number_vars) override {
        std::memcpy(x_unscaled, x_scaled, number_vars * sizeof(f64));
    }

    void unscale_g(const f64* g_scaled, f64* g_unscaled, int number_constraints) override {
        std::memcpy(g_unscaled, g_scaled, number_constraints * sizeof(f64));
    }

    void unscale_f(const f64* f_scaled, f64* f_unscaled) override {
        *f_unscaled = *f_scaled;
    }

    void scale_f(const f64* f_unscaled, f64* f_scaled) override {
        *f_scaled = *f_unscaled;
    }

    void scale_g(const f64* g_unscaled, f64* g_scaled, int number_constraints) override {
        std::memcpy(g_scaled, g_unscaled, number_constraints * sizeof(f64));
    }

    void scale_grad_f(const f64* grad_unscaled, f64* grad_scaled, int number_vars) override {
        std::memcpy(grad_scaled, grad_unscaled, number_vars * sizeof(f64));
    }

    void scale_jac(const f64* jac_unscaled, f64* jac_scaled, int* i_row_jac, int* j_col_jac, int jac_nnz) override {
        std::memcpy(jac_scaled, jac_unscaled, jac_nnz * sizeof(f64));
    }

    void scale_hes(const f64* hes_unscaled, f64* hes_scaled, int* i_row_hes, int* j_col_hes, int hes_nnz) override {
        std::memcpy(hes_scaled, hes_unscaled, hes_nnz * sizeof(f64));
    }
};

class NominalScaling : public Scaling {
public:
    FixedVector<f64> x_scaling; // 1 / nominal
    FixedVector<f64> g_scaling; // 1 / nominal
    f64 f_scaling;              // 1 / nominal

    FixedVector<f64> grad_scaling;
    FixedVector<f64> jac_scaling;
    FixedVector<f64> hes_scaling;

    NominalScaling(FixedVector<f64>&& x_nominal, 
                   FixedVector<f64>&& g_nominal,
                   f64 f_nominal);

    void create_grad_scaling(int number_vars);
    void create_jac_scaling(int* i_row_jac, int* j_col_jac, int jac_nnz);
    void create_hes_scaling(int* i_row_hes, int* j_col_hes, int hes_nnz);

    void inplace_scale_x(f64* x_unscaled) override;
    void inplace_scale_g(f64* g_unscaled) override;

    void unscale_x(const f64* x_scaled, f64* x_unscaled, int number_vars) override;
    void unscale_g(const f64* g_scaled, f64* g_unscaled, int number_constraints) override;
    void unscale_f(const f64* f_scaled, f64* f_unscaled) override;

    void scale_f(const f64* f_unscaled, f64* f_scaled) override;
    void scale_g(const f64* g_unscaled, f64* g_scaled, int number_constraints) override;
    void scale_grad_f(const f64* grad_unscaled, f64* grad_scaled, int number_vars) override;
    void scale_jac(const f64* jac_unscaled, f64* jac_scaled,
                   int* i_row_jac, int* j_col_jac, int jac_nnz) override;
    void scale_hes(const f64* hes_unscaled, f64* hes_scaled,
                   int* i_row_hes, int* j_col_hes, int hes_nnz) override;
};

} // namespace NLP

#endif // MOO_NLP_SCALING_H
