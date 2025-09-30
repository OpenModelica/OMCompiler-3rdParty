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

// TODO: make scaling accept NLP, only call after init x has been set?
//       => create more sophisticated scalings : `https://elib.dlr.de/93327/1/Performance_analysis_of_linear_and_nonlinear_techniques.pdf`

#include "nlp.h"

namespace NLP {

void NLP::solver_get_info(
    int& solver_number_vars,
    int& solver_number_constraints,
    int& solver_nnz_jac,
    int& solver_nnz_hes)
{
    // user queries
    get_sizes(number_vars, number_constraints);
    get_nnz(nnz_jac, nnz_hes);

    // allocate NLP buffers
    allocate_buffers();
    allocate_sparsity_buffers();

    // set solver data
    solver_number_vars        = number_vars;
    solver_number_constraints = number_constraints;
    solver_nnz_jac            = nnz_jac;
    solver_nnz_hes            = nnz_hes;

    // user query
    auto user_scaling = get_scaling();

    // create scaling
    if (user_scaling) {
        scaling = user_scaling;
    } else {
        Log::warning("Overloaded method get_scaling() returned nullptr - defaulting to NoScaling.");
        scaling = std::make_shared<NoScaling>();
    }
}

std::shared_ptr<Scaling> NLP::get_scaling()
{
    return std::make_shared<NoScaling>();
}

void NLP::allocate_buffers()
{
    // current iterates
    curr_x      = FixedVector<f64>(number_vars);
    curr_grad   = FixedVector<f64>(number_vars);
    curr_lambda = FixedVector<f64>(number_constraints);
    curr_g      = FixedVector<f64>(number_constraints);

    // problem bounds
    x_lb        = FixedVector<f64>(number_vars);
    x_ub        = FixedVector<f64>(number_vars);
    g_lb        = FixedVector<f64>(number_constraints);
    g_ub        = FixedVector<f64>(number_constraints);

    // bound multipliers (filled when optimal)
    z_lb        = FixedVector<f64>(number_vars);
    z_ub        = FixedVector<f64>(number_vars);
}

void NLP::allocate_sparsity_buffers()
{
    // COO sparsity patterns
    i_row_jac = FixedVector<int>(nnz_jac); // row COO of the Jacobian
    j_col_jac = FixedVector<int>(nnz_jac); // column COO of the Jacobian
    i_row_hes = FixedVector<int>(nnz_hes); // row COO of the Hessian
    j_col_hes = FixedVector<int>(nnz_hes); // column COO of the Hessian

    // values of sparse matrices
    curr_jac  = FixedVector<f64>(nnz_jac); // current NLP jacobian of the constraints
    curr_hes  = FixedVector<f64>(nnz_hes); // current NLP hessian of the lagrangian
}

void NLP::solver_get_bounds(
    f64* solver_x_lb,
    f64* solver_x_ub,
    f64* solver_g_lb,
    f64* solver_g_ub)
{
    // user query
    get_bounds(x_lb, x_ub, g_lb, g_ub);

    // copy
    x_lb.write_to(solver_x_lb);
    x_ub.write_to(solver_x_ub);
    g_lb.write_to(solver_g_lb);
    g_ub.write_to(solver_g_ub);

    // scale
    scaling->inplace_scale_x(solver_x_lb);
    scaling->inplace_scale_x(solver_x_ub);
    scaling->inplace_scale_g(solver_g_lb);
    scaling->inplace_scale_g(solver_g_ub);
}

void NLP::solver_get_initial_guess(
    bool init_x,
    f64* solver_x_init,
    bool init_lambda,
    f64* solver_lambda_init,
    bool init_z,
    f64* solver_z_lb_init,
    f64* solver_z_ub_init)
{
    // user query
    get_initial_guess(init_x, curr_x, init_lambda, curr_lambda, init_z, z_lb, z_ub);

    // scale + write to solver buffer
    if (init_x) {
        curr_x.write_to(solver_x_init);
        scaling->inplace_scale_x(solver_x_init);
    }
    if (init_lambda) {
        unscale_curr_lambda(solver_lambda_init);
    }
    if (init_z) {
        z_lb.write_to(solver_z_lb_init);
        z_ub.write_to(solver_z_ub_init);
        scaling->inplace_scale_x(solver_z_lb_init);
        scaling->inplace_scale_x(solver_z_ub_init);
    }
}

void NLP::solver_get_jac_sparsity(
    int* solver_i_row_jac,
    int* solver_j_col_jac)
{
    // user query
    get_jac_sparsity(i_row_jac, j_col_jac);

    i_row_jac.write_to(solver_i_row_jac);
    j_col_jac.write_to(solver_j_col_jac);
}

void NLP::solver_get_hes_sparsity(
    int* solver_i_row_hes,
    int* solver_j_col_hes)
{
    // user query
    get_hes_sparsity(i_row_hes, j_col_hes);

    i_row_hes.write_to(solver_i_row_hes);
    j_col_hes.write_to(solver_j_col_hes);
}

void NLP::solver_eval_f(
    bool new_x,
    const f64* solver_x,
    f64& solver_obj_value)
{
    update_unscale_curr_x(new_x, solver_x);

    // user query
    eval_f(new_x, curr_x, curr_obj);

    scaling->scale_f(&curr_obj, &solver_obj_value);
}

void NLP::solver_eval_grad_f(
    bool new_x,
    const f64* solver_x,
    f64* solver_grad_f)
{
    update_unscale_curr_x(new_x, solver_x);

    // user query
    eval_grad_f(new_x, curr_x, curr_grad);

    scaling->scale_grad_f(curr_grad.raw(), solver_grad_f, number_vars);
}

void NLP::solver_eval_g(
    bool new_x,
    const f64* solver_x,
    f64* solver_g)
{
    update_unscale_curr_x(new_x, solver_x);

    // user query
    eval_g(new_x, curr_x, curr_g);

    scaling->scale_g(curr_g.raw(), solver_g, number_constraints);
}

void NLP::solver_eval_jac(
    bool new_x,
    const f64* solver_x,
    f64* solver_jac)
{
    update_unscale_curr_x(new_x, solver_x);

    // user query
    eval_jac_g(new_x, curr_x, i_row_jac, j_col_jac, curr_jac);

    scaling->scale_jac(curr_jac.raw(), solver_jac, i_row_jac.raw(), j_col_jac.raw(), nnz_jac);
}

void NLP::solver_eval_hes(
    bool new_x,
    const f64* solver_x,
    bool new_lambda,
    const f64* solver_lambda,
    const f64 solver_obj_factor,
    f64* solver_hes)
{
    update_unscale_curr_x(new_x, solver_x);
    update_unscale_curr_lambda(new_lambda, solver_lambda);
    unscale_curr_sigma_f(&solver_obj_factor);

    // user query
    eval_hes(new_x, curr_x, new_lambda, curr_lambda, curr_sigma_f, i_row_hes, j_col_hes, curr_hes);

    scaling->scale_hes(curr_hes.raw(), solver_hes, i_row_hes.raw(), j_col_hes.raw(), nnz_hes);
}

void NLP::solver_finalize_solution(
    const f64  solver_obj_value,
    const f64* solver_x,
    const f64* solver_lambda,
    const f64* solver_z_lb,
    const f64* solver_z_ub)
{
    unscale_objective(&solver_obj_value);             // unscaled optimal objective
    update_unscale_curr_x(true, solver_x);            // unscaled optimal x
    update_unscale_curr_lambda(true, solver_lambda);  // unscaled optimal duals
    unscale_dual_bounds(solver_z_lb, solver_z_ub);    // unscaled optimal dual bound multipliers

    // user defined callback to extract info
    finalize_solution(curr_obj, curr_x, curr_lambda, z_lb, z_ub);
}

} // namespace NLP
