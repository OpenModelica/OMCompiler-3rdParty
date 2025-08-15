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

#include "adapter.h"

namespace IpoptSolver {

bool IpoptAdapter::get_nlp_info(
    Ipopt::Index& n,
    Ipopt::Index& m,
    Ipopt::Index& nnz_jac_g,
    Ipopt::Index& nnz_h_lag,
    IndexStyleEnum& index_style)
{
    nlp.solver_get_info(n, m, nnz_jac_g, nnz_h_lag);
    index_style = IndexStyleEnum::C_STYLE;
    return true;
};

bool IpoptAdapter::get_bounds_info(
    Ipopt::Index n,
    Ipopt::Number* x_l,
    Ipopt::Number* x_u,
    Ipopt::Index m,
    Ipopt::Number* g_l,
    Ipopt::Number* g_u)
{
    nlp.solver_get_bounds(x_l, x_u, g_l, g_u);
    return true;
};

bool IpoptAdapter::get_starting_point(
    Ipopt::Index n,
    bool init_x,
    Ipopt::Number* x,
    bool init_z,
    Ipopt::Number* z_L,
    Ipopt::Number* z_U,
    Ipopt::Index m,
    bool init_lambda,
    Ipopt::Number* lambda)
{
    nlp.solver_get_initial_guess(init_x, x, init_lambda, lambda, init_z, z_L, z_U);
    return true;
};

bool IpoptAdapter::get_scaling_parameters(
    Ipopt::Number& obj_scaling,
    bool& use_x_scaling,
    Ipopt::Index n,
    Ipopt::Number* x_scaling,
    bool& use_g_scaling,
    Ipopt::Index m,
    Ipopt::Number* g_scaling)
{
    use_x_scaling = false;
    use_g_scaling = false;
    return true;
};

bool IpoptAdapter::eval_f(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool new_x,
    Ipopt::Number& obj_value)
{
    nlp.solver_eval_f(new_x, x, obj_value);
    return true;
};

bool IpoptAdapter::eval_grad_f(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool new_x,
    Ipopt::Number* grad_f)
{
    nlp.solver_eval_grad_f(new_x, x, grad_f);
    return true;
};

bool IpoptAdapter::eval_g(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool new_x,
    Ipopt::Index m,
    Ipopt::Number* g)
{
    nlp.solver_eval_g(new_x, x, g);
    return true;
};

bool IpoptAdapter::eval_jac_g(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool new_x,
    Ipopt::Index m,
    Ipopt::Index nele_jac,
    Ipopt::Index* iRow,
    Ipopt::Index* jCol,
    Ipopt::Number* values)
{
    if (!values) {
        nlp.solver_get_jac_sparsity(iRow, jCol);
    }
    else {
        nlp.solver_eval_jac(new_x, x, values);
    }
    return true;
};

bool IpoptAdapter::eval_h(
    Ipopt::Index n,
    const Ipopt::Number* x,
    bool new_x,
    Ipopt::Number obj_factor,
    Ipopt::Index m,
    const Ipopt::Number* lambda,
    bool new_lambda,
    Ipopt::Index nele_hess,
    Ipopt::Index* iRow,
    Ipopt::Index* jCol,
    Ipopt::Number* values)
{
    if (!values) {
        nlp.solver_get_hes_sparsity(iRow, jCol);
    }
    else {
        nlp.solver_eval_hes(new_x, x, new_lambda, lambda, obj_factor, values);
    }
    return true;
};

void IpoptAdapter::finalize_solution(
    Ipopt::SolverReturn status,
    Ipopt::Index n,
    const Ipopt::Number* x,
    const Ipopt::Number* z_L,
    const Ipopt::Number* z_U,
    Ipopt::Index m,
    const Ipopt::Number* g,
    const Ipopt::Number* lambda,
    Ipopt::Number obj_value,
    const Ipopt::IpoptData* ip_data,
    Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    nlp.solver_finalize_solution(obj_value, x, lambda, z_L, z_U);
};

bool IpoptAdapter::intermediate_callback(
    Ipopt::AlgorithmMode mode,
    Ipopt::Index iter,
    Ipopt::Number obj_value,
    Ipopt::Number inf_pr,
    Ipopt::Number inf_du,
    Ipopt::Number mu, 
    Ipopt::Number d_norm,
    Ipopt::Number regularization_size,
    Ipopt::Number alpha_du,
    Ipopt::Number alpha_pr,
    Ipopt::Index ls_trials,
    const Ipopt::IpoptData* ip_data,
    Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    // TODO: to implement
    return true;
};


} // namespace IpoptSolver
