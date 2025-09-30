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

#include "gdop.h"

// TODO: add timings for each component

namespace GDOP {

void GDOP::update(std::shared_ptr<const Mesh> new_mesh) {
    // update mesh
    mesh = new_mesh;

    // set the new mesh and update the callback buffers with new sizes
    problem.update_mesh(new_mesh);
}

void GDOP::create_acc_offset_xu(int off_x, int off_xu) {
    off_acc_xu = FixedField<int, 2>(mesh->intervals);
    int off = off_x;
    for (int i = 0; i < mesh->intervals; i++) {
        off_acc_xu[i] = FixedVector<int>(mesh->nodes[i]);
        for (int j = 0; j < mesh->nodes[i]; j++) {
            off_acc_xu[i][j] = off;
            off += off_xu;
        }
    }
}

void GDOP::create_acc_offset_fg(int off_fg) {
    off_acc_fg = FixedField<int, 2>(mesh->acc_nodes.size());
    for (int i = 0; i < mesh->intervals; i++) {
        off_acc_fg[i] = FixedVector<int>(mesh->acc_nodes[i]);
        for (int j = 0; j < mesh->nodes[i]; j++) {
            off_acc_fg[i][j] = mesh->acc_nodes[i][j] * off_fg;
        }
    }
}

// === overload ===
void GDOP::get_sizes(
    int& number_vars,
    int& number_constraints)
{
    off_x = problem.pc->x_size;
    off_u = problem.pc->u_size;
    off_p = problem.pc->p_size;
    off_xu = off_x + off_u;
    create_acc_offset_xu(off_x, off_xu);                           // variables  x_ij offset
    off_last_xu = off_acc_xu.back().back();                        // variables final grid point x_ij
    off_xu_total = off_last_xu + off_xu;                           // first parameter
    number_vars = off_xu_total + problem.pc->p_size;
    create_acc_offset_fg(problem.pc->fg_size);                     // constraint f_ij offset
    off_fg_total = mesh->node_count * problem.pc->fg_size;         // constraint r_0 offset
    number_constraints = problem.pc->r_size + off_fg_total;
}

void GDOP::set_scaling_factory(std::shared_ptr<ScalingFactory> factory) {
    scaling_factory = factory;
}

// === overload ===
std::shared_ptr<Scaling> GDOP::get_scaling() {
    if (scaling_factory) {
        return (*scaling_factory)(*this);
    } else {
        return nullptr;
    }
}

// === overload ===
void GDOP::get_bounds(
    FixedVector<f64>& x_lb,
    FixedVector<f64>& x_ub,
    FixedVector<f64>& g_lb,
    FixedVector<f64>& g_ub)
{
    // standard bounds, but checking for x0_fixed or xf_fixed
    for (int x_index = 0; x_index < off_x; x_index++) {
        x_lb[x_index] = problem.pc->x0_fixed[x_index] ? *problem.pc->x0_fixed[x_index] : problem.pc->x_bounds[x_index].lb;
        x_ub[x_index] = problem.pc->x0_fixed[x_index] ? *problem.pc->x0_fixed[x_index] : problem.pc->x_bounds[x_index].ub;
    }

    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (i == mesh->intervals - 1 && j == mesh->nodes[i] - 1) {
                for (int x_index = 0; x_index < off_x; x_index++) {
                    x_lb[off_acc_xu[i][j] + x_index] = problem.pc->xf_fixed[x_index] ? *problem.pc->xf_fixed[x_index] : problem.pc->x_bounds[x_index].lb;
                    x_ub[off_acc_xu[i][j] + x_index] = problem.pc->xf_fixed[x_index] ? *problem.pc->xf_fixed[x_index] : problem.pc->x_bounds[x_index].ub;
                }
            } 
            else {
                for (int x_index = 0; x_index < off_x; x_index++) {
                    x_lb[off_acc_xu[i][j] + x_index] = problem.pc->x_bounds[x_index].lb;
                    x_ub[off_acc_xu[i][j] + x_index] = problem.pc->x_bounds[x_index].ub;
                }
            }
            for (int u_index = 0; u_index < off_u; u_index++) {
                x_lb[off_acc_xu[i][j] + off_x + u_index] = problem.pc->u_bounds[u_index].lb;
                x_ub[off_acc_xu[i][j] + off_x + u_index] = problem.pc->u_bounds[u_index].ub;
            }
        }
    }

    for (int p_index = 0; p_index < off_p; p_index++) {
        x_lb[off_xu_total + p_index] = problem.pc->p_bounds[p_index].lb;
        x_ub[off_xu_total + p_index] = problem.pc->p_bounds[p_index].ub;
    }

    // standard constraint bounds
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem.pc->f_size; f_index++) {
                g_lb[off_acc_fg[i][j] + f_index] = 0;
                g_ub[off_acc_fg[i][j] + f_index] = 0;
            }
            for (int g_index = 0; g_index < problem.pc->g_size; g_index++) {
                g_lb[off_acc_fg[i][j] + problem.pc->f_size + g_index] = problem.pc->g_bounds[g_index].lb;
                g_ub[off_acc_fg[i][j] + problem.pc->f_size + g_index] = problem.pc->g_bounds[g_index].ub;
            }
        }
    }

    for (int r_index = 0; r_index < problem.pc->r_size; r_index++) {
        g_lb[off_fg_total + r_index] = problem.pc->r_bounds[r_index].lb;
        g_ub[off_fg_total + r_index] = problem.pc->r_bounds[r_index].ub;
    }
}

void GDOP::set_initial_guess(std::unique_ptr<PrimalDualTrajectory> initial_trajectory) {
    initial_guess = std::move(initial_trajectory);
}

// TODO: add flag: with apply_threshold_floor() or without - should be in the new yaml config file where strategies and stuff are given

// === overload ===
void GDOP::get_initial_guess(
      bool init_x,
      FixedVector<f64>& x_init,
      bool init_lambda,
      FixedVector<f64>& lambda_init,
      bool init_z,
      FixedVector<f64>& z_lb_init,
      FixedVector<f64>& z_ub_init)
{
    auto& initial_guess_primal         = initial_guess->primals;
    auto& initial_guess_costate        = initial_guess->costates;
    auto& initial_guess_lower_costates = initial_guess->lower_costates;
    auto& initial_guess_upper_costates = initial_guess->upper_costates;

    if (initial_guess_primal) {
        // check compatibility
        if (initial_guess_primal->inducing_mesh.get() != mesh.get()) {
            initial_guess_primal = std::make_unique<Trajectory>(initial_guess_primal->interpolate_onto_mesh(*mesh));
        }
        assert(check_time_compatibility(initial_guess_primal->t, {initial_guess_primal->x, initial_guess_primal->u}, *mesh)); // debug only

        flatten_trajectory_to_layout(*initial_guess_primal, x_init);
    }
    else {
        Log::error("No primal initial guess supplied in GDOP::init_starting_point().");
    }

    if (initial_guess_costate) {
        // check compatibility
        if (initial_guess_costate->inducing_mesh.get() != mesh.get()) {
            initial_guess_costate = std::make_unique<CostateTrajectory>(initial_guess_costate->interpolate_onto_mesh(*mesh));
        }
        assert(check_time_compatibility(initial_guess_costate->t, {initial_guess_costate->costates_f, initial_guess_costate->costates_g}, *mesh)); // debug only

        int index = 1; // ignore interpolated costates at t = 0
        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                for (int f_index = 0; f_index < problem.pc->f_size; f_index++) {
                    lambda_init[off_acc_fg[i][j] + f_index] = initial_guess_costate->costates_f[f_index][index];
                }
                for (int g_index = 0; g_index < problem.pc->g_size; g_index++) {
                    // apply threshold below 1e-10 to prevent oscillations
                    // lambda_init[off_acc_fg[i][j] + problem.pc->f_size + g_index] = apply_threshold_floor(value, 1e-10, 1e-12);

                    lambda_init[off_acc_fg[i][j] + problem.pc->f_size + g_index] = initial_guess_costate->costates_g[g_index][index];

                }
                index++;
            }
        }
        for (int r_index = 0; r_index < problem.pc->r_size; r_index++) {
            lambda_init[off_fg_total + r_index] = initial_guess_costate->costates_r[r_index];
        }

        transform_duals_costates(lambda_init, false);
    }

    if (initial_guess_lower_costates && initial_guess_upper_costates) {
        // check compatibility
        if (initial_guess_lower_costates->inducing_mesh.get() != mesh.get()) {
            initial_guess_lower_costates = std::make_unique<Trajectory>(initial_guess_lower_costates->interpolate_onto_mesh(*mesh));
        }

        if (initial_guess_upper_costates->inducing_mesh.get() != mesh.get()) {
            initial_guess_upper_costates = std::make_unique<Trajectory>(initial_guess_upper_costates->interpolate_onto_mesh(*mesh));
        }

        assert(check_time_compatibility(initial_guess_lower_costates->t, {initial_guess_lower_costates->x, initial_guess_lower_costates->u}, *mesh));
        assert(check_time_compatibility(initial_guess_lower_costates->t, {initial_guess_lower_costates->x, initial_guess_lower_costates->u}, *mesh));

        flatten_trajectory_to_layout(*initial_guess_lower_costates, z_lb_init);
        flatten_trajectory_to_layout(*initial_guess_upper_costates, z_ub_init);

        transform_duals_costates_bounds(z_lb_init, false);
        transform_duals_costates_bounds(z_ub_init, false);

        // apply threshold below 1e-10 to prevent oscillations
        // for (int idx = 0; idx < get_number_vars(); idx++) {
        //     z_lb_init[idx] = apply_threshold_floor(z_lb_init[idx], 1e-10, 1e-12);
        //     z_ub_init[idx] = apply_threshold_floor(z_ub_init[idx], 1e-10, 1e-12);
        // }
    }
}


void GDOP::init_jacobian_nonzeros(int& nnz_jac) {
    // stage 1: calculate nnz of blocks and number of collisions, where df_k / dx_k != 0. these are contained by default because of the D-Matrix
    int nnz_f = 0;
    int nnz_g = 0;
    int nnz_r = 0;
    int diagonal_collisions = 0;
    for (int f_index = 0; f_index < problem.pc->f_size; f_index++) {
        for (const auto& df_k_dx : problem.full->layout.f[f_index].jac.dx) {
            if (df_k_dx.col == f_index) {
                diagonal_collisions++;
            }
        }
        nnz_f += problem.full->layout.f[f_index].jac.nnz();
    }
    for (int g_index = 0; g_index < problem.pc->g_size; g_index++) {
            nnz_g += problem.full->layout.g[g_index].jac.nnz();
    }

    // nnz of block i can be calculated as m_i * ((m_i + 2) * #f + #g - coll(df_i, dx_i)), where m_i is the number of nodes on that interval
    off_acc_jac_fg = FixedVector<int>(mesh->intervals + 1);
    for (int i = 0; i < mesh->intervals; i++) {
        off_acc_jac_fg[i+1] = off_acc_jac_fg[i] + mesh->nodes[i] * ((mesh->nodes[i] + 1) * off_x + nnz_f + nnz_g - diagonal_collisions);
    }

    for (int r_index = 0; r_index < problem.pc->r_size; r_index++) {
        nnz_r += problem.boundary->layout.r[r_index].jac.nnz();
    }

    nnz_jac = off_acc_jac_fg.back() + nnz_r;

    // allocate memory for constant part
    const_der_jac = FixedVector<f64>(nnz_jac);
}

void GDOP::init_hessian_nonzeros(int& nnz_hes) {
    // takes O(nnz(A) + nnz(B) + ...+ nnz(H)) for creation of ** Maps and O(nnz(Hessian)) for creation of Hessian sparsity pattern

    // reset block sparsities
    hes_a_block = BlockSparsity::create_lower_triangular(problem.pc->x_size, BlockType::Exact);
    hes_b_block = BlockSparsity::create_lower_triangular(problem.pc->x_size + problem.pc->u_size, BlockType::Offset);
    hes_c_block = BlockSparsity::create_rectangular(problem.pc->x_size + problem.pc->u_size, problem.pc->x_size, BlockType::Exact);
    hes_d_block = BlockSparsity::create_lower_triangular(problem.pc->x_size + problem.pc->u_size, BlockType::Exact);
    hes_e_block = BlockSparsity::create_rectangular(problem.pc->p_size, problem.pc->x_size, BlockType::Exact);
    hes_f_block = BlockSparsity::create_rectangular(problem.pc->p_size, problem.pc->x_size + problem.pc->u_size, BlockType::RowOffset);
    hes_g_block = BlockSparsity::create_rectangular(problem.pc->p_size, problem.pc->x_size + problem.pc->u_size, BlockType::Exact);
    hes_h_block = BlockSparsity::create_lower_triangular(problem.pc->p_size, BlockType::Exact);

    // clear previous sparsity (may be reused)
    hes_A_set.clear(); hes_B_set.clear(); hes_C_set.clear(); hes_D_set.clear();
    hes_E_set.clear(); hes_F_set.clear(); hes_G_set.clear(); hes_H_set.clear();

    auto& boundary_hes = problem.boundary->layout.hes;
    auto& full_hes     = problem.full->layout.hes;
    auto& full_pp_hes  = problem.full->layout.pp_hes;

    // calculate IndexSet and nnz
    hes_A_set.insert_sparsity(boundary_hes.dx0_dx0,     0,     0);
    hes_C_set.insert_sparsity(boundary_hes.dxf_dx0,     0,     0);
    hes_C_set.insert_sparsity(boundary_hes.duf_dx0, off_x,     0);
    hes_D_set.insert_sparsity(boundary_hes.dxf_dxf,     0,     0);
    hes_D_set.insert_sparsity(boundary_hes.duf_dxf, off_x,     0);
    hes_D_set.insert_sparsity(boundary_hes.duf_duf, off_x, off_x);
    hes_E_set.insert_sparsity(boundary_hes.dp_dx0,      0,     0);
    hes_G_set.insert_sparsity(boundary_hes.dp_dxf,      0,     0);
    hes_G_set.insert_sparsity(boundary_hes.dp_duf,      0, off_x);
    hes_H_set.insert_sparsity(boundary_hes.dp_dp,       0,     0);

    hes_B_set.insert_sparsity(full_hes.dx_dx,           0,     0);
    hes_B_set.insert_sparsity(full_hes.du_dx,       off_x,     0);
    hes_B_set.insert_sparsity(full_hes.du_du,       off_x, off_x);
    hes_D_set.insert_sparsity(full_hes.dx_dx,           0,     0);
    hes_D_set.insert_sparsity(full_hes.du_dx,       off_x,     0);
    hes_D_set.insert_sparsity(full_hes.du_du,       off_x, off_x);
    hes_F_set.insert_sparsity(full_hes.dp_dx,           0,     0);
    hes_F_set.insert_sparsity(full_hes.dp_du,           0, off_x);
    hes_G_set.insert_sparsity(full_hes.dp_dx,           0,     0);
    hes_G_set.insert_sparsity(full_hes.dp_du,           0, off_x);
    hes_H_set.insert_sparsity(full_pp_hes.dp_dp,        0,     0);

    // calculate nnz from block sparsity
    nnz_hes = (hes_B_set.size() + hes_F_set.size()) * (mesh->node_count - 1)
             + hes_A_set.size() + hes_C_set.size() + hes_D_set.size() + hes_E_set.size() + hes_G_set.size() + hes_H_set.size();
}

// === overload ===
void GDOP::get_nnz(
    int& nnz_jac,
    int& nnz_hes)
{
    init_jacobian_nonzeros(nnz_jac);
    init_hessian_nonzeros(nnz_hes);
}

// === overload ===
void GDOP::get_jac_sparsity(
    FixedVector<int>& i_row_jac,
    FixedVector<int>& j_col_jac)
{
    // calculate the sparsity pattern i_row_jac, j_col_jac and the constant differentiation matrix part der_jac
    for (int i = 0; i < mesh->intervals; i++) {
        int nnz_index = off_acc_jac_fg[i]; // make local var: possible block parallelization
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem.pc->f_size; f_index++) {
                int eqn_index = off_acc_fg[i][j] + f_index;

                // dColl / dx for x_{i-1, m_{i-1}} base point states (k = -1, prev state)
                i_row_jac[nnz_index]     = eqn_index;
                j_col_jac[nnz_index]     = (i == 0 ? 0 : off_acc_xu[i - 1][mesh->nodes[i - 1] - 1]) + f_index;
                const_der_jac[nnz_index] = fLGR::get_D(mesh->nodes[i], j + 1, 0);
                nnz_index++;

                // dColl / dx for x_{i, j} collocation point states up to the collision
                // this means the derivative matrix linear combinations x_ik have k < j
                for (int k = 0; k < j; k++) {
                    i_row_jac[nnz_index]     = eqn_index;
                    j_col_jac[nnz_index]     = off_acc_xu[i][k] + f_index;
                    const_der_jac[nnz_index] = fLGR::get_D(mesh->nodes[i], j + 1, k + 1);
                    nnz_index++;
                }

                // df / dx
                int df_dx_counter = 0;
                std::vector<JacobianSparsity>* df_dx = &problem.full->layout.f[f_index].jac.dx;

                for (int x_elem = 0; x_elem < off_x; x_elem++) {
                    if (x_elem == f_index) {
                        // case for collocation block
                        i_row_jac[nnz_index]     = eqn_index;
                        j_col_jac[nnz_index]     = off_acc_xu[i][j] + f_index; // here df_dx and f_index collide!! (could also be written with dx_dx = f_index), which is nz element of the derivaitve matrix part
                        const_der_jac[nnz_index] = fLGR::get_D(mesh->nodes[i], j + 1, j + 1);

                        // handle the diagonal collision of the diagonal jacobian block
                        // this means the derivative matrix linear combinations x_ik have k == j and df_k / dx_k != 0
                        if (df_dx_counter < int_size(*df_dx) && (*df_dx)[df_dx_counter].col == x_elem) {
                            df_dx_counter++;
                        }

                        nnz_index++;
                    }
                    else if (df_dx_counter < int_size(*df_dx) && (*df_dx)[df_dx_counter].col == x_elem){
                        // no collision between collocation block and df / dx
                        i_row_jac[nnz_index] = eqn_index;
                        j_col_jac[nnz_index] = off_acc_xu[i][j] + x_elem;
                        nnz_index++;
                        df_dx_counter++;
                    }
                }

                // df / du
                for (auto& df_du : problem.full->layout.f[f_index].jac.du) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][j] + off_x + df_du.col;
                    nnz_index++;
                }

                // dColl / dx for x_{i, j} collocation point states after the collision / diagonal block in block jacobian
                // this means the derivative matrix linear combinations x_ik have k > j
                for (int k = j + 1; k < mesh->nodes[i]; k++) {
                    i_row_jac[nnz_index]     = eqn_index;
                    j_col_jac[nnz_index]     = off_acc_xu[i][k] + f_index;
                    const_der_jac[nnz_index] = fLGR::get_D(mesh->nodes[i], j + 1, k + 1);
                    nnz_index++;
                }

                // df / dp
                for (auto& df_dp : problem.full->layout.f[f_index].jac.dp) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_xu_total + df_dp.col;
                    nnz_index++;
                }
            }

            for (int g_index = 0; g_index < problem.pc->g_size; g_index++) {
                int eqn_index = off_acc_fg[i][j] + problem.pc->f_size + g_index;

                // dg / dx
                for (auto& dg_dx : problem.full->layout.g[g_index].jac.dx) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][j] + dg_dx.col;
                    nnz_index++;
                }

                // dg / du
                for (auto& dg_du : problem.full->layout.g[g_index].jac.du) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][j] + off_x + dg_du.col;
                    nnz_index++;
                }

                // dg / dp
                for (auto& dg_dp : problem.full->layout.g[g_index].jac.dp) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_xu_total + dg_dp.col;
                    nnz_index++;
                }
            }
        }
        assert(nnz_index == get_off_acc_jac_fg()[i + 1]);
    }

    int nnz_index = off_acc_jac_fg.back();
    for (int r_index = 0; r_index < problem.pc->r_size; r_index++) {
        int eqn_index = off_fg_total + r_index;

        // dr / dx0
        for (auto& dr_dx0 : problem.boundary->layout.r[r_index].jac.dx0) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = dr_dx0.col;
            nnz_index++;
        }

        // dr / dxf
        for (auto& dr_dxf : problem.boundary->layout.r[r_index].jac.dxf) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = off_last_xu + dr_dxf.col;
            nnz_index++;
        }

        // dr / duf
        for (auto& dr_duf : problem.boundary->layout.r[r_index].jac.duf) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = off_last_xu + off_x + dr_duf.col;
            nnz_index++;
        }

        // dr / dp
        for (auto& dr_dp: problem.boundary->layout.r[r_index].jac.dp) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = off_xu_total + dr_dp.col;
            nnz_index++;
        }
    }
    assert(nnz_index == get_nnz_jac());
}

// === overload ===
void GDOP::get_hes_sparsity(
    FixedVector<int>& i_row_hes,
    FixedVector<int>& j_col_hes)
{
    // allocate buffers for dual transformations

    // 1. get memory for Lagrange object factors
    if (problem.pc->has_lagrange) {
        lagrange_obj_factors = FixedField<f64, 2>(mesh->intervals);
        for (int i = 0; i < mesh->intervals; i++) {
            lagrange_obj_factors[i] = FixedVector<f64>(mesh->nodes[i]);
        }
    }

    // 2. get memory for transformed lambda
    transformed_lambda = FixedVector<f64>(get_number_constraints());


    // build int** (row, col) -> int index for all block structures
    // can be exact (no offset needed) or non exact (offset for full block or even rowwise needed)
    // also init sparsity pattern (i_row_hes, j_col_hes) for all exact blocks
    int hes_nnz_counter = 0;

    // A: exact
    for (auto& [row, col] : hes_A_set.set) {
        i_row_hes[hes_nnz_counter] = row; // x_0
        j_col_hes[hes_nnz_counter] = col; // x_0
        hes_a_block.insert(row, col, hes_nnz_counter++);
    }

    // B: non exact, thus local counter
    int block_b_nnz = 0;
    for (auto& [row, col] : hes_B_set.set) {
        hes_b_block.insert(row, col, block_b_nnz++);
    }
    hes_b_block.off_prev = hes_a_block.nnz;  // set size of A block as offset

    // init B hessian pattern O(node_count * nnz(L_{xu, xu} ∪ f_{xu, xu} ∪ g_{xu, xu})) - expensive, parallel execution should be possible
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (!(i == mesh->intervals - 1 && j == mesh->nodes[mesh->intervals - 1] - 1)) {
                for (auto [row, col] : hes_B_set.set) {
                    int xu_hes_index = hes_b_block.access(row, col, mesh->acc_nodes[i][j]);
                    i_row_hes[xu_hes_index] = off_acc_xu[i][j] + row; // xu_{ij}
                    j_col_hes[xu_hes_index] = off_acc_xu[i][j] + col; // xu_{ij}
                }
            }
        }
    }

    hes_nnz_counter += block_b_nnz * (mesh->node_count - 1);

    // C, D: exact with row dependence
    int c_index = 0;
    int d_index = 0;
    FixedVector<std::pair<int, int>> C_flat(hes_C_set.set.begin(), hes_C_set.set.end());
    FixedVector<std::pair<int, int>> D_flat(hes_D_set.set.begin(), hes_D_set.set.end());
    for (int xu_index = 0; xu_index < off_xu; xu_index++) {
        while (c_index < C_flat.int_size() && C_flat[c_index].first == xu_index) {
            i_row_hes[hes_nnz_counter] = off_last_xu + C_flat[c_index].first; // xu_{nm}
            j_col_hes[hes_nnz_counter] = C_flat[c_index].second;              // x_0
            hes_c_block.insert(C_flat[c_index].first, C_flat[c_index].second, hes_nnz_counter++);
            c_index++;
        }
        while (d_index < D_flat.int_size() && D_flat[d_index].first == xu_index) {
            i_row_hes[hes_nnz_counter] = off_last_xu + D_flat[d_index].first;  // xu_{nm}
            j_col_hes[hes_nnz_counter] = off_last_xu + D_flat[d_index].second; // xu_{nm}
            hes_d_block.insert(D_flat[d_index].first, D_flat[d_index].second, hes_nnz_counter++);
            d_index++;
        }
    }

    // E, F, G, H: partially exact (all except F) with row dependence
    int e_index = 0;
    int f_index = 0;
    int g_index = 0;
    int h_index = 0;
    FixedVector<std::pair<int, int>> E_flat(hes_E_set.set.begin(), hes_E_set.set.end());
    FixedVector<std::pair<int, int>> F_flat(hes_F_set.set.begin(), hes_F_set.set.end());
    FixedVector<std::pair<int, int>> G_flat(hes_G_set.set.begin(), hes_G_set.set.end());
    FixedVector<std::pair<int, int>> H_flat(hes_H_set.set.begin(), hes_H_set.set.end());
    for (int p_index = 0; p_index < problem.pc->p_size; p_index++) {
        while (e_index < E_flat.int_size() && E_flat[e_index].first == p_index) {
            i_row_hes[hes_nnz_counter] = off_xu_total + E_flat[e_index].first; // p
            j_col_hes[hes_nnz_counter] = E_flat[e_index].second;               // x_0
            hes_e_block.insert(E_flat[e_index].first, E_flat[e_index].second, hes_nnz_counter++);
            e_index++;
        }

        int row_f_nnz = 0;
        hes_f_block.row_offset_prev[p_index] = hes_nnz_counter; // E_{p_index, :} offset
        while (f_index < F_flat.int_size() && F_flat[f_index].first == p_index) {
            hes_f_block.insert(F_flat[f_index].first, F_flat[f_index].second, row_f_nnz++);
            f_index++;
        }
        /* F_{p_index, :} size -> offset for next F blocks */
        hes_f_block.row_size[p_index] = row_f_nnz;
        hes_nnz_counter += (mesh->node_count - 1) * row_f_nnz;

        while (g_index < G_flat.int_size() && G_flat[g_index].first == p_index) {
            i_row_hes[hes_nnz_counter] = off_xu_total + G_flat[g_index].first; // p
            j_col_hes[hes_nnz_counter] = off_last_xu + G_flat[g_index].second; // xu_{nm}
            hes_g_block.insert(G_flat[g_index].first, G_flat[g_index].second, hes_nnz_counter++);
            g_index++;
        }

        while (h_index < H_flat.int_size() && H_flat[h_index].first == p_index) {
            i_row_hes[hes_nnz_counter] = off_xu_total + G_flat[h_index].first;  // p
            j_col_hes[hes_nnz_counter] = off_xu_total + G_flat[h_index].second; // p
            hes_h_block.insert(H_flat[h_index].first, H_flat[h_index].second, hes_nnz_counter++);
            h_index++;
        }
    }

    // init F hessian pattern O(node_count * nnz(L_{p, xu} ∪ f_{p, xu} ∪ g_{p, xu})) - expensive, parallel execution should be possible
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (!(i == mesh->intervals - 1 && j == mesh->nodes[mesh->intervals - 1] - 1)) {
                for (auto& [row, col] : hes_F_set.set) {
                    int xu_hes_index = hes_f_block.access(row, col, mesh->acc_nodes[i][j]);
                    i_row_hes[xu_hes_index] = off_xu_total + row;     // p
                    j_col_hes[xu_hes_index] = off_acc_xu[i][j] + col; // xu_{ij}
                }
            }
        }
    }
}

/*
 * NLP function evaluations happen in two stages:
 * 
 * 1. Callback / continuous problem computation:
 *    - Fill input buffers (x, lambda, sigma).
 *    - Run callbacks (evaluation, Jacobian, or, Hessian) to compute all necessary information.
 * 
 * 2. NLP computation:
 *    - Evaluate the NLP function using pre-filled buffers and structures.
 *
 * Because callbacks are only in stage 1, and all data is ready before stage 2, 
 * the stage 2 computations are always thread-safe and parallelizable.
 *
 * Every NLP evaluation must call `check_new_x()` (and possibly lambda/sigma) to trigger callbacks if needed.
 *
 * Callback overview:
 *
 *            | Function      | `new_x` | `new_lambda` | Triggers Callback       | Needs Jacobian |
 *            | ------------- | ------- | ------------ | ----------------------- | -------------- |
 *            | `eval_f`      | +       | -            | `callback_evaluation()` | -              |
 *            | `eval_g`      | +       | -            | `callback_evaluation()` | -              |
 *            | `eval_grad_f` | +       | -            | `callback_jacobian()`   | +              |
 *            | `eval_jac_g`  | +       | -            | `callback_jacobian()`   | +              |
 *            | `eval_hes`    | +       | +            | `callback_hessian()`    | +              |
 */

// ========= virtuals in NLP =========

// === overload ===
// evaluate objective function
void GDOP::eval_f(
    bool new_x,
    const FixedVector<f64>& curr_x,
    f64& curr_obj)
{
    check_new_x(new_x);
    if (!evaluation_state.eval_f) {
        callback_evaluation(curr_x);
    }
    eval_f_internal(curr_obj);
}

// === overload ===
// evaluate constraints
void GDOP::eval_g(
    bool new_x,
    const FixedVector<f64>& curr_x,
    FixedVector<f64>& curr_g)
{
    check_new_x(new_x);
    if (!evaluation_state.eval_g) {
        callback_evaluation(curr_x);
    }
    eval_g_internal(curr_x, curr_g);
}

// === overload ===
// evaluate gradient of objective
void GDOP::eval_grad_f(
    bool new_x,
    const FixedVector<f64>& curr_x,
    FixedVector<f64>& curr_grad_f)
{
    check_new_x(new_x);
    if (!evaluation_state.grad_f) {
        callback_jacobian(curr_x);
    }
    eval_grad_f_internal(curr_grad_f);
}

// === overload ===
// evaluate Jacobian of constraints
void GDOP::eval_jac_g(
    bool new_x,
    const FixedVector<f64>& curr_x,
    const FixedVector<int>& i_row_jac,
    const FixedVector<int>& j_col_jac,
    FixedVector<f64>& curr_jac)
{
    check_new_x(new_x);
    if (!evaluation_state.jac_g) {
        callback_jacobian(curr_x);
    }
    eval_jac_g_internal(curr_jac);
}

// === overload ===
// evaluate Hessian of the Lagrangian
void GDOP::eval_hes(
    bool new_x,
    const FixedVector<f64>& curr_x,
    bool new_lambda,
    const FixedVector<f64>& curr_lambda,
    f64 curr_obj_factor,
    const FixedVector<int>& i_row_hes,
    const FixedVector<int>& j_col_hes,
    FixedVector<f64>& curr_hes)
{
    check_new_x(new_x);
    check_new_lambda(new_lambda);

    if (!evaluation_state.hes) {
        // ensure Jacobian is available (required for numerical Hessian)
        if (!evaluation_state.jac_g) {
            callback_jacobian(curr_x);
        }
        callback_hessian(curr_x, curr_lambda, curr_obj_factor);
    }
    eval_hes_internal(curr_hes);
}

// ========= callbacks and internal evaluation =========

// check if a new x was received; if so, reset evaluation state and update current x.
void GDOP::check_new_x(bool new_x) {
    evaluation_state.check_reset_x(new_x);
}

// similar check for lambda (dual variables).
void GDOP::check_new_lambda(bool new_lambda) {
    evaluation_state.check_reset_lambda(new_lambda);
}

void GDOP::callback_evaluation(const FixedVector<f64>& curr_x) {
    problem.full->callback_eval(get_x_xu(curr_x), get_x_p(curr_x));
    problem.boundary->callback_eval(get_x_x0(curr_x), get_x_xuf(curr_x), get_x_p(curr_x));
    evaluation_state.eval_f = true;
    evaluation_state.eval_g = true;
}

void GDOP::callback_jacobian(const FixedVector<f64>& curr_x) {
    problem.full->callback_jac(get_x_xu(curr_x), get_x_p(curr_x));
    problem.boundary->callback_jac(get_x_x0(curr_x), get_x_xuf(curr_x), get_x_p(curr_x));
    evaluation_state.grad_f = true;
    evaluation_state.jac_g = true;
}

// perform update of dual variables, such that callback can use the exact multiplier
void GDOP::update_curr_lambda_obj_factors(const FixedVector<f64>& curr_lambda, f64 curr_sigma_f) {
    transformed_lambda = curr_lambda; // copy curr_lambda

    for (int i = 0; i < mesh->intervals; i++) {
        f64 delta_t = mesh->delta_t[i];
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (problem.pc->has_lagrange) {
                lagrange_obj_factors[i][j] = curr_sigma_f * fLGR::get_b(mesh->nodes[i], j) * delta_t;
            }
            for (int f = 0; f < problem.pc->f_size; f++) {
                transformed_lambda[off_acc_fg[i][j] + f] *= -delta_t;
            }
        }
    }
}

void GDOP::callback_hessian(const FixedVector<f64> x, const FixedVector<f64>& curr_lambda, f64 curr_sigma_f) {
    update_curr_lambda_obj_factors(curr_lambda, curr_sigma_f);

    problem.full->callback_hes(get_x_xu(x), get_x_p(x), lagrange_obj_factors, get_lmbd_fg(transformed_lambda));
    problem.boundary->callback_hes(get_x_x0(x), get_x_xuf(x), get_x_p(x), curr_sigma_f, get_lmbd_r(transformed_lambda));
    evaluation_state.hes = true;
}

void GDOP::eval_f_internal(f64& curr_obj) {
    f64 mayer = 0;
    if (problem.pc->has_mayer) {
        mayer = problem.mr_eval_M();
    };

    f64 lagrange = 0;
    if (problem.pc->has_lagrange) {
        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                lagrange += mesh->delta_t[i] * fLGR::get_b(mesh->nodes[i], j) * problem.lfg_eval_L(i, j);
            }
        }
    }

    curr_obj = mayer + lagrange;
}

void GDOP::eval_grad_f_internal(FixedVector<f64>& curr_grad) {
    curr_grad.fill_zero();
    if (problem.pc->has_lagrange) {
        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                for (auto& dL_dx : problem.full->layout.L->jac.dx) {
                    curr_grad[off_acc_xu[i][j] + dL_dx.col] = mesh->delta_t[i] * fLGR::get_b(mesh->nodes[i], j) * problem.lfg_jac(dL_dx.buf_index, i, j);
                }
                for (auto& dL_du : problem.full->layout.L->jac.du) {
                    curr_grad[off_acc_xu[i][j] + off_x + dL_du.col] = mesh->delta_t[i] * fLGR::get_b(mesh->nodes[i], j) * problem.lfg_jac(dL_du.buf_index, i, j);
                }
                for (auto& dL_dp : problem.full->layout.L->jac.dp) {
                    curr_grad[off_xu_total + dL_dp.col] += mesh->delta_t[i] * fLGR::get_b(mesh->nodes[i], j) * problem.lfg_jac(dL_dp.buf_index, i, j);
                }
            }
        }
    }
    if (problem.pc->has_mayer) {
        for (auto& dM_dx0 : problem.boundary->layout.M->jac.dx0) {
            curr_grad[dM_dx0.col] += problem.mr_jac(dM_dx0.buf_index);
        }
        for (auto& dM_dxf : problem.boundary->layout.M->jac.dxf) {
            curr_grad[off_last_xu + dM_dxf.col] += problem.mr_jac(dM_dxf.buf_index);
        }
        for (auto& dM_duf : problem.boundary->layout.M->jac.duf) {
            curr_grad[off_last_xu + off_x + dM_duf.col] += problem.mr_jac(dM_duf.buf_index);
        }
        for (auto& dM_dp : problem.boundary->layout.M->jac.dp) {
            curr_grad[off_xu_total + dM_dp.col] += problem.mr_jac(dM_dp.buf_index);
        }
    }
};

void GDOP::eval_g_internal(const FixedVector<f64>& curr_x, FixedVector<f64>& curr_g) {
    curr_g.fill_zero();
    for (int i = 0; i < mesh->intervals; i++) {
        fLGR::diff_matrix_multiply_block_strided(mesh->nodes[i], off_x, off_xu, problem.pc->fg_size,
                                                       &curr_x[i == 0 ? 0 : off_acc_xu[i - 1][mesh->nodes[i - 1] - 1]],  // x_{i-1, m_{i-1}} base point states
                                                       &curr_x[off_acc_xu[i][0]],                                       // collocation point states
                                                       &curr_g[off_acc_fg[i][0]]);                                      // constraint start index 
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem.pc->f_size; f_index++) {
                curr_g[off_acc_fg[i][j] + f_index] -= mesh->delta_t[i] * problem.lfg_eval_f(f_index, i, j);
            }
            for (int g_index = 0; g_index < problem.pc->g_size; g_index++) {
                curr_g[off_acc_fg[i][j] + problem.pc->f_size + g_index] += problem.lfg_eval_g(g_index, i, j);
            }
        }
    }
    for (int r_index = 0; r_index < problem.pc->r_size; r_index++) {
        curr_g[off_fg_total + r_index] = problem.mr_eval_r(r_index);
    }
}

void GDOP::eval_jac_g_internal(FixedVector<f64>& curr_jac) {
    // copy constant part into curr_jac
    curr_jac = const_der_jac;

    for (int i = 0; i < mesh->intervals; i++) {
        int nnz_index = off_acc_jac_fg[i]; // make local var: possible block parallelization
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem.pc->f_size; f_index++) {
                // offset for all leading diagonal matrix blocks: d_{jk} * I at t_ik with j < k
                nnz_index += j + 1;

                int df_dx_counter = 0;
                std::vector<JacobianSparsity>& df_dx = problem.full->layout.f[f_index].jac.dx;

                for (int x_elem = 0; x_elem < off_x; x_elem++) {
                    if (x_elem == f_index) {
                        // case for collocation block
                        // handle the diagonal collision of the diagonal jacobian block
                        // this means the derivative matrix linear combinations x_ik have k == j and df_k / dx_k != 0
                        if (df_dx_counter < int_size(df_dx) && (df_dx)[df_dx_counter].col == x_elem) {
                            curr_jac[nnz_index] -= mesh->delta_t[i] * problem.lfg_jac((df_dx)[df_dx_counter].buf_index, i, j);
                            df_dx_counter++;
                        }
                        // even if df_k / dx_k == 0 => increment nnz from the collocation block nonzero
                        nnz_index++;
                    }
                    else if (df_dx_counter < int_size(df_dx) && (df_dx)[df_dx_counter].col == x_elem){
                        // no collision between collocation block and df / dx
                        curr_jac[nnz_index++] -= mesh->delta_t[i] * problem.lfg_jac((df_dx)[df_dx_counter].buf_index, i, j);
                        df_dx_counter++;
                    }
                }

                // df / du
                for (auto& df_du : problem.full->layout.f[f_index].jac.du) {
                    curr_jac[nnz_index++] = -mesh->delta_t[i] * problem.lfg_jac(df_du.buf_index, i, j);
                }

                // offset for all remaining diagonal matrix blocks: d_{jk} * I at t_ik with k > j
                 nnz_index += mesh->nodes[i] - j - 1;

                // df / dp
                for (auto& df_dp : problem.full->layout.f[f_index].jac.dp) {
                    curr_jac[nnz_index++] = -mesh->delta_t[i] * problem.lfg_jac(df_dp.buf_index, i, j);
                }
            }

            for (int g_index = 0; g_index < problem.pc->g_size; g_index++) {
                // dg / dx
                for (auto& dg_dx : problem.full->layout.g[g_index].jac.dx) {
                    curr_jac[nnz_index++] = problem.lfg_jac(dg_dx.buf_index, i, j);
                }

                // dg / du
                for (auto& dg_du : problem.full->layout.g[g_index].jac.du) {
                    curr_jac[nnz_index++] = problem.lfg_jac(dg_du.buf_index, i, j);
                }

                // dg / dp
                for (auto& dg_dp : problem.full->layout.g[g_index].jac.dp) {
                    curr_jac[nnz_index++] = problem.lfg_jac(dg_dp.buf_index, i, j);
                }
            }
        }
        assert(nnz_index == off_acc_jac_fg[i + 1]);
    }

    int nnz_index = off_acc_jac_fg.back();
    for (int r_index = 0; r_index < problem.pc->r_size; r_index++) {

        // dr / dx0
        for (auto& dr_dx0 : problem.boundary->layout.r[r_index].jac.dx0) {
            curr_jac[nnz_index++] = problem.mr_jac(dr_dx0.buf_index);
        }

        // dr / dxf
        for (auto& dr_dxf : problem.boundary->layout.r[r_index].jac.dxf) {
            curr_jac[nnz_index++] = problem.mr_jac(dr_dxf.buf_index);
        }

        // dr / duf
        for (auto& dr_duf : problem.boundary->layout.r[r_index].jac.duf) {
            curr_jac[nnz_index++] = problem.mr_jac(dr_duf.buf_index);
        }

        // dr / dp
        for (auto& dr_dp: problem.boundary->layout.r[r_index].jac.dp) {
            curr_jac[nnz_index++] = problem.mr_jac(dr_dp.buf_index);
        }
    }
    assert(nnz_index == get_nnz_jac());
};

void GDOP::eval_hes_internal(FixedVector<f64>& curr_hes) {
    curr_hes.fill_zero();

    for (int i = 0; i < mesh->intervals; i++) {
        // insert buffer for parallel execution here / pass the buffer in updateHessian(f64*) calls
        for (int j = 0; j < mesh->nodes[i]; j++) {
            // make sure to be in the right ptr_map region, B and F are only valid for i,j != n,m
            //                                              D and G are only valid for i,j == n,m
            const BlockSparsity* ptr_map_xu_xu;
            const BlockSparsity* ptr_map_p_xu;
            if (!(i == mesh->intervals - 1 && j == mesh->nodes[mesh->intervals - 1] - 1)) {
                ptr_map_xu_xu = &hes_b_block;
                ptr_map_p_xu  = &hes_f_block;
            }
            else {
                ptr_map_xu_xu = &hes_d_block;
                ptr_map_p_xu  = &hes_g_block;
            }
            update_hessian_lfg(problem.full->layout.hes, i, j, ptr_map_xu_xu, ptr_map_p_xu, curr_hes);
        }
    }
    update_parameter_hessian_lfg(problem.full->layout.pp_hes, curr_hes);
    update_hessian_mr(problem.boundary->layout.hes, curr_hes);
}

void GDOP::update_hessian_lfg(const HessianLFG& hes, const int i, const int j,
                                        const BlockSparsity* ptr_map_xu_xu, const BlockSparsity* ptr_map_p_xu,
                                        FixedVector<f64>& curr_hes) {
    const int block_count = mesh->acc_nodes[i][j];
    for (const auto& dx_dx : hes.dx_dx) {
        curr_hes[ptr_map_xu_xu->access(dx_dx.row, dx_dx.col, block_count)] += problem.lfg_hes(dx_dx.buf_index, i, j);
    }
    for (const auto& du_dx : hes.du_dx) {
        curr_hes[ptr_map_xu_xu->access(off_x + du_dx.row, du_dx.col, block_count)] += problem.lfg_hes(du_dx.buf_index, i, j);
    }
    for (const auto& du_du : hes.du_du) {
        curr_hes[ptr_map_xu_xu->access(off_x + du_du.row, off_x + du_du.col, block_count)] += problem.lfg_hes(du_du.buf_index, i, j);
    }
    for (const auto& dp_dx : hes.dp_dx) {
        curr_hes[ptr_map_p_xu->access(dp_dx.row, dp_dx.col, block_count)] += problem.lfg_hes(dp_dx.buf_index, i, j);
    }
    for (const auto& dp_du : hes.dp_du) {
        curr_hes[ptr_map_p_xu->access(dp_du.row, off_x + dp_du.col, block_count)] += problem.lfg_hes(dp_du.buf_index, i, j);
    }
}

void GDOP::update_parameter_hessian_lfg(const ParameterHessian& pp_hes, FixedVector<f64>& curr_hes) {
    for (const auto& dp_dp : pp_hes.dp_dp) {
        curr_hes[hes_h_block.access(dp_dp.row, dp_dp.col)] += problem.lfg_pp_hes(dp_dp.buf_index);
    }
}

void GDOP::update_hessian_mr(const HessianMR& hes, FixedVector<f64>& curr_hes) {
    for (const auto& dx0_dx0 : hes.dx0_dx0) {
        curr_hes[hes_a_block.access(dx0_dx0.row, dx0_dx0.col)] += problem.mr_hes(dx0_dx0.buf_index);
    }
    for (const auto& dxf_dx0 : hes.dxf_dx0) {
        curr_hes[hes_c_block.access(dxf_dx0.row, dxf_dx0.col)] += problem.mr_hes(dxf_dx0.buf_index);
    }
    for (const auto& dxf_dxf : hes.dxf_dxf) {
        curr_hes[hes_d_block.access(dxf_dxf.row, dxf_dxf.col)] += problem.mr_hes(dxf_dxf.buf_index);
    }
    for (const auto& duf_dx0 : hes.duf_dx0) {
        curr_hes[hes_c_block.access(off_x + duf_dx0.row, duf_dx0.col)] += problem.mr_hes(duf_dx0.buf_index);
    }
    for (const auto& duf_dxf : hes.duf_dxf) {
        curr_hes[hes_d_block.access(off_x + duf_dxf.row, duf_dxf.col)] += problem.mr_hes(duf_dxf.buf_index);
    }
    for (const auto& duf_duf : hes.duf_duf) {
        curr_hes[hes_d_block.access(off_x + duf_duf.row, off_x + duf_duf.col)] += problem.mr_hes(duf_duf.buf_index);
    }
    for (const auto& dp_dx0 : hes.dp_dx0) {
        curr_hes[hes_e_block.access(dp_dx0.row, dp_dx0.col)]  += problem.mr_hes(dp_dx0.buf_index);
    }
    for (const auto& dp_dxf : hes.dp_dxf) {
        curr_hes[hes_g_block.access(dp_dxf.row, dp_dxf.col)]  += problem.mr_hes(dp_dxf.buf_index);
    }
    for (const auto& dp_duf : hes.dp_duf) {
        curr_hes[hes_g_block.access(dp_duf.row, off_x + dp_duf.col)]  += problem.mr_hes(dp_duf.buf_index);
    }
    for (const auto& dp_dp : hes.dp_dp) {
        curr_hes[hes_h_block.access(dp_dp.row, dp_dp.col)]   += problem.mr_hes(dp_dp.buf_index);
    }
}


// === Optimal Solution Retrieval and Costate Estimations ===

/**
 * Dual Transformation in Direct Collocation for Dynamic Optimization
 *
 * When solving a dynamic optimization problem using direct collocation with flipped Legendre-Gauss-Radau (fLGR)
 * quadrature, the optimizer returns Karush-Kuhn-Tucker (KKT) multipliers.
 * These discrete multipliers represent the sensitivity of the objective function to changes in the
 * constraints at specific collocation points. To obtain **smooth, continuous-time trajectories**, and also
 * continuous costate estimates, and other dual variables, these discrete multipliers must be transformed.
 *
 * IPOPT returns KKT multipliers for:
 * - **ODE constraints** ($\lambda_f$): Duals associated with the system dynamics.
 * - **Path constraints** ($\lambda_g$): Duals for inequality constraints that apply over the trajectory.
 * - **Boundary constraints** ($\lambda_r$): Duals for inequality constraints at the initial/final time.
 * - **Bound constraints on x, u, p** ($z_{lb}, z_{ub}$): Multipliers for lower and upper bounds on states, controls, and parameters.
 *
 * The transformation from discrete NLP multipliers to continuous-time trajectories for each
 * collocation point $t_{ij}$ with quadrature weight $b_j$ is as follows:
 *
 * $\lambda_f(t_{ij})  = -\lambda_{f\_NLP}(t_{ij}) / b_j$     // Dynamics: Requires a sign flip and scaling by the quadrature weight.
 * $\lambda_g(t_{ij})  = -\lambda_{g\_NLP}(t_{ij}) / b_j$     // Path constraints: Requires scaling by the quadrature weight; sign convention may vary based on problem formulation.
 * $\zeta_{xu}(t_{ij}) = -z_{NLP}(t_{ij}) / b_j$              // Bound multipliers ($\zeta$): Requires scaling by the quadrature weight; sign convention may vary based on problem formulation.
 *
 * Note: Multipliers for static parameters `p` and boundary constraints `r` are generally not
 * scaled by quadrature weights as they are static and not distributed across collocation points.
 *
 * **Important Considerations:**
 * - **Sign Conventions:** Be vigilant about sign conventions, as they can differ between optimization solvers and theoretical definitions.
 *          The negative sign in the transformation for $\lambda_f$ is crucial for obtaining the correct costates.
 * - **Numerical Noise:** Very small bound multipliers (e.g., below $1 \times 10^{-10}$ after scaling) likely exhibit numerical noise. These often arise
 *          from the interior-point nature of the solver and loose thresholds for constraint activity. Such values typically do not exhibit the
 *          expected polynomial structure of collocation and can often be safely ignored or thresholded (e.g., set to fixed eps when below a cutoff).
 * - **Smoothness of Bound Multipliers:** Bound multipliers $\zeta_x(t)$ and $\zeta_u(t)$ may appear smooth (exhibiting the expected piecewise polynomial structure)
 *          if the corresponding bound constraint is active over a full collocation interval.
 */

void GDOP::flatten_trajectory_to_layout(const Trajectory& trajectory, FixedVector<f64>& flat_buffer) {
    for (int x_index = 0; x_index < off_x; x_index++) {
        flat_buffer[x_index] = trajectory.x[x_index][0];
    }

    int index = 1;
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int x_index = 0; x_index < off_x; x_index++) {
                flat_buffer[off_acc_xu[i][j] + x_index] = trajectory.x[x_index][index];
            }
            for (int u_index = 0; u_index < off_u; u_index++) {
                flat_buffer[off_acc_xu[i][j] + off_x + u_index] = trajectory.u[u_index][index];
            }
            index++;
        }
    }
    for (int p_index = 0; p_index < off_p; p_index++) {
        flat_buffer[off_xu_total + p_index] = trajectory.p[p_index];
    }
}

 /**
 * @brief Finalizes and returns the optimal primal trajectories (states and controls).
 *
 * This function extracts the optimal primal variables (state `x` and control `u`)
 * from the `curr_x` member, which holds the solution returned by IPOPT, and
 * organizes them into a `Trajectory` object. It populates the time points,
 * state trajectories, control trajectories, and static parameters.
 *
 * The initial state `x(0)` and control `u(0)` are handled specifically, with `u(0)`
 * being interpolated if necessary. Subsequent points correspond to the collocation
 * nodes.
 *
 * @note The Trajectories are unscaled, as they are accessible from the NLP and not from the Solver.
 *
 * @return A `std::unique_ptr<Trajectory>` containing the time, state, and control trajectories.
 */
std::unique_ptr<Trajectory> GDOP::finalize_optimal_primals(const FixedVector<f64>& opt_x) {
    auto optimal_primals = std::make_unique<Trajectory>();

    optimal_primals->t.reserve(mesh->node_count + 1);
    optimal_primals->x.resize(off_x);
    optimal_primals->u.resize(off_u);
    optimal_primals->inducing_mesh = mesh->shared_from_this();
    optimal_primals->interpolation = InterpolationMethod::POLYNOMIAL;

    for (auto& v : optimal_primals->x) { v.reserve(mesh->node_count + 1); }
    for (auto& v : optimal_primals->u) { v.reserve(mesh->node_count + 1); }

    for (int x_index = 0; x_index < off_x; x_index++) {
        optimal_primals->x[x_index].push_back(opt_x[x_index]);
    }

    for (int u_index = 0; u_index < off_u; u_index++) {
        f64 u0 = fLGR::interpolate(mesh->nodes[0], false, &opt_x[2 * off_x + u_index], off_xu, mesh->t[0][0], mesh->grid[1], 0.0);
        optimal_primals->u[u_index].push_back(u0);
    }

    optimal_primals->t.push_back(0.0);

    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int x_index = 0; x_index < off_x; x_index++) {
                optimal_primals->x[x_index].push_back(opt_x[off_acc_xu[i][j] + x_index]);
            }

            for (int u_index = 0; u_index < off_u; u_index++) {
                optimal_primals->u[u_index].push_back(opt_x[off_acc_xu[i][j] + off_x + u_index]);
            }

            optimal_primals->t.push_back(mesh->t[i][j]);
        }
    }

    if (off_p > 0) {
        optimal_primals->p = std::vector(get_x_p(opt_x), opt_x.end());
    }

    return optimal_primals;
}

/**
 * @brief Transforms dual multipliers (from IPOPT) to continuous-time costates, or vice-versa.
 *
 * This function applies the necessary scaling and sign flips to convert NLP KKT multipliers
 * for ODE constraints ($\lambda_f$) and path constraints ($\lambda_g$) into their
 * continuous-time costate equivalents, or to convert costates back into NLP duals.
 * The transformation involves dividing (or multiplying) by the negative of the
 * collocation quadrature weights ($b_j$).
 *
 * The transformation is applied in-place to the provided `lambda` vector.
 *
 * @param lambda A `FixedVector<f64>` containing the dual multipliers or costates to be transformed.
 * @param to_costate A boolean flag. If `true`, transforms NLP duals to costates. If `false`, transforms costates to NLP duals.
 */
void GDOP::transform_duals_costates(FixedVector<f64>& lambda, bool to_costate) {
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int fg_index = 0; fg_index < problem.pc->fg_size; fg_index++) {
                if (to_costate) {
                    // NLP dual lambda -> costates lambda
                    lambda[off_acc_fg[i][j] + fg_index] /= -fLGR::get_b(mesh->nodes[i], j);
                }
                else {
                    // costates lambda -> NLP dual lambda
                    lambda[off_acc_fg[i][j] + fg_index] *= -fLGR::get_b(mesh->nodes[i], j);
                }
            }
        }
    }
}

/**
 * @brief Finalizes and returns the optimal costate trajectories.
 *
 * This function processes the KKT multipliers for ODE path constraints, and also boundary
 * constraints (`curr_lambda`) obtained from IPOPT, transforms them into continuous-time
 * costates using `transform_duals_costates`, and then organizes them into a
 * `CostateTrajectory` object.
 *
 * The costates at $t=0$ are obtained through interpolation. Costates for boundary
 * constraints ($\lambda_r$) are also included and just given by the KKT multiplier.
 *
 * @return A `std::unique_ptr<CostateTrajectory>` containing the time, costates for dynamics,
 * costates for path constraints, and costates for boundary constraints.
 */
std::unique_ptr<CostateTrajectory> GDOP::finalize_optimal_costates(const FixedVector<f64>& opt_lambda) {
    auto optimal_costates = std::make_unique<CostateTrajectory>();

    const int f_size = problem.pc->f_size;
    const int g_size = problem.pc->g_size;

    // transform curr_lambda from NLP duals -> costates
    FixedVector<f64> costates(opt_lambda);
    transform_duals_costates(costates, true);

    optimal_costates->t.reserve(mesh->node_count + 1);
    optimal_costates->costates_f.resize(f_size);
    optimal_costates->costates_g.resize(g_size);
    optimal_costates->inducing_mesh = mesh->shared_from_this();
    optimal_costates->interpolation = InterpolationMethod::POLYNOMIAL;

    for (auto& v : optimal_costates->costates_f) { v.reserve(mesh->node_count + 1); }
    for (auto& v : optimal_costates->costates_g) { v.reserve(mesh->node_count + 1); }

    // interpolate lambda at t = 0
    const int inp_stride = f_size + g_size;

    for (int f_index = 0; f_index < f_size; f_index++) {
        f64 lambda_f_0 = fLGR::interpolate(mesh->nodes[0], false, &costates[f_index], inp_stride, mesh->t[0][0], mesh->grid[1], 0.0);
        optimal_costates->costates_f[f_index].push_back(lambda_f_0);
    }

    for (int g_index = 0; g_index < g_size; g_index++) {
        f64 lambda_g_0 = fLGR::interpolate(mesh->nodes[0], false, &costates[f_size + g_index], inp_stride, mesh->t[0][0], mesh->grid[1], 0.0);
        optimal_costates->costates_g[g_index].push_back(lambda_g_0);
    }

    optimal_costates->t.push_back(0.0);

    // use exact values for the others
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < f_size; f_index++) {
                optimal_costates->costates_f[f_index].push_back(costates[off_acc_fg[i][j] + f_index]);
            }

            for (int g_index = 0; g_index < g_size; g_index++) {
                optimal_costates->costates_g[g_index].push_back(costates[off_acc_fg[i][j] + f_size + g_index]);
            }

            optimal_costates->t.push_back(mesh->t[i][j]);
        }
    }

    optimal_costates->costates_r.resize(problem.pc->r_size);

    for (int r_index = 0; r_index < problem.pc->r_size; r_index++) {
        optimal_costates->costates_r[r_index] = costates[off_fg_total + r_index];
    }

    return optimal_costates;
}


/**
 * @brief Transforms dual multipliers for variable bounds to their continuous-time costate equivalents, or vice-versa.
 *
 * This function performs the transformation of IPOPT's KKT multipliers for lower and upper
 * bounds on states (`x`) and controls (`u`) into their continuous-time dual counterparts
 * ($\zeta_x(t)$ and $\zeta_u(t)$). The transformation involves dividing (or multiplying)
 * by the negative of the collocation quadrature weights ($b_j$).
 *
 * The transformation is applied in-place to the provided `zeta` vector.
 *
 * @param zeta A `FixedVector<f64>` containing the bound multipliers or their transformed costate equivalents.
 * @param to_costate A boolean flag. If `true`, transforms NLP duals to costates. If `false`, transforms costates to NLP duals.
 *
 * @attention The sign convention for non-ODE costates is not clear right now.
 */
void GDOP::transform_duals_costates_bounds(FixedVector<f64>& zeta, bool to_costate) {
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int xu_index = 0; xu_index < off_xu; xu_index++) {
                if (to_costate) {
                    // NLP dual lambda -> costates lambda
                    zeta[off_acc_xu[i][j] + xu_index] /= -fLGR::get_b(mesh->nodes[i], j);
                }
                else {
                    // costates lambda -> NLP dual lambda
                    zeta[off_acc_xu[i][j] + xu_index] *= -fLGR::get_b(mesh->nodes[i], j);
                }
            }
        }
    }
}

/**
 * @brief Finalizes and returns the optimal dual trajectories for lower and upper variable bounds.
 *
 * This function takes the KKT multipliers for lower (`z_lb`) and upper (`z_ub`) bounds
 * on states and controls from IPOPT, transforms them into continuous-time bound duals
 * using `transform_duals_costates_bounds`, and then organizes them into a pair of
 * `Trajectory` objects (one for lower bounds, one for upper bounds).
 *
 * These bound duals, denoted as $\zeta_x(t)$ and $\zeta_u(t)$, represent the
 * sensitivity of the optimal cost to changes in the bounds.
 *
 * @return A `std::pair<std::unique_ptr<Trajectory>, std::unique_ptr<Trajectory>>`
 * where the first element is the trajectory of lower bound duals and the
 * second is the trajectory of upper bound duals.
 *
 * @note Very small bound multipliers (e.g., below 1e-10) may be numerical
 * noise and can often be ignored or thresholded.
 */
std::pair<std::unique_ptr<Trajectory>, std::unique_ptr<Trajectory>> GDOP::finalize_optimal_bound_duals(
    const FixedVector<f64>& opt_z_lb,
    const FixedVector<f64>& opt_z_ub)
{
    auto optimal_bound_duals = std::make_pair<std::unique_ptr<Trajectory>, std::unique_ptr<Trajectory>>(
        std::make_unique<Trajectory>(),
        std::make_unique<Trajectory>()
    );

    // transform z_lb, z_ub from NLP bound duals -> costates
    FixedVector<f64> costates_z_lb(opt_z_lb);
    FixedVector<f64> costates_z_ub(opt_z_ub);

    transform_duals_costates_bounds(costates_z_lb, true);
    transform_duals_costates_bounds(costates_z_ub, true);

    auto& lower_traj = *(optimal_bound_duals.first);
    auto& upper_traj = *(optimal_bound_duals.second);

    std::vector<Trajectory*> bound_trajs = { &lower_traj, &upper_traj };
    std::vector<FixedVector<f64>*> bound_vals = { &costates_z_lb, &costates_z_ub };

    for (short bound_idx = 0; bound_idx < 2; bound_idx++) {
        Trajectory& traj         = *bound_trajs[bound_idx];
        FixedVector<f64>& z_dual = *bound_vals[bound_idx];

        traj.t.reserve(mesh->node_count + 1);
        traj.x.resize(off_x);
        traj.u.resize(off_u);
        traj.p.resize(off_p);
        traj.inducing_mesh = mesh->shared_from_this();
        traj.interpolation = InterpolationMethod::POLYNOMIAL;

        for (auto& v : traj.x) { v.reserve(mesh->node_count + 1); }
        for (auto& v : traj.u) { v.reserve(mesh->node_count + 1); }

        for (int x_index = 0; x_index < off_x; x_index++) {
            traj.x[x_index].push_back(z_dual[x_index]);
        }

        for (int u_index = 0; u_index < off_u; u_index++) {
            f64 u0 = fLGR::interpolate(mesh->nodes[0], false, &z_dual[2 * off_x + u_index], off_xu, mesh->t[0][0], mesh->grid[1], 0.0);
            traj.u[u_index].push_back(u0);
        }

        traj.t.push_back(0.0);

        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                for (int x_index = 0; x_index < off_x; x_index++) {
                    traj.x[x_index].push_back(z_dual[off_acc_xu[i][j] + x_index]);
                }

                for (int u_index = 0; u_index < off_u; u_index++) {
                    traj.u[u_index].push_back(z_dual[off_acc_xu[i][j] + off_x + u_index]);
                }

                traj.t.push_back(mesh->t[i][j]);
            }
        }

        for (int p_idx = 0; p_idx < off_p; p_idx++) {
            traj.p[p_idx] = z_dual[off_xu_total + p_idx];
        }
    }

    return optimal_bound_duals;
}

/**
 * @brief Finalizes the complete optimal solution, including primals, costates, and bound duals.
 *
 * This function orchestrates the finalization process by calling `finalize_optimal_primals()`,
 * `finalize_optimal_costates()`, and `finalize_optimal_bound_duals()` to gather all
 * relevant parts of the optimal solution. The results are then combined into a
 * `PrimalDualTrajectory` object and stored in the `optimal_solution` member.
 */
void GDOP::finalize_solution(
    f64 opt_obj,
    const FixedVector<f64>& opt_x,
    const FixedVector<f64>& opt_lambda,
    const FixedVector<f64>& opt_z_lb,
    const FixedVector<f64>& opt_z_ub)
{
    auto optimal_primals     = finalize_optimal_primals(opt_x);
    auto optimal_costates    = finalize_optimal_costates(opt_lambda);
    auto optimal_lu_costates = finalize_optimal_bound_duals(opt_z_lb, opt_z_ub);
    optimal_solution         = std::make_unique<PrimalDualTrajectory>(std::move(optimal_primals),
                                                                      std::move(optimal_costates),
                                                                      std::move(optimal_lu_costates.first),
                                                                      std::move(optimal_lu_costates.second));
}

} // namespace GDOP
