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

GDOP::GDOP(Problem& problem)
    : NLP::NLP(),
      mesh(problem.pc->mesh),
      spectral_mesh(problem.pc->free_time ? std::dynamic_pointer_cast<SpectralMesh>(mesh) : nullptr),
      problem(problem)
{
    if (problem.pc->free_time && spectral_mesh == nullptr) {
        Log::error("Problem has free initial or final time, but provided Mesh is not SpectralMesh.");
        std::abort();
    }
}

void GDOP::update(std::shared_ptr<Mesh> new_mesh) {
    // update mesh
    mesh = new_mesh;

    // update spectral mesh
    spectral_mesh = problem.pc->free_time ? std::dynamic_pointer_cast<SpectralMesh>(mesh) : nullptr;
    if (problem.pc->free_time && spectral_mesh == nullptr) {
        Log::error("Problem has free initial or final time, but provided Mesh is not SpectralMesh.");
        std::abort();
    }

    // set the new mesh and update the callback buffers with new sizes
    problem.update_mesh(new_mesh);
}

void GDOP::create_acc_offset_xu(int off_xu) {
    off_acc_xu = FixedField<int, 2>(mesh->intervals);
    int off = off_xu;
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
    create_acc_offset_xu(off_xu);                                    // variables  x_ij offset
    off_last_xu = off_acc_xu.back().back();                          // variables final grid point x_ij
    off_xu_total = off_last_xu + off_xu;                             // first parameter
    off_xup_total = off_xu_total + off_p;                            // first time variable
    number_vars = off_xup_total + (spectral_mesh ? 2 : 0);           // if free time -> += 2 for t0 and tf
    create_acc_offset_fg(problem.pc->fg_size);                       // constraint f_ij offset
    off_fg_total = mesh->node_count * problem.pc->fg_size;           // constraint r_0 offset
    off_fgr_total = off_fg_total + problem.pc->r_size;               // constraint a offset (initial control constraint)
    number_constraints = off_fgr_total + off_u;                      // f, r and a
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
    // standard bounds, but checking for xu0_fixed or xuf_fixed
    for (int x_index = 0; x_index < off_x; x_index++) {
        x_lb[x_index] = problem.pc->xu0_fixed[x_index] ? *problem.pc->xu0_fixed[x_index] : problem.pc->x_bounds[x_index].lb;
        x_ub[x_index] = problem.pc->xu0_fixed[x_index] ? *problem.pc->xu0_fixed[x_index] : problem.pc->x_bounds[x_index].ub;
    }

    for (int u_index = 0; u_index < off_u; u_index++) {
        int xu_idx = off_x + u_index;
        x_lb[xu_idx] = problem.pc->xu0_fixed[xu_idx] ? *problem.pc->xu0_fixed[xu_idx] : problem.pc->u_bounds[u_index].lb;
        x_ub[xu_idx] = problem.pc->xu0_fixed[xu_idx] ? *problem.pc->xu0_fixed[xu_idx] : problem.pc->u_bounds[u_index].ub;
    }

    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (i == mesh->intervals - 1 && j == mesh->nodes[i] - 1) {
                for (int x_index = 0; x_index < off_x; x_index++) {
                    x_lb[off_acc_xu[i][j] + x_index] = problem.pc->xuf_fixed[x_index] ? *problem.pc->xuf_fixed[x_index] : problem.pc->x_bounds[x_index].lb;
                    x_ub[off_acc_xu[i][j] + x_index] = problem.pc->xuf_fixed[x_index] ? *problem.pc->xuf_fixed[x_index] : problem.pc->x_bounds[x_index].ub;
                }

                for (int u_index = 0; u_index < off_u; u_index++) {
                    int xu_idx = off_x + u_index;
                    x_lb[off_acc_xu[i][j] + xu_idx] = problem.pc->xuf_fixed[xu_idx] ? *problem.pc->xuf_fixed[xu_idx] : problem.pc->u_bounds[u_index].lb;
                    x_ub[off_acc_xu[i][j] + xu_idx] = problem.pc->xuf_fixed[xu_idx] ? *problem.pc->xuf_fixed[xu_idx] : problem.pc->u_bounds[u_index].ub;
                }
            }
            else {
                for (int x_index = 0; x_index < off_x; x_index++) {
                    x_lb[off_acc_xu[i][j] + x_index] = problem.pc->x_bounds[x_index].lb;
                    x_ub[off_acc_xu[i][j] + x_index] = problem.pc->x_bounds[x_index].ub;
                }

                for (int u_index = 0; u_index < off_u; u_index++) {
                    x_lb[off_acc_xu[i][j] + off_x + u_index] = problem.pc->u_bounds[u_index].lb;
                    x_ub[off_acc_xu[i][j] + off_x + u_index] = problem.pc->u_bounds[u_index].ub;
                }
            }
        }
    }

    for (int p_index = 0; p_index < off_p; p_index++) {
        x_lb[off_xu_total + p_index] = problem.pc->p_bounds[p_index].lb;
        x_ub[off_xu_total + p_index] = problem.pc->p_bounds[p_index].ub;
    }

    if (spectral_mesh) {
        for (int t_index = 0; t_index < 2; t_index++) {
            x_lb[off_xup_total + t_index] = problem.pc->T_bounds[t_index].lb;
            x_ub[off_xup_total + t_index] = problem.pc->T_bounds[t_index].ub;
        }
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

    for (int a_index = 0; a_index < off_u; a_index++) {
        g_lb[off_fgr_total + a_index] = 0;
        g_ub[off_fgr_total + a_index] = 0;
    }
}

void GDOP::set_initial_guess(std::unique_ptr<PrimalDualTrajectory> initial_trajectory) {
    initial_guess = std::move(initial_trajectory);
}

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

        /**
         * debug only: first check tests if traj has inducing mesh, then the provided mesh should be equal to it
         *            second check tests if time grid are compatible
         * @note one of both should succeed, else we have a problem!
         */
        assert((!initial_guess_primal->inducing_mesh || initial_guess_primal->inducing_mesh == mesh)
                || check_time_compatibility(initial_guess_primal->t, {initial_guess_primal->x, initial_guess_primal->u}, *mesh));

        flatten_trajectory_to_layout(*initial_guess_primal, x_init, false);
    }
    else {
        Log::error("No primal initial guess supplied in GDOP::init_starting_point().");
    }

    if (initial_guess_costate) {
        // check compatibility
        if (initial_guess_costate->inducing_mesh.get() != mesh.get()) {
            initial_guess_costate = std::make_unique<CostateTrajectory>(initial_guess_costate->interpolate_onto_mesh(*mesh));
        }

        /**
         * debug only: first check tests if traj has inducing mesh, then the provided mesh should be equal to it
         *            second check tests if time grid are compatible
         * @note one of both should succeed, else we have a problem!
         */
        assert((!initial_guess_costate->inducing_mesh || initial_guess_costate->inducing_mesh == mesh)
                || check_time_compatibility(initial_guess_costate->t, {initial_guess_costate->costates_f, initial_guess_costate->costates_g}, *mesh));

        int index = 1; // ignore interpolated costates at t = 0
        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                for (int f_index = 0; f_index < problem.pc->f_size; f_index++) {
                    lambda_init[off_acc_fg[i][j] + f_index] = initial_guess_costate->costates_f[f_index][index];
                }
                for (int g_index = 0; g_index < problem.pc->g_size; g_index++) {
                    lambda_init[off_acc_fg[i][j] + problem.pc->f_size + g_index] = initial_guess_costate->costates_g[g_index][index];

                }
                index++;
            }
        }

        for (int r_index = 0; r_index < problem.pc->r_size; r_index++) {
            lambda_init[off_fg_total + r_index] = initial_guess_costate->costates_r[r_index];
        }

        for (int a_index = 0; a_index < off_u; a_index++) {
            lambda_init[off_fgr_total + a_index] = initial_guess_costate->costates_r[problem.pc->r_size + a_index];
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

        /**
         * debug only: first check tests if traj has inducing mesh, then the provided mesh should be equal to it
         *            second check tests if time grid are compatible
         * @note one of both should succeed, else we have a problem!
         */
        assert((!initial_guess_lower_costates->inducing_mesh || initial_guess_lower_costates->inducing_mesh == mesh)
                || check_time_compatibility(initial_guess_lower_costates->t, {initial_guess_lower_costates->x, initial_guess_lower_costates->u}, *mesh));

        assert((!initial_guess_upper_costates->inducing_mesh || initial_guess_upper_costates->inducing_mesh == mesh)
                || check_time_compatibility(initial_guess_upper_costates->t, {initial_guess_upper_costates->x, initial_guess_upper_costates->u}, *mesh));

        // get lower and upper costates (assume that if free time optimization -> guess contains duals for time vars at end of parameter vector!)
        flatten_trajectory_to_layout(*initial_guess_lower_costates, z_lb_init, true);
        flatten_trajectory_to_layout(*initial_guess_upper_costates, z_ub_init, true);

        transform_duals_costates_bounds(z_lb_init, false);
        transform_duals_costates_bounds(z_ub_init, false);
    }
}


void GDOP::init_jacobian_nonzeros(int& nnz_jac) {
    // stage 1: calculate nnz of blocks and number of collisions, where df_k / dx_k != 0. these are contained by default because of the D-Matrix
    int nnz_f = 0;
    int nnz_g = 0;
    int nnz_r = 0;
    int nnz_a = 0;
    int diagonal_collisions = 0;
    for (int f_index = 0; f_index < problem.pc->f_size; f_index++) {
        for (const auto& df_k_dx : problem.full->layout.f[f_index].jac.dx) {
            if (df_k_dx.col == f_index) {
                diagonal_collisions++;
            }
        }
        nnz_f += problem.full->layout.f[f_index].jac.nnz();

        if (spectral_mesh) {
            nnz_f += 2; // every dynamic constraint has non-zero derivative w.r.t. t0 and tf, since D * x - deltaT * f() = 0 depends on deltaT
        }
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
        nnz_r += problem.boundary->layout.r[r_index].jac.nnz_no_time();

        /**
         * @note 1 - Information about fixed time optimizations, when dM_dT or dr_dT have non-zero elements
         *
         * @attention: if the problem is optimized on a fixed time horizon, but the Boundary (Mayer M + Boundary Constraint r) is time-dependent,
         *             then we simply ignore these constraints, as the time variables (t0, tf) are not included / relevant for the NLP
         * @note: the callback buffers may still contain these derivatives since the user may supply these, therefore the buffers
         *        must be sufficiently large but are not used further
         *        thus, we guard every place that assumes / uses the time derivatives / time variables with if (spectral_mesh) == if (free time optimization)
         */
        if (spectral_mesh) {
            nnz_r += problem.boundary->layout.r[r_index].jac.nnz_time();
        }
    }

    nnz_a = off_u * (mesh->nodes[0] + 1);

    nnz_jac = off_acc_jac_fg.back() + nnz_r + nnz_a;

    // allocate memory for constant part
    const_der_jac = FixedVector<f64>(nnz_jac);
}

void GDOP::init_hessian_nonzeros(int& nnz_hes) {
    // takes O(nnz(A) + nnz(B) + ...+ nnz(H)) for creation of ** Maps and O(nnz(Hessian)) for creation of Hessian sparsity pattern

    // reset block sparsities
    hes_a_block = BlockSparsity::create_lower_triangular(problem.pc->xu_size, BlockType::Exact);
    hes_b_block = BlockSparsity::create_lower_triangular(problem.pc->x_size + problem.pc->u_size, BlockType::Offset);
    hes_c_block = BlockSparsity::create_rectangular(problem.pc->xu_size, problem.pc->xu_size, BlockType::Exact);
    hes_d_block = BlockSparsity::create_lower_triangular(problem.pc->xu_size, BlockType::Exact);
    hes_e_block = BlockSparsity::create_rectangular(problem.pc->p_size, problem.pc->xu_size, BlockType::Exact);
    hes_f_block = BlockSparsity::create_rectangular(problem.pc->p_size, problem.pc->xu_size, BlockType::RowOffset);
    hes_g_block = BlockSparsity::create_rectangular(problem.pc->p_size, problem.pc->xu_size, BlockType::Exact);
    hes_h_block = BlockSparsity::create_lower_triangular(problem.pc->p_size, BlockType::Exact);

    if (spectral_mesh) {
        hes_i_block = BlockSparsity::create_rectangular(2, problem.pc->xu_size, BlockType::Exact);
        hes_j_block = DenseRectangularBlockSparsity::create(2, problem.pc->xu_size);
        hes_k_block = BlockSparsity::create_rectangular(2, problem.pc->p_size, BlockType::Exact);
        hes_l_block = BlockSparsity::create_lower_triangular(2, BlockType::Exact);
    }

    // clear previous sparsity (may be reused) + setting boolean for lower triangular violation asserts
    hes_A_set.clear(true); hes_B_set.clear(true); hes_C_set.clear(false); hes_D_set.clear(true);
    hes_E_set.clear(false); hes_F_set.clear(false); hes_G_set.clear(false); hes_H_set.clear(true);

    if (spectral_mesh) {
        hes_I_set.clear(false);
        // J is dense
        hes_K_set.clear(false);
        hes_L_set.clear(true);
    }

    auto& boundary_hes = problem.boundary->layout.hes;
    auto& full_hes     = problem.full->layout.hes;
    auto& full_pp_hes  = problem.full->layout.pp_hes;

    // calculate IndexSet and nnz
    hes_A_set.insert_sparsity(boundary_hes.dx0_dx0,     0,     0);
    hes_A_set.insert_sparsity(boundary_hes.du0_dx0, off_x,     0);
    hes_A_set.insert_sparsity(boundary_hes.du0_du0, off_x, off_x);
    hes_C_set.insert_sparsity(boundary_hes.dxf_dx0,     0,     0);
    hes_C_set.insert_sparsity(boundary_hes.dxf_du0,     0, off_x);
    hes_D_set.insert_sparsity(boundary_hes.dxf_dxf,     0,     0);
    hes_C_set.insert_sparsity(boundary_hes.duf_dx0, off_x,     0);
    hes_C_set.insert_sparsity(boundary_hes.duf_du0, off_x, off_x);
    hes_D_set.insert_sparsity(boundary_hes.duf_dxf, off_x,     0);
    hes_D_set.insert_sparsity(boundary_hes.duf_duf, off_x, off_x);
    hes_E_set.insert_sparsity(boundary_hes.dp_dx0,      0,     0);
    hes_E_set.insert_sparsity(boundary_hes.dp_du0,      0, off_x);
    hes_G_set.insert_sparsity(boundary_hes.dp_dxf,      0,     0);
    hes_G_set.insert_sparsity(boundary_hes.dp_duf,      0, off_x);
    hes_H_set.insert_sparsity(boundary_hes.dp_dp,       0,     0);

    if (spectral_mesh) {
        hes_I_set.insert_sparsity(boundary_hes.dT_dx0,      0,     0);
        hes_I_set.insert_sparsity(boundary_hes.dT_du0,      0, off_x);
        hes_K_set.insert_sparsity(boundary_hes.dT_dp,       0,     0);
        hes_L_set.insert_sparsity(boundary_hes.dT_dT,       0,     0);
    }

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

    /**
     * @note K block can be defined as struct(K) = [0, 1] X [struct(L_p) + \sum_{f \in Dynamics} struct(f_p)]
     *       these stem from the dynamic constraints (same for Lagrange term): F = D * x - deltat * f(x, u, p)
     *        -> writing it with functions we get L_1(x) - L_2(tf, t0) * f(x, u, p), where L denotes a linear function.
     *        Because of the product rule, we get nabla² F_{(t0, tf), p} = dL_2(tf, t0) / d{t0, tf} * df(x, u, p) / dp
     *        -> for the Lagrangian Hessian we sum over all these f and L, so we get the above formula for the sparsity struct
     * @note we assume that dL_2(tf, t0) / d{t0, tf} * df(x, u, p) / d{x, u} is fully dense anyway,
     *       so derivatives w.r.t. x and u are not inserted into a block sparsity and are given by a DenseRectangularBlockSparsity (Block J)
     */
    if (spectral_mesh) {
        std::vector<int> rows = { 0, 1 };

        if (problem.pc->has_lagrange) {
            hes_K_set.insert_sparsity(rows, problem.full->layout.L->jac.dp, 0, 0);
        }

        for (auto const& f : problem.full->layout.f) {
            hes_K_set.insert_sparsity(rows, f.jac.dp, 0, 0);
        }
    }

    // calculate nnz from block sparsity
    nnz_hes = (hes_B_set.size() + hes_F_set.size()) * (mesh->node_count - 1)
              + hes_A_set.size() + hes_C_set.size() + hes_D_set.size() + hes_E_set.size() + hes_G_set.size() + hes_H_set.size();

    if (spectral_mesh) {
        nnz_hes += hes_I_set.size() + (2 * off_xu) * mesh->node_count + hes_K_set.size() + hes_L_set.size();
    }
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

                // df / dT
                if (spectral_mesh) {
                    // df / dt0
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_xup_total;
                    nnz_index++;

                    // df / dtf
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_xup_total + 1;
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

        // dr / du0
        for (auto& dr_du0 : problem.boundary->layout.r[r_index].jac.du0) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = off_x + dr_du0.col;
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

        // dr / dT
        if (spectral_mesh) {
            // only use the derivative w.r.t. time if free time optimization (see "@note 1")
            for (auto& dr_dT: problem.boundary->layout.r[r_index].jac.dT) {
                i_row_jac[nnz_index] = eqn_index;
                j_col_jac[nnz_index] = off_xup_total + dr_dT.col;
                nnz_index++;
            }
        }
    }

    for (int a_index = 0; a_index < off_u; a_index++) {
        int eqn_index = off_fgr_total + a_index;

        for (int j = 0; j < mesh->nodes[0] + 1; j++) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = off_x + a_index + j * off_xu;
            // This interpolation identity somehow holds: u(t0) = sum_(j >= 1) c[j] * D[0, j] * u(t_0j)
            const_der_jac[nnz_index] = (j == 0) ? 1.0 : -1.0 * fLGR::get_c0(mesh->nodes[0], j) * fLGR::get_D(mesh->nodes[0], 0, j);
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
        i_row_hes[hes_nnz_counter] = row; // xu_0
        j_col_hes[hes_nnz_counter] = col; // xu_0
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
            j_col_hes[hes_nnz_counter] = C_flat[c_index].second;              // xu_0
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
            j_col_hes[hes_nnz_counter] = E_flat[e_index].second;               // xu_0
            hes_e_block.insert(E_flat[e_index].first, E_flat[e_index].second, hes_nnz_counter++);
            e_index++;
        }

        int row_f_nnz = 0;
        hes_f_block.row_offset_prev[p_index] = hes_nnz_counter; // E_{p_index, :} offset
        while (f_index < F_flat.int_size() && F_flat[f_index].first == p_index) {
            hes_f_block.insert(F_flat[f_index].first, F_flat[f_index].second, row_f_nnz++);
            f_index++;
        }
        /* F_{p_index, :} size -> offset for next F blocks in same row of full sparsity */
        hes_f_block.row_size[p_index] = row_f_nnz;
        hes_nnz_counter += (mesh->node_count - 1) * row_f_nnz;

        while (g_index < G_flat.int_size() && G_flat[g_index].first == p_index) {
            i_row_hes[hes_nnz_counter] = off_xu_total + G_flat[g_index].first; // p
            j_col_hes[hes_nnz_counter] = off_last_xu + G_flat[g_index].second; // xu_{nm}
            hes_g_block.insert(G_flat[g_index].first, G_flat[g_index].second, hes_nnz_counter++);
            g_index++;
        }

        while (h_index < H_flat.int_size() && H_flat[h_index].first == p_index) {
            i_row_hes[hes_nnz_counter] = off_xu_total + H_flat[h_index].first;  // p
            j_col_hes[hes_nnz_counter] = off_xu_total + H_flat[h_index].second; // p
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

    if (spectral_mesh) {
        // I, J, K, L: partially exact (all except J) with dense row dependence
        int i_index = 0;
        int k_index = 0;
        int l_index = 0;
        FixedVector<std::pair<int, int>> I_flat(hes_I_set.set.begin(), hes_I_set.set.end());
        FixedVector<std::pair<int, int>> K_flat(hes_K_set.set.begin(), hes_K_set.set.end());
        FixedVector<std::pair<int, int>> L_flat(hes_L_set.set.begin(), hes_L_set.set.end());
        for (int t_index = 0; t_index < 2; t_index++) {
            while (i_index < I_flat.int_size() && I_flat[i_index].first == t_index) {
                i_row_hes[hes_nnz_counter] = off_xup_total + I_flat[i_index].first; // T
                j_col_hes[hes_nnz_counter] = I_flat[i_index].second;                // xu_0
                hes_i_block.insert(I_flat[i_index].first, I_flat[i_index].second, hes_nnz_counter++);
                i_index++;
            }

            // init J hessian pattern O(node_count * (#x + #u))
            hes_j_block.row_offset_prev[t_index] = hes_nnz_counter; // I_{t_index, :} offset
            const int len = off_xu * mesh->node_count;
            std::fill(&i_row_hes[hes_nnz_counter], /* care: too far ptr */ &i_row_hes[hes_nnz_counter] + len, off_xup_total + t_index);
            std::iota(&j_col_hes[hes_nnz_counter], /* care: too far ptr */ &j_col_hes[hes_nnz_counter] + len, get_off_first_xu());
            hes_nnz_counter += len;

            while (k_index < K_flat.int_size() && K_flat[k_index].first == t_index) {
                i_row_hes[hes_nnz_counter] = off_xup_total + K_flat[k_index].first; // T
                j_col_hes[hes_nnz_counter] = off_xu_total + K_flat[k_index].second; // p
                hes_k_block.insert(K_flat[k_index].first, K_flat[k_index].second, hes_nnz_counter++);
                k_index++;
            }

            while (l_index < L_flat.int_size() && L_flat[l_index].first == t_index) {
                i_row_hes[hes_nnz_counter] = off_xup_total + L_flat[l_index].first;  // T
                j_col_hes[hes_nnz_counter] = off_xup_total + L_flat[l_index].second; // T
                hes_l_block.insert(L_flat[l_index].first, L_flat[l_index].second, hes_nnz_counter++);
                l_index++;
            }
        }

        assert(hes_nnz_counter == get_nnz_hes());
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
 *   | Function      | `new_x` | `new_lambda` | Triggers Callback       | Needs callback_jacobian | Needs callback_evaluation |
 *   | ------------- | ------- | ------------ | ----------------------- | ------------------------| --------------------------|
 *   | `eval_f`      | +       | -            | `callback_evaluation()` | -                       | +                         |
 *   | `eval_g`      | +       | -            | `callback_evaluation()` | -                       | +                         |
 *   | `eval_grad_f` | +       | -            | `callback_jacobian()`   | +                       | +                         |
 *   | `eval_jac_g`  | +       | -            | `callback_jacobian()`   | +                       | +                         |
 *   | `eval_hes`    | +       | +            | `callback_hessian()`    | +                       | +  (transitive)           |
 *
 *  Because of free time optimizations and possible user numerical evaluations of callback functions, the solver callbacks
 *  eval_f, eval_g, eval_grad_f, eval_jac_g and eval_hes do not only require the evaluation itself, but also all prior callbacks to be available.
 *                - eval_f and eval_g only require callback_evaluation()
 *                - eval_grad_f and eval_jac_g require callback_evaluation() and callback_jacobian()
 *                - eval_hes requires callback_jacobian() (and thus transitively callback_evaluation()) and callback_hessian() (so all)
 *
 */

// ========= virtuals in NLP =========

// === overload ===
// evaluate objective function
void GDOP::eval_f(
    bool new_x,
    const FixedVector<f64>& curr_x,
    f64& curr_obj)
{
    check_new_x(new_x, curr_x);
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
    check_new_x(new_x, curr_x);
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
    check_new_x(new_x, curr_x);
    if (!evaluation_state.eval_f) {
        callback_evaluation(curr_x);
    }

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
    check_new_x(new_x, curr_x);
    if (!evaluation_state.eval_g) {
        callback_evaluation(curr_x);
    }

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
    check_new_x(new_x, curr_x);
    check_new_lambda(new_lambda);

    // ensure Jacobian is available (which itself needs the evaluation) (required for numerical Hessian and esp. free time optimization)
    if (!evaluation_state.hes) {
        if (!(evaluation_state.jac_g && evaluation_state.grad_f)) {
            if (!(evaluation_state.eval_f && evaluation_state.eval_g)) {
                callback_evaluation(curr_x);
            }
            callback_jacobian(curr_x);
        }
        callback_hessian(curr_x, curr_lambda, curr_obj_factor);
    }
    eval_hes_internal(curr_hes, curr_lambda, curr_obj_factor);
}

// ========= callbacks and internal evaluation =========

// check if a new x was received; if so, reset evaluation state and update time variables for spectral meshes
void GDOP::check_new_x(bool new_x, const FixedVector<f64>& curr_x) {
    evaluation_state.check_reset_x(new_x);
    if (problem.pc->free_time) {
        spectral_mesh->update_physical_from_spectral(*get_x_t0(curr_x),*get_x_tf(curr_x));
    }
}

// similar check for lambda (dual variables).
void GDOP::check_new_lambda(bool new_lambda) {
    evaluation_state.check_reset_lambda(new_lambda);
}

void GDOP::callback_evaluation(const FixedVector<f64>& curr_x) {
    problem.full->callback_eval(get_x_xu(curr_x), get_x_p(curr_x));
    problem.boundary->callback_eval(get_x_xu0(curr_x), get_x_xuf(curr_x), get_x_p(curr_x), mesh->t0, mesh->tf);
    evaluation_state.eval_f = true;
    evaluation_state.eval_g = true;
}

void GDOP::callback_jacobian(const FixedVector<f64>& curr_x) {
    problem.full->callback_jac(get_x_xu(curr_x), get_x_p(curr_x));
    problem.boundary->callback_jac(get_x_xu0(curr_x), get_x_xuf(curr_x), get_x_p(curr_x), mesh->t0, mesh->tf);
    evaluation_state.grad_f = true;
    evaluation_state.jac_g = true;
}

// perform update of dual variables, such that callback can use the exact multiplier
void GDOP::update_curr_lambda_obj_factors(const FixedVector<f64>& curr_lambda, f64 curr_obj_factor) {
    transformed_lambda = curr_lambda; // copy curr_lambda

    for (int i = 0; i < mesh->intervals; i++) {
        f64 delta_t = mesh->delta_t[i];
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (problem.pc->has_lagrange) {
                lagrange_obj_factors[i][j] = curr_obj_factor * fLGR::get_b(mesh->nodes[i], j) * delta_t;
            }
            for (int f = 0; f < problem.pc->f_size; f++) {
                transformed_lambda[off_acc_fg[i][j] + f] *= -delta_t;
            }
        }
    }
}

void GDOP::callback_hessian(const FixedVector<f64> x, const FixedVector<f64>& curr_lambda, f64 curr_obj_factor) {
    update_curr_lambda_obj_factors(curr_lambda, curr_obj_factor);

    problem.full->callback_hes(get_x_xu(x), get_x_p(x), lagrange_obj_factors, get_lmbd_fg(transformed_lambda));
    problem.boundary->callback_hes(get_x_xu0(x), get_x_xuf(x), get_x_p(x), mesh->t0, mesh->tf, curr_obj_factor, get_lmbd_r(transformed_lambda));
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

                // guard see "@note 1"
                if (spectral_mesh) {
                    curr_grad[off_xup_total]     -= spectral_mesh->delta_tau(i) * fLGR::get_b(mesh->nodes[i], j) * problem.lfg_eval_L(i, j);
                    curr_grad[off_xup_total + 1] += spectral_mesh->delta_tau(i) * fLGR::get_b(mesh->nodes[i], j) * problem.lfg_eval_L(i, j);
                }
            }
        }
    }
    if (problem.pc->has_mayer) {
        for (auto& dM_dx0 : problem.boundary->layout.M->jac.dx0) {
            curr_grad[dM_dx0.col] += problem.mr_jac(dM_dx0.buf_index);
        }
        for (auto& dM_du0 : problem.boundary->layout.M->jac.du0) {
            curr_grad[off_x + dM_du0.col] += problem.mr_jac(dM_du0.buf_index);
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

        // guard see "@note 1"
        if (spectral_mesh) {
            for (auto& dM_dT : problem.boundary->layout.M->jac.dT) {
                curr_grad[off_xup_total + dM_dT.col] += problem.mr_jac(dM_dT.buf_index);
            }
        }
    }
};

void GDOP::eval_g_internal(const FixedVector<f64>& curr_x, FixedVector<f64>& curr_g) {
    curr_g.fill_zero();
    for (int i = 0; i < mesh->intervals; i++) {
        fLGR::diff_matrix_multiply_block_strided(mesh->nodes[i], off_x, off_xu, problem.pc->fg_size,
                                                 &curr_x[i == 0 ? 0 : off_acc_xu[i - 1][mesh->nodes[i - 1] - 1]],  // x_{i-1, m_{i-1}} base point states
                                                 &curr_x[off_acc_xu[i][0]],                                        // collocation point states
                                                 &curr_g[off_acc_fg[i][0]]);                                       // constraint start index
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

    const f64 *u0 = &get_x_xu0(curr_x)[off_x];
    for (int a_index = 0; a_index < off_u; a_index++) {
        for (int j = 0; j < mesh->nodes[0] + 1; j++) {
            if (j == 0) {
                curr_g[off_fgr_total + a_index] += u0[a_index];
            }
            else {
                curr_g[off_fgr_total + a_index] -= fLGR::get_c0(mesh->nodes[0], j) * fLGR::get_D(mesh->nodes[0], 0, j) * u0[a_index + j * off_xu];
            }
        }
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

                // df / dT
                if (spectral_mesh) {
                    // df / dt0
                    curr_jac[nnz_index++] = spectral_mesh->delta_tau(i) * problem.lfg_eval_f(f_index, i, j);

                    // df / dtf
                    curr_jac[nnz_index++] = -spectral_mesh->delta_tau(i) * problem.lfg_eval_f(f_index, i, j);
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

        // dr / du0
        for (auto& dr_du0 : problem.boundary->layout.r[r_index].jac.du0) {
            curr_jac[nnz_index++] = problem.mr_jac(dr_du0.buf_index);
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

        // dr / dT
        if (spectral_mesh) {
            // only use the derivative w.r.t. time if free time optimization (see "@note 1")
            for (auto& dr_dT: problem.boundary->layout.r[r_index].jac.dT) {
                curr_jac[nnz_index++] = problem.mr_jac(dr_dT.buf_index);
            }
        }
    }

    // offset for a (artificial initial cotrol constraints)
    nnz_index += off_u * (mesh->nodes[0] + 1);

    assert(nnz_index == get_nnz_jac());
};

void GDOP::eval_hes_internal(FixedVector<f64>& curr_hes, const FixedVector<f64>& curr_lambda, f64 curr_obj_factor) {
    curr_hes.fill_zero();

    for (int i = 0; i < mesh->intervals; i++) {
        // insert buffer for parallel execution here / pass the buffer in updateHessian(f64*) calls
        for (int j = 0; j < mesh->nodes[i]; j++) {
            // make sure to be in the right ptr_map region, B and F are only valid for i,j != n,m
            //                                              D and G are only valid for i,j == n,m
            if (!(i == mesh->intervals - 1 && j == mesh->nodes[mesh->intervals - 1] - 1)) {
                accumulate_hessian_lfg(problem.full->layout.hes, i, j, hes_b_block, hes_f_block, curr_hes);

            } else {
                accumulate_hessian_lfg(problem.full->layout.hes, i, j, hes_d_block, hes_g_block, curr_hes);
            }

            if (spectral_mesh) {
                // d{Lf} / { d{t0, tf} d{x_ij, u_ij, p}}
                accumulate_hessian_from_lagrangian_gradient_lf(i, j, curr_hes, curr_lambda, curr_obj_factor);
            }
        }
    }
    accumulate_hessian_parameter_lfg(problem.full->layout.pp_hes, curr_hes);
    accumulate_hessian_mr(problem.boundary->layout.hes, curr_hes);
}

void GDOP::accumulate_hessian_lfg(const HessianLFG& hes, int interval_i, int node_j,
                                  const BlockSparsity& ptr_map_xu_xu, const BlockSparsity& ptr_map_p_xu,
                                  FixedVector<f64>& curr_hes) {
    const int block_count = mesh->acc_nodes[interval_i][node_j];
    for (const auto& dx_dx : hes.dx_dx) {
        curr_hes[ptr_map_xu_xu.access(dx_dx.row, dx_dx.col, block_count)] += problem.lfg_hes(dx_dx.buf_index, interval_i, node_j);
    }
    for (const auto& du_dx : hes.du_dx) {
        curr_hes[ptr_map_xu_xu.access(off_x + du_dx.row, du_dx.col, block_count)] += problem.lfg_hes(du_dx.buf_index, interval_i, node_j);
    }
    for (const auto& du_du : hes.du_du) {
        curr_hes[ptr_map_xu_xu.access(off_x + du_du.row, off_x + du_du.col, block_count)] += problem.lfg_hes(du_du.buf_index, interval_i, node_j);
    }
    for (const auto& dp_dx : hes.dp_dx) {
        curr_hes[ptr_map_p_xu.access(dp_dx.row, dp_dx.col, block_count)] += problem.lfg_hes(dp_dx.buf_index, interval_i, node_j);
    }
    for (const auto& dp_du : hes.dp_du) {
        curr_hes[ptr_map_p_xu.access(dp_du.row, off_x + dp_du.col, block_count)] += problem.lfg_hes(dp_du.buf_index, interval_i, node_j);
    }
}

void GDOP::accumulate_hessian_parameter_lfg(const ParameterHessian& pp_hes, FixedVector<f64>& curr_hes) {
    for (const auto& dp_dp : pp_hes.dp_dp) {
        curr_hes[hes_h_block.access(dp_dp.row, dp_dp.col)] += problem.lfg_pp_hes(dp_dp.buf_index);
    }
}

void GDOP::accumulate_hessian_mr(const HessianMR& hes, FixedVector<f64>& curr_hes) {
    for (const auto& dx0_dx0 : hes.dx0_dx0) {
        curr_hes[hes_a_block.access(dx0_dx0.row, dx0_dx0.col)] += problem.mr_hes(dx0_dx0.buf_index);
    }
    for (const auto& du0_dx0 : hes.du0_dx0) {
        curr_hes[hes_a_block.access(off_x + du0_dx0.row, du0_dx0.col)] += problem.mr_hes(du0_dx0.buf_index);
    }
    for (const auto& du0_du0 : hes.du0_du0) {
        curr_hes[hes_a_block.access(off_x + du0_du0.row, off_x + du0_du0.col)] += problem.mr_hes(du0_du0.buf_index);
    }
    for (const auto& dxf_dx0 : hes.dxf_dx0) {
        curr_hes[hes_c_block.access(dxf_dx0.row, dxf_dx0.col)] += problem.mr_hes(dxf_dx0.buf_index);
    }
    for (const auto& dxf_du0 : hes.dxf_du0) {
        curr_hes[hes_c_block.access(dxf_du0.row, off_x + dxf_du0.col)] += problem.mr_hes(dxf_du0.buf_index);
    }
    for (const auto& dxf_dxf : hes.dxf_dxf) {
        curr_hes[hes_d_block.access(dxf_dxf.row, dxf_dxf.col)] += problem.mr_hes(dxf_dxf.buf_index);
    }
    for (const auto& duf_dx0 : hes.duf_dx0) {
        curr_hes[hes_c_block.access(off_x + duf_dx0.row, duf_dx0.col)] += problem.mr_hes(duf_dx0.buf_index);
    }
    for (const auto& duf_du0 : hes.duf_du0) {
        curr_hes[hes_c_block.access(off_x + duf_du0.row, off_x + duf_du0.col)] += problem.mr_hes(duf_du0.buf_index);
    }
    for (const auto& duf_dxf : hes.duf_dxf) {
        curr_hes[hes_d_block.access(off_x + duf_dxf.row, duf_dxf.col)] += problem.mr_hes(duf_dxf.buf_index);
    }
    for (const auto& duf_duf : hes.duf_duf) {
        curr_hes[hes_d_block.access(off_x + duf_duf.row, off_x + duf_duf.col)] += problem.mr_hes(duf_duf.buf_index);
    }
    for (const auto& dp_dx0 : hes.dp_dx0) {
        curr_hes[hes_e_block.access(dp_dx0.row, dp_dx0.col)] += problem.mr_hes(dp_dx0.buf_index);
    }
    for (const auto& dp_du0 : hes.dp_du0) {
        curr_hes[hes_e_block.access(dp_du0.row, off_x + dp_du0.col)] += problem.mr_hes(dp_du0.buf_index);
    }
    for (const auto& dp_dxf : hes.dp_dxf) {
        curr_hes[hes_g_block.access(dp_dxf.row, dp_dxf.col)] += problem.mr_hes(dp_dxf.buf_index);
    }
    for (const auto& dp_duf : hes.dp_duf) {
        curr_hes[hes_g_block.access(dp_duf.row, off_x + dp_duf.col)] += problem.mr_hes(dp_duf.buf_index);
    }
    for (const auto& dp_dp : hes.dp_dp) {
        curr_hes[hes_h_block.access(dp_dp.row, dp_dp.col)] += problem.mr_hes(dp_dp.buf_index);
    }

    if (spectral_mesh) {
        for (const auto& dT_dx0 : hes.dT_dx0) {
            curr_hes[hes_i_block.access(dT_dx0.row, dT_dx0.col)] += problem.mr_hes(dT_dx0.buf_index);
        }
        for (const auto& dT_du0 : hes.dT_du0) {
            curr_hes[hes_i_block.access(dT_du0.row, off_x + dT_du0.col)] += problem.mr_hes(dT_du0.buf_index);
        }
        for (auto& dT_dxf : hes.dT_dxf) {
            curr_hes[hes_j_block.access(dT_dxf.row, dT_dxf.col, mesh->node_count - 1)] += problem.mr_hes(dT_dxf.buf_index);
        }
        for (const auto& dT_duf : hes.dT_duf) {
            curr_hes[hes_j_block.access(dT_duf.row, off_x + dT_duf.col, mesh->node_count - 1)] += problem.mr_hes(dT_duf.buf_index);
        }
        for (const auto& dT_dp : hes.dT_dp) {
            curr_hes[hes_k_block.access(dT_dp.row, dT_dp.col)] += problem.mr_hes(dT_dp.buf_index);
        }
        for (const auto& dT_dT : hes.dT_dT) {
            curr_hes[hes_l_block.access(dT_dT.row, dT_dT.col)] += problem.mr_hes(dT_dT.buf_index);
        }
    }
}

/**
 * @brief Adds the cross terms (d{t0, tf} d{x, u, p}) of the dynamic constraints + Lagrange term to the Lagrangian Hessian
 *
 * @note The GDOP contains L and f that have the structure: some_constant * deltaT * {L, f}(x, u, p),
 *       we get the cross term Hessian of these constraints as sum_{phi_i = all functions f and L} dual_i * constant_i * delta_T' * d{phi_i} / d_{x, u, p}
 *       So w.r.t. tf and t0 we get: sum_{phi_i = all functions f and L} +- dual_i * tau_i * constant_i * d{phi_i} / d_{x, u, p}
 *       Clearly, this is a scaled gradient of the (we call) `partial Lagrangian`: obj_factor * L(x, u, p) + lambda^T * f(x, u, p)!
 * @note g is not multiplied with deltaT, so it doesnt appear in this partial Lagrangian
 */
void GDOP::accumulate_hessian_from_lagrangian_gradient_lf(int interval_i,
                                                          int node_j,
                                                          FixedVector<f64>& curr_hes,
                                                          const FixedVector<f64>& curr_lambda,
                                                          f64 curr_obj_factor)
{
    const int block_count = mesh->acc_nodes[interval_i][node_j];
    const auto& layout = problem.full->layout;

    if (problem.pc->has_lagrange) {
        f64 factor_tf = curr_obj_factor * spectral_mesh->delta_tau(interval_i) * fLGR::get_b(mesh->nodes[interval_i], node_j); // must be without minus sign!!
        accumulate_hessian_from_lagrangian_gradient_lf_jac(block_count, interval_i, node_j, factor_tf, layout.L->jac, curr_hes);
    }

    for (int f_index = 0; f_index < problem.pc->f_size; f_index++) {
        f64 factor_tf = -curr_lambda[off_acc_fg[interval_i][node_j] + f_index] * spectral_mesh->delta_tau(interval_i); // must be with minus sign!!
        accumulate_hessian_from_lagrangian_gradient_lf_jac(block_count, interval_i, node_j, factor_tf, layout.f[f_index].jac, curr_hes);
    }
}

void GDOP::accumulate_hessian_from_lagrangian_gradient_lf_jac(int block_count,
                                                              int interval_i,
                                                              int node_j,
                                                              f64 factor_tf, // factor_t0 = -factor_tf
                                                              const JacobianLFG& jac,
                                                              FixedVector<f64>& curr_hes)
{
    for (auto& ddx : jac.dx) {
        f64 res_tf = factor_tf * problem.lfg_jac(ddx.buf_index, interval_i, node_j);
        curr_hes[hes_j_block.access(/* t0 */ 0, ddx.col, block_count)] -= res_tf;
        curr_hes[hes_j_block.access(/* tf */ 1, ddx.col, block_count)] += res_tf;
    }

    for (auto& ddu : jac.du) {
        f64 res_tf = factor_tf * problem.lfg_jac(ddu.buf_index, interval_i, node_j);
        curr_hes[hes_j_block.access(/* t0 */ 0, ddu.col + off_x, block_count)] -= res_tf;
        curr_hes[hes_j_block.access(/* tf */ 1, ddu.col + off_x, block_count)] += res_tf;
    }

    for (auto& ddp : jac.dp) {
        f64 res_tf = factor_tf * problem.lfg_jac(ddp.buf_index, interval_i, node_j);
        curr_hes[hes_k_block.access(/* t0 */ 0, ddp.col, block_count)] -= res_tf;
        curr_hes[hes_k_block.access(/* tf */ 1, ddp.col, block_count)] += res_tf;
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

void GDOP::flatten_trajectory_to_layout(const Trajectory& trajectory, FixedVector<f64>& flat_buffer, bool from_costates) {
    for (int x_index = 0; x_index < off_x; x_index++) {
        flat_buffer[x_index] = trajectory.x[x_index][0];
    }

    for (int u_index = 0; u_index < off_u; u_index++) {
        flat_buffer[off_x + u_index] = trajectory.u[u_index][0];
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

    if (spectral_mesh) {
        if (!from_costates) {
            flat_buffer[off_xup_total] = trajectory.t.front();
            flat_buffer[off_xup_total + 1] = trajectory.t.back();
        }
        else {
            /** @note we include the z-duals of the time variables in the parameter vector at the end
             * so duals are still nicely visible in the CSV export
             * and we can use the real time variables for trajectory.t */
            assert(int_size(trajectory.p) == off_p + 2);
            flat_buffer[off_xup_total] = trajectory.p[off_p];
            flat_buffer[off_xup_total + 1] = trajectory.p[off_p + 1];
        }
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
        optimal_primals->u[u_index].push_back(opt_x[off_x + u_index]);
    }

    optimal_primals->t.push_back(mesh->t0);

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

    if (spectral_mesh) {
        // also write t0 and tf to parameter array
        optimal_primals->p.push_back(*get_x_t0(opt_x));
        optimal_primals->p.push_back(*get_x_tf(opt_x));
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
        f64 lambda_f_0 = fLGR::interpolate(mesh->nodes[0], false, &costates[f_index], inp_stride, mesh->t[0][0], mesh->grid[1], mesh->t0);
        optimal_costates->costates_f[f_index].push_back(lambda_f_0);
    }

    for (int g_index = 0; g_index < g_size; g_index++) {
        f64 lambda_g_0 = fLGR::interpolate(mesh->nodes[0], false, &costates[f_size + g_index], inp_stride, mesh->t[0][0], mesh->grid[1], mesh->t0);
        optimal_costates->costates_g[g_index].push_back(lambda_g_0);
    }

    optimal_costates->t.push_back(mesh->t0);

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

    optimal_costates->costates_r.resize(problem.pc->r_size + off_u);

    for (int r_index = 0; r_index < problem.pc->r_size; r_index++) {
        optimal_costates->costates_r[r_index] = costates[off_fg_total + r_index];
    }

    // TODO ?
    for (int a_index = 0; a_index < off_u; a_index++) {
        optimal_costates->costates_r[problem.pc->r_size + a_index] = costates[off_fgr_total + a_index];
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
        traj.p.resize(off_p + (spectral_mesh ? 2 : 0));
        traj.inducing_mesh = mesh->shared_from_this();
        traj.interpolation = InterpolationMethod::POLYNOMIAL;

        for (auto& v : traj.x) { v.reserve(mesh->node_count + 1); }
        for (auto& v : traj.u) { v.reserve(mesh->node_count + 1); }

        // x0, u0
        for (int x_index = 0; x_index < off_x; x_index++) {
            traj.x[x_index].push_back(z_dual[x_index]);
        }

        for (int u_index = 0; u_index < off_u; u_index++) {
            traj.u[u_index].push_back(z_dual[off_x + u_index]);
        }

        traj.t.push_back(mesh->t0);

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

        // include duals of time variables in the parameter vector for latter reinit from costates
        if (spectral_mesh) {
            traj.p[off_p] = z_dual[off_xup_total];
            traj.p[off_p + 1] = z_dual[off_xup_total + 1];
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
