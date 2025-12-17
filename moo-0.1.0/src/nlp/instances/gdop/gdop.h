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

#ifndef MOO_GDOP_H
#define MOO_GDOP_H

#include <cassert>

#include <base/block_sparsity.h>
#include <base/fLGR.h>
#include <base/trajectory.h>
#include <base/fixed_vector.h>
#include <base/linalg.h>
#include <base/nlp_structs.h>
#include <base/mesh.h>
#include <base/util.h>
#include <nlp/nlp.h>

#include "problem.h"
#include "strategies.h"

namespace GDOP {

using NLP::Scaling;

class MOO_EXPORT GDOP : public NLP::NLP {
public:
    GDOP(Problem& problem);

    // === API + external calls ===

    void update(std::shared_ptr<Mesh> new_mesh);

    void set_initial_guess(std::unique_ptr<PrimalDualTrajectory> initial_trajectory);

    void set_scaling_factory(std::shared_ptr<ScalingFactory> factory);

    // === Constant Reference Getters for Private Members ===

    // main objects
    inline const Mesh&               get_mesh()           const { return *mesh; }
    inline const Problem&            get_problem()        const { return problem; }

    // offsets and sizes
    inline int                       get_off_x()          const { return off_x; }
    inline int                       get_off_u()          const { return off_u; }
    inline int                       get_off_p()          const { return off_p; }
    inline int                       get_off_xu()         const { return off_xu; }
    inline int                       get_off_first_xu()   const { return off_xu; }
    inline int                       get_off_last_xu()    const { return off_last_xu; }
    inline int                       get_off_xu_total()   const { return off_xu_total; }
    inline int                       get_off_fg_total()   const { return off_fg_total; }
    inline int                       get_off_fgr_total()  const { return off_fgr_total; }
    inline const FixedField<int, 2>& get_off_acc_xu()     const { return off_acc_xu; }
    inline const FixedField<int, 2>& get_off_acc_fg()     const { return off_acc_fg; }
    inline const FixedVector<int>&   get_off_acc_jac_fg() const { return off_acc_jac_fg; }

    // optimal solution
    inline const PrimalDualTrajectory* get_optimal_solution() const { return optimal_solution.get(); }

    // === Virtuals in NLP ===

    void get_sizes(
      int& number_vars,
      int& number_constraints) override;

    std::shared_ptr<Scaling> get_scaling() override;

    void get_bounds(
        FixedVector<f64>& x_lb,
        FixedVector<f64>& x_ub,
        FixedVector<f64>& g_lb,
        FixedVector<f64>& g_ub) override;

    void get_initial_guess(
        bool init_x,
        FixedVector<f64>& x_init,
        bool init_lambda,
        FixedVector<f64>& lambda_init,
        bool init_z,
        FixedVector<f64>& z_lb_init,
        FixedVector<f64>& z_ub_init) override;

    void get_nnz(
        int& nnz_jac,
        int& nnz_hes) override;

    void get_jac_sparsity(
        FixedVector<int>& i_row_jac,
        FixedVector<int>& j_col_jac) override;

    void get_hes_sparsity(
        FixedVector<int>& i_row_hes,
        FixedVector<int>& j_col_hes) override;

    void eval_f(
        bool new_x,
        const FixedVector<f64>& curr_x,
        f64& curr_obj) override;

    void eval_g(
        bool new_x,
        const FixedVector<f64>& curr_x,
        FixedVector<f64>& curr_g) override;

    void eval_grad_f(
        bool new_x,
        const FixedVector<f64>& curr_x,
        FixedVector<f64>& curr_grad_f) override;

    void eval_jac_g(
        bool new_x,
        const FixedVector<f64>& curr_x,
        const FixedVector<int>& i_row_jac,
        const FixedVector<int>& j_col_jac,
        FixedVector<f64>& curr_jac) override;

    void eval_hes(
        bool new_x,
        const FixedVector<f64>& curr_x,
        bool new_lambda,
        const FixedVector<f64>& curr_lambda,
        f64 curr_obj_factor,
        const FixedVector<int>& i_row_hes,
        const FixedVector<int>& j_col_hes,
        FixedVector<f64>& curr_hes) override;

    void finalize_solution(
        f64 opt_obj,
        const FixedVector<f64>& opt_x,
        const FixedVector<f64>& opt_lambda,
        const FixedVector<f64>& opt_z_lb,
        const FixedVector<f64>& opt_z_ub) override;

private:

    // === private structures ===

    std::shared_ptr<Mesh> mesh;                   // stndard physical grid / mesh with interval [t0, tf]
    std::shared_ptr<SpectralMesh> spectral_mesh;  // extended grid / mesh including the nominal interval [0, 1] (necessary for free t0, tf optimizations)
    Problem& problem;                             // continuous GDOP
    NLPState evaluation_state;                    // simple state to check which callbacks are performed for an iteration

    // scaling
    std::shared_ptr<ScalingFactory> scaling_factory;

    // initial guess
    std::unique_ptr<PrimalDualTrajectory> initial_guess; // set by initialization strategy

    // optimal solution
    std::unique_ptr<PrimalDualTrajectory> optimal_solution;

    // constant NLP derivative matrix part of the jacobian
    FixedVector<f64> const_der_jac;

    // offsets
    int off_x;         // offset #xVars
    int off_u;         // offset #uVars
    int off_p;         // offset #pVars
    int off_xu;        // number of vars for one collocation grid point
    int off_last_xu;   // last collocation grid point *x_nm, u_nm
    int off_xu_total;  // first parameter variable index
    int off_xup_total; // first time t0 variable index (tf += 1)
    int off_fg_total;  // first boundary constraint index
    int off_fgr_total; // first artificial initial control constraint index

    // note off_acc_xu[0][0] = off_xu for time t = first collocation node (as we also have xu(t0) before)
    FixedField<int, 2>   off_acc_xu; // offset to NLP_X first index of (x, u)(t_ij), i.e. NLP_X[off_acc_xu[i][j]] = x[i][j], u[i][j]
    FixedField<int, 2>   off_acc_fg; // offset to NLP_G first index of (f, g)(t_ij), i.e. NLP_G[off_acc_fg[i][j]] = f[i][j], g[i][j]
    FixedVector<int> off_acc_jac_fg; // offset to NLP_JAC_G first index of nabla (f, g)(t_ij)

    /* transforming dual variables (curr_obj_factor and lambda) passed to the callbacks
     * as callback is an Hessian we must do this ourselves before the callback
     * callback should already produce:
     *    ∑_{ij} curr_obj_factor * deltaT[i] * b[i][j] * L[i][j] + λ^T ∑_{ij} [-deltaT[i] * f_ij, g_ij] */

    /* scaling factors for lagrange terms in Hessian callback -
     * absorb "b[i][j] * delta_t[i]" in Lagrange term */
    FixedField<f64, 2> lagrange_obj_factors; // = curr_obj_factor * fLGR::b(mesh.intervals[i], mesh.nodes[j]) * mesh.delta_t[i]

    /* transformed dual variables passed to the callbacks - absorb "-deltaT[i]" in dynamics only */
    FixedVector<f64> transformed_lambda; // = lambda_NLP[i][j] * (- mesh.deltaT[i])

    // hessian sparsity helpers, O(1/2 * (x + u)² + p * (p + x + u)) memory, but no need for hashmaps, these are still fairly cheap
    // for further info see hessian layout at the bottom
    BlockSparsity hes_a_block, hes_b_block, hes_c_block, hes_d_block, hes_e_block, hes_f_block, hes_g_block, hes_h_block, hes_i_block, hes_k_block, hes_l_block;
    DenseRectangularBlockSparsity hes_j_block; // fully dense block
    OrderedIndexSet hes_A_set, hes_B_set, hes_C_set, hes_D_set, hes_E_set, hes_F_set, hes_G_set, hes_H_set, hes_I_set, hes_K_set, hes_L_set; // rename

    // === private / internal methods ===

    // inline methods for getting and providing current variable / dual addresses in callback
    // xu0 => xu(t0), xu => xu(t_01), xuf => xu(t_f), p => p, lamb_fg => fg(t_01), lamb_r => r
    inline const f64* get_x_xu0(const FixedVector<f64>& x) { return off_x        != 0 ?  x.raw()              : nullptr; }
    inline const f64* get_x_xu(const FixedVector<f64>& x)  { return off_xu       != 0 ? &x[off_xu]            : nullptr; }
    inline const f64* get_x_xuf(const FixedVector<f64>& x) { return off_xu       != 0 ? &x[off_last_xu]       : nullptr; }
    inline const f64* get_x_p(const FixedVector<f64>& x)   { return off_p        != 0 ? &x[off_xu_total]      : nullptr; }
    inline const f64* get_x_t0(const FixedVector<f64>& x)  { return spectral_mesh     ? &x[off_xup_total]     : nullptr; }
    inline const f64* get_x_tf(const FixedVector<f64>& x)  { return spectral_mesh     ? &x[off_xup_total + 1] : nullptr; }
    inline f64* get_lmbd_fg(FixedVector<f64>& lambda)      { return off_fg_total != 0 ?  lambda.raw()         : nullptr; }
    inline f64* get_lmbd_r(FixedVector<f64>& lambda) { return problem.pc->r_size != 0 ? &lambda[off_fg_total] : nullptr; }

    // helpers for initialize offsets 
    void create_acc_offset_xu(int off_xu);
    void create_acc_offset_fg(int off_fg);

    // init nlp and sparsity
    void init_jacobian_nonzeros(int& nnz_jac);
    void init_hessian_nonzeros(int& nnz_hes);

    // mutiply lambda (dual) with mesh factors => callbacks (except Lagrange) can use exact multipliers
    void update_curr_lambda_obj_factors(const FixedVector<f64>& curr_lambda, f64 curr_obj_factor);

    // hessian updates
    void accumulate_hessian_lfg(const HessianLFG& hes, int interval_i, int node_j,
                                const BlockSparsity& ptr_map_xu_xu, const BlockSparsity& ptr_map_p_xu,
                                FixedVector<f64>& curr_hes);
    void accumulate_hessian_parameter_lfg(const ParameterHessian& pp_hes, FixedVector<f64>& curr_hes); // sum of all weighted Hessian(Lfg)_pp
    void accumulate_hessian_mr(const HessianMR& hes, FixedVector<f64>& curr_hes);
    void accumulate_hessian_from_lagrangian_gradient_lf(int interval_i, int node_j, FixedVector<f64>& curr_hes,
                                                        const FixedVector<f64>& curr_lambda, f64 curr_obj_factor);
    void accumulate_hessian_from_lagrangian_gradient_lf_jac(int block_count, int interval_i, int node_j,
                                                            f64 factor_tf, const JacobianLFG& jac, FixedVector<f64>& curr_hes);
    // === callbacks to continuous problem (fill problem buffers) ===

    void callback_evaluation(const FixedVector<f64>& curr_x);
    void callback_jacobian(const FixedVector<f64>& curr_x);
    void callback_hessian(const FixedVector<f64> x, const FixedVector<f64>& curr_lambda, f64 curr_obj_factor);

    // === internal evaluations ===

    void check_new_x(bool new_x, const FixedVector<f64>& curr_x);
    void check_new_lambda(bool new_lambda);
    void eval_f_internal(f64& curr_obj);
    void eval_g_internal(const FixedVector<f64>& curr_x, FixedVector<f64>& curr_g);
    void eval_grad_f_internal(FixedVector<f64>& curr_grad);
    void eval_jac_g_internal(FixedVector<f64>& curr_jac);
    void eval_hes_internal(FixedVector<f64>& curr_hes, const FixedVector<f64>& curr_lambda, f64 curr_obj_factor);

    // === results + finalize + transform costates ===

    void flatten_trajectory_to_layout(
        const Trajectory& Trajectory,
        FixedVector<f64>& flat_buffer,
        bool from_costates);

    void transform_duals_costates(
        FixedVector<f64>& lambda,
        bool to_costate);

    void transform_duals_costates_bounds(
        FixedVector<f64>& zeta,
        bool to_costate);

    std::unique_ptr<Trajectory> finalize_optimal_primals(const FixedVector<f64>& opt_x);

    std::unique_ptr<CostateTrajectory> finalize_optimal_costates(const FixedVector<f64>& opt_lambda);

    std::pair<std::unique_ptr<Trajectory>, std::unique_ptr<Trajectory>> finalize_optimal_bound_duals(
        const FixedVector<f64>& opt_z_lb,
        const FixedVector<f64>& opt_z_ub);
};

} // namespace GDOP

#endif // MOO_GDOP_H

/*
Hessian Sparsity Layout (lower triangle):
    L: lower triangular matrix with diagonal
    X: square / rectangular matrix
    Note that blocks [[L, 0], [X, L]] are also triangular

                                            {n,m-1}
     | x00 u00 | x01 u01 | x02 u02 | x** u** | xnm1 unm1 | xnm unm |  p  | t0 tf |
---------------------------------------------------------------------------------|
 x00 |    L    |         |         |         |           |         |     |       |
 u00 |         |         |         |         |           |         |     |       |
---------------------------------------------------------------------------------|
 x01 |         |  L   0  |         |         |           |         |     |       |
 u01 |         |  X   L  |         |         |           |         |     |       |
 --------------------------------------------------------------------------------|
 x02 |         |         |  L   0  |         |           |         |     |       |
 u02 |         |         |  X   L  |         |           |         |     |       |
 --------------------------------------------------------------------------------|
 x** |         |         |         |  L   0  |           |         |     |       |
 u** |         |         |         |  X   L  |           |         |     |       |
 --------------------------------------------------------------------------------|
 xnm1|         |         |         |         |   L   0   |         |     |       |
 unm1|         |         |         |         |   X   L   |         |     |       |
---------------------------------------------------------------------------------|
 xnm |    X    |         |         |         |           |  L   0  |     |       |
 unm |    X    |         |         |         |           |  X   L  |     |       |
---------------------------------------------------------------------------------|
  p  |    X    |  X   X  |  X   X  |  X   X  |   X   X   |  X   X  |  L  |       |
---------------------------------------------------------------------------------|
  t0 |    X    |  X   X  |  X   X  |  X   X  |   X   X   |  X   X  |  X  |   L   |
  tf |    X    |  X   X  |  X   X  |  X   X  |   X   X   |  X   X  |  X  |       |
---------------------------------------------------------------------------------*

Block Sparsity Patterns: A - H
    where A=triang(x + u) B=triang(x + u), C=sq(x + u), D=triang(x + u), E=rect(p, x + u),
          F=rect(p, x + u), G=rect(p, x + u), H=triang(p, p)

Incl. SpectralMesh: I - K
    where I=rect(T, x + u), J=rect_dense(T, x + u), K=rect(T, p), L=triang(T, T) (sizeof(T) == 2)
    and   J are fully dense blocks
                                            {n,m-1}
     | x00 u00 | x01 u01 | x02 u02 | x** u** | xnm1 unm1 | xnm unm |  p  | t0 tf |
-------------------------------------------------------------------------|-------|
 x00 |    A    |         |         |         |           |         |     |       |
 u00 |         |         |         |         |           |         |     |       |
---------------------------------------------------------------------------------|
 x01 |         |    B    |         |         |           |         |     |       |
 u01 |         |         |         |         |           |         |     |       |
 --------------------------------------------------------------------------------|
 x02 |         |         |    B    |         |           |         |     |       |
 u02 |         |         |         |         |           |         |     |       |
 --------------------------------------------------------------------------------|
 x** |         |         |         |    B    |           |         |     |       |
 u** |         |         |         |         |           |         |     |       |
 --------------------------------------------------------------------------------|
 xnm1|         |         |         |         |     B     |         |     |       |
 unm1|         |         |         |         |           |         |     |       |
---------------------------------------------------------------------------------|
 xnm |    C    |         |         |         |           |    D    |     |       |
 unm |         |         |         |         |           |         |     |       |
---------------------------------------------------------------------------------|
  p  |    E    |    F    |    F    |    F    |     F     |    G    |  H  |       |
---------------------------------------------------------------------------------|
  t0 |    I    |    J    |    J    |    J    |     J     |    J    |  K  |   L   |
  tf |         |         |         |         |           |         |     |       |
---------------------------------------------------------------------------------*
*/
