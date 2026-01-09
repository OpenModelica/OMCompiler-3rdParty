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

/**
 * A simple DAE example
 *
 * @minimize tf
 *
 * s.t.   x' = -x + y + u
 *        0  = u² - sin(y) - x/4
 * fixed:
 *        x(t0) = 1.5
 *        x(tf) = 1.0
 *           t0 = -0.5
 *
 * variable:
 *     -inf <= x  <= inf
 *     -inf <= y  <= inf
 *     -1.0 <= u  <= 1.0
 *     -0.5 <= tf <= 5.0
 */

#include <interfaces/c/structures.h>
#include <interfaces/gdopt/main_gdopt.h>
#include <generated.h>

#include <math.h>

// === problem sizes (compile const) ===

#define X_SIZE 1
#define U_SIZE 2
#define P_SIZE 0

#define RP_SIZE 0

#define R_SIZE 0
#define G_SIZE 1

#define HAS_MAYER true
#define HAS_LAGRANGE false


#define FILE_COUNT 0

// === flat indices for sparsity

#define X_LFG_IDX(i)  (0 + (i))
#define U_LFG_IDX(i)  (X_SIZE + (i))
#define P_LFG_IDX(i)  (X_SIZE + U_SIZE + (i))

#define X0_MR_IDX(i)  (0 + (i))
#define U0_MR_IDX(i)  (X_SIZE + (i))
#define XF_MR_IDX(i)  (X_SIZE + U_SIZE + (i))
#define UF_MR_IDX(i)  (2 * X_SIZE + U_SIZE + (i))
#define P_MR_IDX(i)   (2 * (X_SIZE + U_SIZE) + (i))
#define T0_MR_IDX     (2 * (X_SIZE + U_SIZE) + P_SIZE)
#define TF_MR_IDX     (2 * (X_SIZE + U_SIZE) + P_SIZE + 1)

#define XU_SIZE (X_SIZE + U_SIZE)

// === declare global variables (values can be influenced by runtime parameters) ===

static bounds_t globl_x_bounds[X_SIZE] = { { -DBL_MAX, DBL_MAX } };
static bounds_t globl_u_bounds[U_SIZE] = { { -1, 1}, {-5, 5 } };
static bounds_t globl_p_bounds[P_SIZE];
static bounds_t globl_T_bounds[2]      = { { -0.5, -0.5 }, { -0.5, 5.0 } };

static bounds_t globl_g_bounds[G_SIZE] = { { 0, 0 } };
static bounds_t globl_r_bounds[R_SIZE];

static optional_value_t globl_xu0_fixed[XU_SIZE] = { {1.5, true}, {0, false}, {0, false} };
static optional_value_t globl_xuf_fixed[XU_SIZE] = { {1,   true}, {0, false}, {0, false} };

static f64 globl_x_nominal[X_SIZE];
static f64 globl_u_nominal[U_SIZE];
static f64 globl_p_nominal[P_SIZE];

static f64 globl_obj_nominal;
static f64 globl_f_nominal[X_SIZE];
static f64 globl_g_nominal[G_SIZE];
static f64 globl_r_nominal[R_SIZE];

static f64 globl_rp[RP_SIZE];

// === include data from csv-like files ===

static const char* data[FILE_COUNT]; // = { "inputpath.csv" };

// === optimization sparsity and evaluation structures (compile const) ===

static eval_structure_t globl_lfg_eval = {
    .buf_index = (int[]){0, 1}
};

static coo_t globl_lfg_jac = {
    .row = (int[]){0, 0, 0, 1, 1, 1},
    .col = (int[]){X_LFG_IDX(0), U_LFG_IDX(0), U_LFG_IDX(1), X_LFG_IDX(0), U_LFG_IDX(0), U_LFG_IDX(1)},
    .buf_index = (int[]){0, 1, 2, 3, 4, 5},
    .nnz = 6
};

static coo_t globl_lfg_lt_hes = {
    .row = (int[]){U_LFG_IDX(0), U_LFG_IDX(1)},
    .col = (int[]){U_LFG_IDX(0), U_LFG_IDX(1)},
    .buf_index = (int[]){0, 1},
    .nnz = 2
};

static eval_structure_t globl_mr_eval = {
    .buf_index = (int[]){0}
};

static coo_t globl_mr_jac = {
    .row = (int[]){0},
    .col = (int[]){TF_MR_IDX},
    .buf_index = (int[]){},
    .nnz = 1
};

static coo_t globl_mr_lt_hes = {
    .row = (int[]){},
    .col = (int[]){},
    .buf_index = (int[]){},
    .nnz = 0
};

// === simulation sparsity & functions ===

static coo_t globl_ode_jac = {
    .row = (int[1]){0},
    .col = (int[1]){X_LFG_IDX(0)},
    .buf_index = (int[1]){0},
    .nnz = 1
};

static void ode_eval_f(const f64* x, const f64* u, const f64* p, f64 t, const f64* data, f64* f, void* user_data) {
    f[0] = -x[0] + u[0] + u[1];
}

static void ode_jac_f(const f64* x, const f64* u, const f64* p, f64 t, const f64* data, f64* dfdx, void* user_data) {
    dfdx[0] = -1;
}

// optimization functions

// [L, f, g]
static void eval_lfg(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data) {
    const f64* x = xu;
    const f64* u = xu + X_SIZE;

    out[0] /* f */ = -x[0] + u[0] + u[1];
    out[1] /* g */ = u[0] * u[0] - sin(u[1]) - 0.25 * x[0];
}

// ∇ [L, f, g]
static void jac_lfg(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data) {
    const f64* x = xu;
    const f64* u = xu + X_SIZE;

    out[0] = -1; /* f_x  */
    out[1] = 1; /* f_u */
    out[2] = 1; /* f_y */
    out[3] = -0.25; /* g_x */
    out[4] = 2 * u[0]; /* g_u */
    out[5] = -cos(u[1]); /* g_y */
}

// σ ∇² L + λ^T ∇² [f, g] (lower triangle)
static void hes_lfg(const f64* xu, const f64* p, const f64* lambda, const f64 obj_factor, f64 t, const f64* data, f64* out, f64* out_pp, void* user_data) {
    const f64* x = xu;
    const f64* u = xu + X_SIZE;

    out[0] = 2 * lambda[1];
    out[1] = lambda[1] * sin(u[1]);
}

// [M, r]
static void eval_mr(const f64* xu0, const f64* xuf, const f64* p, f64 t0, f64 tf, const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {
    const f64* xf = xuf;
    const f64* uf = xuf + X_SIZE;

    out[0] = tf;
}

// ∇ [M, r]
static void jac_mr(const f64* xu0, const f64* xuf, const f64* p, f64 t0, f64 tf,
            const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {
    out[0] = 1.0;
}

// σ ∇² M + λ^T ∇² r (lower triangle)
static void hes_mr(const f64* xu0, const f64* xuf, const f64* p, const f64* lambda, const f64 obj_factor, f64 t0, f64 tf,
            const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {

}

// === objects ===

static c_callbacks_t globl_callbacks = {
    eval_lfg,
    jac_lfg,
    hes_lfg,
    eval_mr,
    jac_mr,
    hes_mr,

    ode_eval_f,
    ode_jac_f
};

static mesh_ref_ctx_t globl_mesh_ctx = {
    .initial_intervals = 10,
    .nodes_per_interval = 5,
    .l2bn_p1_it = 1,
    .l2bn_p2_it = 0,
    .l2bn_p2_lvl = 0.0
};

static solver_ctx_t globl_solver_ctx = {
    .derivative_test = false
};

static c_problem_t globl_c_problem = {
    .x_size = X_SIZE,
    .u_size = U_SIZE,
    .xu_size = X_SIZE + U_SIZE,
    .p_size = P_SIZE,
    .rp_size = RP_SIZE,
    .r_size = R_SIZE,
    .g_size = G_SIZE,
    .has_mayer = HAS_MAYER,
    .has_lagrange = HAS_LAGRANGE,
    .data_filepath = data,
    .data_file_count = FILE_COUNT,
    .rp = globl_rp,
    .x_bounds = globl_x_bounds,
    .u_bounds = globl_u_bounds,
    .p_bounds = globl_p_bounds,
    .T_bounds = globl_T_bounds,
    .r_bounds = globl_r_bounds,
    .g_bounds = globl_g_bounds,
    .xu0_fixed = globl_xu0_fixed,
    .xuf_fixed = globl_xuf_fixed,
    .x_nominal = globl_x_nominal,
    .u_nominal = globl_u_nominal,
    .p_nominal = globl_p_nominal,
    .obj_nominal = &globl_obj_nominal,
    .f_nominal = globl_f_nominal,
    .g_nominal = globl_g_nominal,
    .r_nominal = globl_r_nominal,
    .lfg_eval  = &globl_lfg_eval,
    .lfg_jac  = &globl_lfg_jac,
    .lfg_lt_hes  = &globl_lfg_lt_hes,
    .mr_eval = &globl_mr_eval,
    .mr_jac = &globl_mr_jac,
    .mr_lt_hes = &globl_mr_lt_hes,
    .ode_jac = &globl_ode_jac,
    .callbacks = &globl_callbacks,
    .mesh_ctx = &globl_mesh_ctx,
    .solver_ctx = &globl_solver_ctx,
    .user_data = (void*)0
};


int main_sanity_check(int argc, char** argv) {
    return main_gdopt(argc, argv, &globl_c_problem);
}
