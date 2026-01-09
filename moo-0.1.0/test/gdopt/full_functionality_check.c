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

#include <interfaces/c/structures.h>
#include <interfaces/gdopt/main_gdopt.h>
#include <generated.h>

// === problem sizes (compile const) ===

#define X_SIZE 1
#define U_SIZE 1
#define P_SIZE 1

#define RP_SIZE 1

#define R_SIZE 1
#define G_SIZE 1

#define HAS_MAYER true
#define HAS_LAGRANGE true


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

static bounds_t globl_x_bounds[X_SIZE] = { { -1, 1 } };
static bounds_t globl_u_bounds[U_SIZE] = { { -1, 1} };
static bounds_t globl_p_bounds[P_SIZE] = { { -1, 1} };
static bounds_t globl_T_bounds[2]      = { { -1.0, -0.5 }, { -0.2, 1.0 } };

static bounds_t globl_g_bounds[G_SIZE] = { {  -5, 5} };
static bounds_t globl_r_bounds[R_SIZE] = { { -10, 10}};

static optional_value_t globl_xu0_fixed[XU_SIZE] = { {0.0, false}, {0, false} };
static optional_value_t globl_xuf_fixed[XU_SIZE] = { {0.0, false}, {0, false} };

static f64 globl_x_nominal[X_SIZE];
static f64 globl_u_nominal[U_SIZE];
static f64 globl_p_nominal[P_SIZE];

static f64 globl_obj_nominal;
static f64 globl_f_nominal[X_SIZE];
static f64 globl_g_nominal[G_SIZE];
static f64 globl_r_nominal[R_SIZE];

static f64 globl_rp[RP_SIZE] = { 0.1 };

// === include data from csv-like files ===

static const char* data[FILE_COUNT]; // = { "inputpath.csv" };

// === optimization sparsity and evaluation structures (compile const) ===

static eval_structure_t globl_lfg_eval = {
    .buf_index = (int[]){0, 1, 2}
};

static coo_t globl_lfg_jac = {
    .row = (int[]){0, 0, 0, 1, 1, 1, 2, 2, 2},
    .col = (int[]){0, 1, 2, 0, 1, 2, 0, 1, 2},
    .buf_index = (int[]){0, 1, 2, 3, 4, 5, 6, 7, 8},
    .nnz = 9
};

static coo_t globl_lfg_lt_hes = {
    .row = (int[]){      0, 1, 1, 2, 2,                2},
    .col = (int[]){      0, 0, 1, 0, 1,                2},
    .buf_index = (int[]){0, 1, 2, 3, 4, /* pp block */ 0},
    .nnz = 6
};

static eval_structure_t globl_mr_eval = {
    .buf_index = (int[]){0, 1}
};

static coo_t globl_mr_jac = {
    .row = (int[]){0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
    .col = (int[]){0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6,},
    .buf_index = (int[]){0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
    .nnz = 14
};

static coo_t globl_mr_lt_hes = {
    .row = (int[]){0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6},
    .col = (int[]){0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6},
    .buf_index = (int[]){0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27},
    .nnz = 28
};

// === simulation sparsity & functions ===

static coo_t globl_ode_jac = {
    .row = (int[1]){0},
    .col = (int[1]){X_LFG_IDX(0)},
    .buf_index = (int[1]){0},
    .nnz = 1
};

static f64 lfg_arr[] = { 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9 };
static f64 mr_arr[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8 };

static inline int lt_index(int i, int j) {
    return i * (i + 1) / 2 + j;
}

static void ode_eval_f(const f64* x, const f64* u, const f64* p, f64 t,
                       const f64* data, f64* f, void* user_data)
{
    const int var_size = 3;
    const int func_idx = 1;
    const f64 v[] = { x[0], u[0], p[0] };

    int terms_per_func = var_size * (var_size + 1) / 2;
    const f64* b = lfg_arr + func_idx * terms_per_func;

    f64 sum = 0.0;
    for (int i = 0; i < var_size; i++) {
        int row_off = i * (i + 1) / 2;
        for (int j = 0; j <= i; j++) {
            f64 c = b[row_off + j];
            sum += c * v[i] * v[j];
        }
    }

    f[0] = sum;
}

static void ode_jac_f(const f64* x, const f64* u, const f64* p,
                      f64 t, const f64* data, f64* dfdx, void* user_data)
{
    const int func_idx = 1;
    const int col_idx  = 0;
    const int var_size = 3;
    const f64 v[] = { x[0], u[0], p[0] };

    const int terms_per_func = var_size * (var_size + 1) / 2;
    const f64* b = lfg_arr + func_idx * terms_per_func;

    f64 val = 0.0;
    for (int i = 0; i < var_size; i++) {
        int row_off = i * (i + 1) / 2;
        for (int j = 0; j <= i; j++) {
            f64 c = b[row_off + j];
            if (i == col_idx) val += c * v[j];
            if (j == col_idx) val += c * v[i];
        }
    }

    dfdx[0] = val;
}

// optimization functions

static void fill_eval(const f64* v, f64* out,
                      const f64* buf, int var_size, int func_size)
{
    int terms_per_func = var_size * (var_size + 1) / 2;
    for (int f = 0; f < func_size; f++) {
        const f64* b = buf + f * terms_per_func;
        f64 sum = 0.0;
        for (int i = 0; i < var_size; i++) {
            int row_off = i * (i + 1) / 2;
            for (int j = 0; j <= i; j++) {
                f64 c = b[row_off + j];
                sum += c * v[i] * v[j];
            }
        }
        out[f] = sum;
    }
}

static void fill_jac(const f64* v, f64* out,
                     const f64* buf, int var_size, int func_size)
{
    int terms_per_func = var_size * (var_size + 1) / 2;
    for (int f = 0; f < func_size; f++) {
        const f64* b = buf + f * terms_per_func;
        for (int k = 0; k < var_size; k++) {
            f64 val = 0.0;
            for (int i = 0; i < var_size; ++i) {
                int row_off = i * (i + 1) / 2;
                for (int j = 0; j <= i; j++) {
                    f64 c = b[row_off + j];
                    if (i == k) val += c * v[j];
                    if (j == k) val += c * v[i];
                }
            }
            out[f * var_size + k] = val;
        }
    }
}

static void fill_hes(const f64* v, const f64* lambda,
                     f64* out, f64* out_pp, const f64* buf,
                     int var_size, int func_size, int p_idx)
{
    int terms_per_func = var_size * (var_size + 1) / 2;
    int idx = 0;
    for (int k = 0; k < var_size; k++) {
        for (int l = 0; l <= k; l++) {
            f64 val = 0.0;
            for (int f = 0; f < func_size; f++) {
                const f64* b = buf + f * terms_per_func;
                int pos = k * (k + 1) / 2 + l;
                f64 c = b[pos];
                if (k == l) {
                    val += lambda[f] * (2.0 * c);
                } else {
                    val += lambda[f] * c;
                }
            }

            if (k == p_idx && l == p_idx) {
                if (out_pp) *out_pp += val;
            } else {
                out[idx++] = val;
            }
        }
    }
}

// [L, f, g]
static void eval_lfg(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data) {
    const f64* x = xu;
    const f64* u = xu + X_SIZE;
    f64 v[] = { x[0], u[0], p[0] };

    fill_eval(v, out, lfg_arr, 3, 3);
}

// ∇ [L, f, g]
static void jac_lfg(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data) {
    const f64* x = xu;
    const f64* u = xu + X_SIZE;
    f64 v[] = { x[0], u[0], p[0] };

    fill_jac(v, out, lfg_arr, 3, 3);
}

// σ ∇² L + λ^T ∇² [f, g] (lower triangle)
static void hes_lfg(const f64* xu, const f64* p, const f64* lambda, const f64 obj_factor, f64 t, const f64* data, f64* out, f64* out_pp, void* user_data) {
    const f64* x = xu;
    const f64* u = xu + X_SIZE;
    f64 v[] = { x[0], u[0], p[0] };
    f64 mu[] = { obj_factor, lambda[0], lambda[1] };

    fill_hes(v, mu, out, out_pp, lfg_arr, 3, 3, 2);
}

// [M, r]
static void eval_mr(const f64* xu0, const f64* xuf, const f64* p, f64 t0, f64 tf, const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {
    const f64* x0 = xu0;
    const f64* u0 = xu0 + X_SIZE;

    const f64* xf = xuf;
    const f64* uf = xuf + X_SIZE;

    f64 v[] = { x0[0], u0[0], xf[0], uf[0], p[0], t0, tf };

    fill_eval(v, out, mr_arr, 7, 2);
}

// ∇ [M, r]
static void jac_mr(const f64* xu0, const f64* xuf, const f64* p, f64 t0, f64 tf,
            const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {
    const f64* x0 = xu0;
    const f64* u0 = xu0 + X_SIZE;

    const f64* xf = xuf;
    const f64* uf = xuf + X_SIZE;

    f64 v[] = { x0[0], u0[0], xf[0], uf[0], p[0], t0, tf };

    fill_jac(v, out, mr_arr, 7, 2);
}

// σ ∇² M + λ^T ∇² r (lower triangle)
static void hes_mr(const f64* xu0, const f64* xuf, const f64* p, const f64* lambda, const f64 obj_factor, f64 t0, f64 tf,
            const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {
    const f64* x0 = xu0;
    const f64* u0 = xu0 + X_SIZE;

    const f64* xf = xuf;
    const f64* uf = xuf + X_SIZE;

    f64 mu[] = { obj_factor, lambda[0] };
    f64 v[] = { x0[0], u0[0], xf[0], uf[0], p[0], t0, tf };

    fill_hes(v, mu, out, (void*)0, mr_arr, 7, 2, -1);
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
    .initial_intervals = 3,
    .nodes_per_interval = 5,
    .l2bn_p1_it = 1,
    .l2bn_p2_it = 0,
    .l2bn_p2_lvl = 0.0
};

static solver_ctx_t globl_solver_ctx = {
    .derivative_test = true
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


int main_full_functionality(int argc, char** argv) {
    return main_gdopt(argc, argv, &globl_c_problem);
}
