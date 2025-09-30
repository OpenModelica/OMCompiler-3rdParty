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

#define R_SIZE 0
#define G_SIZE 0

#define HAS_MAYER false
#define HAS_LAGRANGE true


#define FILE_COUNT 0

// === declare global variables (values can be influenced by runtime parameters) ===

bounds_t globl_x_bounds[X_SIZE] = { { -DBL_MAX, DBL_MAX } };
bounds_t globl_u_bounds[U_SIZE] = { { -DBL_MAX, DBL_MAX } };
bounds_t globl_p_bounds[P_SIZE] = { { -DBL_MAX, DBL_MAX } };

bounds_t globl_g_bounds[G_SIZE];
bounds_t globl_r_bounds[R_SIZE];

optional_value_t globl_x0_fixed[X_SIZE] = { {1.5, true} };
optional_value_t globl_xf_fixed[X_SIZE] = { {1,   true} };

f64 globl_x_nominal[X_SIZE];
f64 globl_u_nominal[U_SIZE];
f64 globl_p_nominal[P_SIZE];

f64 globl_obj_nominal;
f64 globl_f_nominal[X_SIZE];
f64 globl_g_nominal[G_SIZE];
f64 globl_r_nominal[R_SIZE];

f64 globl_rp[RP_SIZE];

// TODO: include also csv for file init?

// === include data from csv-like files ===

const char* data[FILE_COUNT]; // = { "inputpath.csv" };

/* data is available in the callback functions as f64* data (data at time of evaluation)
   the control columns of the csvs are extracted and sorted as a flat array
        [CSV1: [col1, col2, ..., colN], CSV2: [col1, col2, ..., colM], ...]
   the constant parameters (static section) are written to the passed runtime parameter array
   these do not change over time

 example CSV like format

# time: 0                      // time column = 0.00000000000000000, ...
# dynamic: "x", [1]            // dynamic (time-varying) column for states = 1.50000000000000000, ...
# dynamic: "u", [2]            // dynamic (time-varying) column for controls = -0.62864318974690347, ...
# static: "p", [3]             // static (time-varying) column for parameters = 0.02500176076771821
t,x_0,u_0,p_0
0.00000000000000000,1.50000000000000000,-0.62864318974690347,0.02500176076771821
0.01445750910019651,1.46989539310990791,-0.61617345328311646
0.07616008346272038,1.34811292027240026,-0.56572950137448774

*/

// === optimization sparsity and evaluation structures (compile const) ===

eval_structure_t globl_lfg_eval = {
    .buf_index = (int[]){0, 1}
};

coo_t globl_lfg_jac = {
    .row = (int[]){0, 0, 1, 1, 1},
    .col = (int[]){0, 1, 0, 1, 2},
    .buf_index = (int[]){0, 1, 2, 3, 4},
    .nnz = 5
};

coo_t globl_lfg_lt_hes = {
    .row = (int[]){0, 1},
    .col = (int[]){0, 1},
    .buf_index = (int[]){0, 1},
    .nnz = 2
};

eval_structure_t globl_mr_eval = {
    .buf_index = (int[]){}
};

coo_t globl_mr_jac = {
    .row = (int[]){},
    .col = (int[]){},
    .buf_index = (int[]){},
    .nnz = 0
};

coo_t globl_mr_lt_hes = {
    .row = (int[]){},
    .col = (int[]){},
    .buf_index = (int[]){},
    .nnz = 0
};

// === simulation sparsity & functions ===

coo_t globl_ode_jac = {
    .row = (int[1]){0},
    .col = (int[1]){0},
    .buf_index = (int[1]){0},
    .nnz = 1
};

void ode_eval_f(const f64* x, const f64* u, const f64* p, f64 t, const f64* data, f64* f, void* user_data) {
    f[0] = -x[0] + u[0] + p[0];
}

void ode_jac_f(const f64* x, const f64* u, const f64* p, f64 t, const f64* data, f64* dfdx, void* user_data) {
    dfdx[0] = -1;
}

// optimization functions

// [L, f, g]
void eval_lfg(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data) {
    const f64* x = xu;
    const f64* u = xu + X_SIZE;

    out[0] /* L */ = 0.5 * (x[0] * x[0] + u[0] * u[0]);
    out[1] /* f */ = -x[0] + u[0] + p[0];
}

// ∇ [L, f, g]
void jac_lfg(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data) {
    const f64* x = xu;
    const f64* u = xu + X_SIZE;

    out[0] = x[0]; /* L_x */
    out[1] = u[0]; /* L_u */
    out[2] = -1; /* f_x  */
    out[3] = 1; /* f_u */
    out[4] = 1; /* f_p */
}

// σ ∇² L + λ^T ∇² [f, g] (lower triangle)
void hes_lfg(const f64* xu, const f64* p, const f64* lambda, const f64 obj_factor, f64 t, const f64* data, f64* out, void* user_data) {
    out[0] = obj_factor; /* xx */
    out[1] = obj_factor; /* uu */
}

// [M, r]
void eval_mr(const f64* x0, const f64* xuf, const f64* p, f64 t0, f64 tf, const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {
    const f64* xf = xuf;
    const f64* uf = xuf + X_SIZE;

}

// ∇ [M, r]
void jac_mr(const f64* x0, const f64* xuf, const f64* p, f64 t0, f64 tf,
            const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {
}

// σ ∇² M + λ^T ∇² r (lower triangle)
void hes_mr(const f64* x0, const f64* xuf, const f64* p, const f64* lambda, const f64 obj_factor, f64 t0, f64 tf,
            const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {

}

// === objects ===

c_callbacks_t globl_callbacks = {
    eval_lfg,
    jac_lfg,
    hes_lfg,
    eval_mr,
    jac_mr,
    hes_mr,

    ode_eval_f,
    ode_jac_f
};

c_problem_t globl_c_problem = {
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
    .r_bounds = globl_r_bounds,
    .g_bounds = globl_g_bounds,
    .x0_fixed = globl_x0_fixed,
    .xf_fixed = globl_xf_fixed,
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
    .user_data = (void*)0
};


int main_generated(int argc, char** argv) {
    main_gdopt(argc, argv, &globl_c_problem);
    return 0;
}
