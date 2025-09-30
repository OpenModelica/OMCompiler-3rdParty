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
#ifndef MOO_C_STRUCTS_H
#define MOO_C_STRUCTS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <float.h>

typedef double f64;

typedef struct bounds_t {
    f64 lb;
    f64 ub;
} bounds_t;

typedef struct optional_value_t {
    f64 value;
    bool is_set;
} optional_value_t;

// eval_structure_t[fn] -> buf_index of out in eval_lfg() - fn is sorted as L -> f -> g
typedef struct eval_structure_t {
    int* buf_index;  // buf_index
} eval_structure_t;

typedef struct coo_t {
    int* row;        // row indices
    int* col;        // col indices
    int* buf_index;  // non-zero index == buf_index
    int nnz;         // total nnz
} coo_t;

typedef struct c_callbacks_t {
    void (*eval_lfg)(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data);
    void (*jac_lfg)(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data);
    void (*hes_lfg)(const f64* xu, const f64* p, const f64* lambda, const f64 obj_factor, f64 t, const f64* data, f64* out, void* user_data);
    void (*eval_mr)(const f64* x0, const f64* xuf, const f64* p, f64 t0, f64 tf,
                    const f64* data_t0, const f64* data_tf, f64* out, void* user_data);
    void (*jac_mr)(const f64* x0, const f64* xuf, const f64* p, f64 t0, f64 tf,
                   const f64* data_t0, const f64* data_tf, f64* out, void* user_data);
    void (*hes_mr)(const f64* x0, const f64* xuf, const f64* p, const f64* lambda, const f64 obj_factor, f64 t0, f64 tf,
                   const f64* data_t0, const f64* data_tf, f64* out, void* user_data);

    void (*ode_f)(const f64* x, const f64* u, const f64* p, f64 t, const f64* data, f64* f, void* user_data);
    void (*ode_jac_f)(const f64* x, const f64* u, const f64* p, f64 t, const f64* data, f64* dfdx, void* user_data) ;
} c_callbacks_t;

typedef struct c_problem_t {
    const int x_size;
    const int u_size;
    const int xu_size;
    const int p_size;
    const int rp_size;

    const int r_size;
    const int g_size;

    const bool has_mayer;
    const bool has_lagrange;

    const char** data_filepath;
    const int data_file_count;
    f64* rp;

    bounds_t* x_bounds;
    bounds_t* u_bounds;
    bounds_t* p_bounds;

    bounds_t* r_bounds;
    bounds_t* g_bounds;

    optional_value_t* x0_fixed;
    optional_value_t* xf_fixed;

    f64* x_nominal;
    f64* u_nominal;
    f64* p_nominal;

    f64* obj_nominal;
    f64* f_nominal;
    f64* g_nominal;
    f64* r_nominal;

    eval_structure_t* lfg_eval;
    coo_t* lfg_jac;
    coo_t* lfg_lt_hes;

    eval_structure_t* mr_eval;
    coo_t* mr_jac;
    coo_t* mr_lt_hes;

    // TODO: this ignores the buf_idx of coo_t - do smth about it?
    coo_t* ode_jac;

    c_callbacks_t* callbacks;
    void* user_data;

// private

    // data read from csvs -> automatically passed to callbacks
    f64* data;
    int data_chunk_size;
} c_problem_t;

#ifdef __cplusplus
}
#endif

#endif // MOO_C_STRUCTS_H
