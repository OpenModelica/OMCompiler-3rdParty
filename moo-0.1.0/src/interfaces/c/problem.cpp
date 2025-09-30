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

#include "problem.h"

namespace C {

FixedVector<Bounds> create_bounds(bounds_t* c_bounds, int size) {
    FixedVector<Bounds> to_fill(size);
    for (int idx = 0; idx < size; idx++) {
        to_fill[idx].lb = c_bounds[idx].lb;
        to_fill[idx].ub = c_bounds[idx].ub;
    }
    return to_fill;
}

FixedVector<std::optional<f64>> create_fixed(optional_value_t* c_optional, int size) {
    FixedVector<std::optional<f64>> to_fill(size);
    for (int idx = 0; idx < size; idx++) {
        if (c_optional[idx].is_set) {
            to_fill[idx] = c_optional[idx].value;
        }
    }
    return to_fill;
}

FixedVector<f64> assign_or_one(f64* c_array, int size) {
    FixedVector<f64> to_fill(size);
    if (!c_array) {
        for (int idx = 0; idx < size; idx++) { c_array[idx] = 1; }
    } else {
        to_fill.assign(c_array);
    }
    return to_fill;
}

void layout_lfg_init_eval(GDOP::FullSweepLayout& layout_lfg, c_problem_t* c_problem) {
    if (layout_lfg.L) {
        layout_lfg.L->buf_index = c_problem->lfg_eval->buf_index[0];
    }

    int offset = (layout_lfg.L ? 1 : 0);
    for (int f_idx = 0; f_idx < c_problem->x_size; f_idx++) {
        layout_lfg.f[f_idx].buf_index = c_problem->lfg_eval->buf_index[offset + f_idx];
    }

    offset += c_problem->x_size;
    for (int g_idx = 0; g_idx < c_problem->g_size; g_idx++) {
        layout_lfg.g[g_idx].buf_index = c_problem->lfg_eval->buf_index[offset + g_idx];
    }
}

void layout_lfg_init_jac(GDOP::FullSweepLayout& layout_lfg, c_problem_t* c_problem){
    for (int nz = 0; nz < c_problem->lfg_jac->nnz; nz++) {
        int row = c_problem->lfg_jac->row[nz];
        int col = c_problem->lfg_jac->col[nz];
        int buf_index = c_problem->lfg_jac->buf_index[nz];

        auto& fn = GDOP::access_Lfg_from_row(layout_lfg, row);
        if (col < c_problem->x_size) {
            fn.jac.dx.push_back(JacobianSparsity{col, buf_index});
        }
        else if (col < c_problem->xu_size) {
            fn.jac.du.push_back(JacobianSparsity{col - c_problem->x_size, buf_index});
        }
        else {
            fn.jac.dp.push_back(JacobianSparsity{col - c_problem->xu_size, buf_index});
        }
    }
}

void layout_lfg_init_hes(GDOP::FullSweepLayout& layout_lfg, c_problem_t* c_problem){
    for (int nz = 0; nz < c_problem->lfg_lt_hes->nnz; nz++) {
        int row = c_problem->lfg_lt_hes->row[nz];
        int col = c_problem->lfg_lt_hes->col[nz];
        int buf_index = c_problem->lfg_lt_hes->buf_index[nz];
        auto& hes = layout_lfg.hes;
        auto& hes_pp = layout_lfg.pp_hes;

        assert(row >= col && "Hessian must be supplied in lower triangular form: row >= col index.");

        if (row < c_problem->x_size && col < c_problem->x_size) {
            hes.dx_dx.push_back({row, col, buf_index});
        }
        else if (row < c_problem->xu_size && col < c_problem->x_size) {
            hes.du_dx.push_back({row - c_problem->x_size, col, buf_index});
        }
        else if (row < c_problem->xu_size && col < c_problem->xu_size) {
            hes.du_du.push_back({row - c_problem->x_size, col - c_problem->x_size, buf_index});
        }
        else if (col < c_problem->x_size) {
            hes.dp_dx.push_back({row - c_problem->xu_size, col, buf_index});
        }
        else if (col < c_problem->xu_size) {
            hes.dp_du.push_back({row - c_problem->xu_size, col - c_problem->x_size, buf_index});
        }
        else {
            hes_pp.dp_dp.push_back({row - c_problem->xu_size, col - c_problem->xu_size, buf_index});
        }
    }
}

GDOP::FullSweepLayout create_fullsweep_layout(c_problem_t* c_problem) {
    auto layout_lfg = GDOP::FullSweepLayout(c_problem->has_lagrange, c_problem->x_size, c_problem->g_size);

    layout_lfg_init_eval(layout_lfg, c_problem);
    layout_lfg_init_jac(layout_lfg, c_problem);
    layout_lfg_init_hes(layout_lfg, c_problem);

    return layout_lfg;
}

void layout_mr_init_eval(GDOP::BoundarySweepLayout& layout_mr, c_problem_t* c_problem) {
    if (layout_mr.M) {
        layout_mr.M->buf_index = c_problem->mr_eval->buf_index[0];
    }

    int offset = (layout_mr.M ? 1 : 0);
    for (int r_idx = 0; r_idx < c_problem->r_size; r_idx++) {
        layout_mr.r[r_idx].buf_index = c_problem->mr_eval->buf_index[offset + r_idx];
    }
}

void layout_mr_init_jac(GDOP::BoundarySweepLayout& layout_mr, c_problem_t* c_problem){
    for (int nz = 0; nz < c_problem->mr_jac->nnz; nz++) {
        int row = c_problem->mr_jac->row[nz];
        int col = c_problem->mr_jac->col[nz];
        int buf_index = c_problem->mr_jac->buf_index[nz];

        auto& fn = (layout_mr.M && row == 0 ? *layout_mr.M : layout_mr.r[row - (layout_mr.M ? 1 : 0)]);
        if (col < c_problem->x_size) {
            fn.jac.dx0.push_back(JacobianSparsity{col, buf_index});
        }
        else if (col < 2 * c_problem->x_size) {
            fn.jac.dxf.push_back(JacobianSparsity{col - c_problem->x_size, buf_index});
        }
        else if (col < 2 * c_problem->x_size + c_problem->u_size) {
            fn.jac.duf.push_back(JacobianSparsity{col - 2 * c_problem->x_size, buf_index});
        }
        else {
            fn.jac.dp.push_back(JacobianSparsity{col - (2 * c_problem->x_size + c_problem->u_size), buf_index});
        }
    }
}

void layout_mr_init_hes(GDOP::BoundarySweepLayout& layout_mr, c_problem_t* c_problem){
    for (int nz = 0; nz < c_problem->lfg_lt_hes->nnz; nz++) {
        int row = c_problem->lfg_lt_hes->row[nz];
        int col = c_problem->lfg_lt_hes->col[nz];
        int buf_index = c_problem->lfg_lt_hes->buf_index[nz];
        auto& hes = layout_mr.hes;

        assert(row >= col && " Hessian must be supplied in lower triangular form: row >= col index.");

        if (row < c_problem->x_size && col < c_problem->x_size) {
            hes.dx0_dx0.push_back({row, col, buf_index});
        }
        else if (row < 2 * c_problem->x_size && col < c_problem->x_size) {
            hes.dxf_dx0.push_back({row - c_problem->x_size, col, buf_index});
        }
        else if (row < 2 * c_problem->x_size + c_problem->u_size && col < c_problem->x_size) {
            hes.duf_dx0.push_back({row - 2 * c_problem->x_size, col, buf_index});
        }
        else if (col < c_problem->x_size) {
            hes.dp_dx0.push_back({row - (2 * c_problem->x_size + c_problem->u_size), col, buf_index});
        }
        else if (row < 2 * c_problem->x_size && col < 2 * c_problem->x_size) {
            hes.dxf_dxf.push_back({row - c_problem->x_size, col - c_problem->x_size, buf_index});
        }
        else if (row < 2 * c_problem->x_size + c_problem->u_size && col < 2 * c_problem->x_size) {
            hes.duf_dxf.push_back({row - 2 * c_problem->x_size, col - c_problem->x_size, buf_index});
        }
        else if (col < 2 * c_problem->x_size) {
            hes.dp_dxf.push_back({row - (2 * c_problem->x_size + c_problem->u_size), col - c_problem->x_size, buf_index});
        }
        else if (row < 2 * c_problem->x_size + c_problem->u_size && col < 2 * c_problem->x_size + c_problem->u_size) {
            hes.duf_duf.push_back({row - 2 * c_problem->x_size, col - 2 * c_problem->x_size, buf_index});
        }
        else if (col < 2 * c_problem->x_size + c_problem->u_size) {
            hes.dp_duf.push_back({row - (2 * c_problem->x_size + c_problem->u_size), col - 2 * c_problem->x_size, buf_index});
        }
        else {
            hes.dp_dp.push_back({row - (2 * c_problem->x_size + c_problem->u_size), col - (2 * c_problem->x_size + c_problem->u_size), buf_index});
        }
    }
}

GDOP::BoundarySweepLayout create_boundarysweep_layout(c_problem_t* c_problem) {
    auto layout_mr = GDOP::BoundarySweepLayout(c_problem->has_mayer, c_problem->r_size);

    layout_mr_init_eval(layout_mr, c_problem);
    layout_mr_init_jac(layout_mr, c_problem);
    layout_mr_init_hes(layout_mr, c_problem);

    return layout_mr;
}

void FullSweep::callback_eval(
    const f64* xu_nlp,
    const f64* p)
{
    for (int i = 0; i < pc.mesh->intervals; i++) {
        for (int j = 0; j < pc.mesh->nodes[i]; j++) {
            const f64* xu_ij = get_xu_ij(xu_nlp, i, j);
            const f64* data_ij = get_data_ij(i, j);
            f64* eval_buf_ij = get_eval_buffer(i, j);
            c_callbacks->eval_lfg(xu_ij, p, pc.mesh->t[i][j], data_ij, eval_buf_ij, c_problem->user_data);
        }
    }
}

void FullSweep::callback_jac(
    const f64* xu_nlp,
    const f64* p)
{
    for (int i = 0; i < pc.mesh->intervals; i++) {
        for (int j = 0; j < pc.mesh->nodes[i]; j++) {
            const f64* xu_ij = get_xu_ij(xu_nlp, i, j);
            const f64* data_ij = get_data_ij(i, j);
            f64* jac_buf_ij = get_jac_buffer(i, j);
            c_callbacks->jac_lfg(xu_ij, p, pc.mesh->t[i][j], data_ij, jac_buf_ij, c_problem->user_data);
        }
    }
}

void FullSweep::callback_hes(
    const f64* xu_nlp,
    const f64* p,
    const FixedField<f64, 2>& lagrange_factors,
    const f64* lambda)
{
    for (int i = 0; i < pc.mesh->intervals; i++) {
        for (int j = 0; j < pc.mesh->nodes[i]; j++) {
            const f64* xu_ij = get_xu_ij(xu_nlp, i, j);
            const f64* data_ij = get_data_ij(i, j);
            const f64* lmbd_ij = get_lambda_ij(lambda, i, j);
            f64* hes_buffer = get_hes_buffer(i, j);
            c_callbacks->hes_lfg(xu_ij, p, lmbd_ij, lagrange_factors[i][j], pc.mesh->t[i][j], data_ij, hes_buffer, c_problem->user_data);
        }
    }
}

void BoundarySweep::callback_eval(
    const f64* x0_nlp,
    const f64* xuf_nlp,
    const f64* p)
{
    const f64* data_t0 = get_data_t0();
    const f64* data_tf = get_data_tf();
    f64* eval_buf = get_eval_buffer();
    c_callbacks->eval_mr(x0_nlp, xuf_nlp, p, 0, pc.mesh->tf, data_t0, data_tf, eval_buf, c_problem->user_data);
};

void BoundarySweep::callback_jac(
    const f64* x0_nlp,
    const f64* xuf_nlp,
    const f64* p)
{
    const f64* data_t0 = get_data_t0();
    const f64* data_tf = get_data_tf();
    f64* jac_buf = get_jac_buffer();
    c_callbacks->jac_mr(x0_nlp, xuf_nlp, p, 0, pc.mesh->tf, data_t0, data_tf, jac_buf, c_problem->user_data);
};

void BoundarySweep::callback_hes(
    const f64* x0_nlp,
    const f64* xuf_nlp,
    const f64* p,
    const f64 mayer_factor,
    const f64* lambda)
{
    const f64* data_t0 = get_data_t0();
    const f64* data_tf = get_data_tf();
    f64* hes_buf = get_hes_buffer();
    c_callbacks->hes_mr(x0_nlp, xuf_nlp, p, lambda, mayer_factor, 0, pc.mesh->tf, data_t0, data_tf, hes_buf, c_problem->user_data);
};

f64* Dynamics::get_data(f64 t) {
    if (c_problem->data_file_count == 0) {
        return nullptr;
    }

    f64* ret_ptr = current_data.raw();
    f64* work_ptr = ret_ptr;
    if (t != current_data_time) {
        // new interpolation -> write to contiguous `current_data` buffer
        for (int file_idx = 0; file_idx < c_problem->data_file_count; file_idx++) {
            auto& ctrl = raw_ctrl_data[file_idx];
            ctrl.interpolate_at(t, work_ptr);
            work_ptr += ctrl.u.size();
        }

        current_data_time = t;
    }

    return ret_ptr;
}

/**
 * @brief Get the sum of all sizes of the user-provided data trajectories.
 */
size_t accumulate_data_size(c_problem_t* c_problem, std::unique_ptr<ControlTrajectory[]>& raw_ctrl_data) {
    size_t size = 0;
    for (int file_idx = 0; file_idx < c_problem->data_file_count; file_idx++) {
        size += raw_ctrl_data[file_idx].u.size();
    }
    return size;
}

/**
 * @brief Extract the controls of the raw data trajectories. Allows for interpolation using builtin
 *        ControlTrajectory utils.
 */
std::unique_ptr<ControlTrajectory[]> extract_ctrl_data(c_problem_t* c_problem, std::shared_ptr<Trajectory[]>& raw_data) {
    auto ctrl = std::make_unique<ControlTrajectory[]>(c_problem->data_file_count);
    for (int file_idx = 0; file_idx < c_problem->data_file_count; file_idx++) {
        ctrl[file_idx] = raw_data[file_idx].copy_extract_controls();
    }
    return ctrl;
}

Dynamics::Dynamics(const GDOP::ProblemConstants& pc,
                   c_problem_t* c_problem_,
                   std::shared_ptr<Trajectory[]> raw_data)
    : GDOP::Dynamics(pc, ::Simulation::Jacobian::sparse(
                             ::Simulation::JacobianFormat::COO,
                             c_problem_->ode_jac->row,
                             c_problem_->ode_jac->col,
                             c_problem_->ode_jac->nnz)
                    ),
      c_callbacks(c_problem_->callbacks),
      c_problem(c_problem_),
      raw_ctrl_data(extract_ctrl_data(c_problem, raw_data)),
      current_data(FixedVector<f64>(accumulate_data_size(c_problem, raw_ctrl_data))),
      current_data_time(MINUS_INFINITY) {};

void Dynamics::eval(const f64* x, const f64* u, const f64* p, f64 t, f64* f, void* user_data) {
    c_callbacks->ode_f(x, u, p, t, get_data(t), f, user_data);
}

void Dynamics::jac(const f64* x, const f64* u, const f64* p, f64 t, f64* dfdx, void* user_data) {
    c_callbacks->ode_jac_f(x, u, p, t, get_data(t), dfdx, user_data);
}

GDOP::Problem create_gdop_problem(c_problem_t* c_problem, const Mesh& mesh, std::shared_ptr<Trajectory[]> raw_data) {
    auto x_bounds = create_bounds(c_problem->x_bounds, c_problem->x_size);
    auto u_bounds = create_bounds(c_problem->u_bounds, c_problem->u_size);
    auto p_bounds = create_bounds(c_problem->p_bounds, c_problem->p_size);
    auto g_bounds = create_bounds(c_problem->g_bounds, c_problem->g_size);
    auto r_bounds = create_bounds(c_problem->r_bounds, c_problem->r_size);

    auto x0_fixed = create_fixed(c_problem->x0_fixed, c_problem->x_size);
    auto xf_fixed = create_fixed(c_problem->xf_fixed, c_problem->x_size);

    auto pc = std::make_unique<GDOP::ProblemConstants>(
        c_problem->has_mayer,
        c_problem->has_lagrange,
        std::move(x_bounds),
        std::move(u_bounds),
        std::move(p_bounds),
        std::move(x0_fixed),
        std::move(xf_fixed),
        std::move(r_bounds),
        std::move(g_bounds),
        mesh
    );

    auto fs  = std::make_unique<FullSweep>(create_fullsweep_layout(c_problem), *pc, c_problem);
    auto bs  = std::make_unique<BoundarySweep>(create_boundarysweep_layout(c_problem), *pc, c_problem);
    auto dyn = std::make_unique<Dynamics>(*pc, c_problem, raw_data);

    return GDOP::Problem(std::move(fs), std::move(bs),std::move(pc), std::move(dyn));
}

/**
 * @brief Read in all trajectories provided by the user
 */
std::shared_ptr<Trajectory[]> create_raw_data(c_problem_t* c_problem) {
    if (c_problem->data_file_count == 0) return nullptr;

    std::shared_ptr<Trajectory[]> raw_data(new Trajectory[c_problem->data_file_count],
                                           std::default_delete<Trajectory[]>());

    for (int i = 0; i < c_problem->data_file_count; i++) {
        raw_data[i] = Trajectory::from_csv(std::string(c_problem->data_filepath[i]));
    }

    return raw_data;
}

/**
 * @brief Extract the parameters of all input files and write them into the
 *        user-provided contiguous runtime parameter buffer
 */
void Problem::fill_runtime_parameters() {
    size_t offset = 0;
    for (int file_idx = 0; file_idx < c_problem->data_file_count; file_idx++) {
        auto const& raw_rp = raw_data[file_idx].p;
        if (offset + raw_rp.size() > static_cast<size_t>(c_problem->rp_size)) {
            Log::error("Runtime parameters out of range.");
            abort();
        }

        std::memcpy(c_problem->rp + offset, raw_rp.data(), raw_rp.size() * sizeof(f64));
        offset += raw_rp.size();
    }
}

/**
 * @brief Interpolate the raw data trajectories onto a given mesh (interpolated_data),
 *        extract the controls and write them into c_problem->data as contiguous memory for callbacks.
 *        -> No more interpolation needed during the optimization callbacks / only lookup.
 */
void Problem::fill_data(const Mesh& mesh) {
    int block_size = 0;
    for (int file_idx = 0; file_idx < c_problem->data_file_count; file_idx++) {
        interpolated_data[file_idx] = raw_data[file_idx].interpolate_onto_mesh(mesh);
        block_size += static_cast<int>(raw_data[file_idx].u.size());
    }

    if (block_size == 0) {
        c_problem->data = nullptr;
        c_problem->data_chunk_size = 0;
        return;
    }

    flat_interpolated_input_data = FixedVector<f64>(block_size * (mesh.node_count + 1));

    c_problem->data_chunk_size = block_size;
    c_problem->data = flat_interpolated_input_data.raw(); // ref to flat_interpolated_input_data

    size_t offset = 0;
    for (int t_idx = 0; t_idx < mesh.node_count + 1; t_idx++) {
        for (int file_idx = 0; file_idx < c_problem->data_file_count; file_idx++) {
            auto& u = interpolated_data[file_idx].u;
            for (size_t u_idx = 0; u_idx < u.size(); u_idx++) {
                flat_interpolated_input_data[offset + u_idx] = u[u_idx][t_idx];
            }
            offset += u.size();
        }
    }
}

Problem::Problem(c_problem_t* c_problem, const Mesh& mesh, std::shared_ptr<Trajectory[]> raw_data)
    : GDOP::Problem(create_gdop_problem(c_problem, mesh, raw_data)),
      c_callbacks(c_problem->callbacks),
      c_problem(c_problem),
      raw_data(raw_data),
      interpolated_data(std::make_unique<Trajectory[]>(c_problem->data_file_count))
{
    fill_runtime_parameters();
    fill_data(mesh);
}

Problem Problem::create(c_problem_t* c_problem, const Mesh& mesh) {
    return Problem(c_problem, mesh, create_raw_data(c_problem));
}


/* here is probably not the correct place for these */
/*
void __used_in_nominal_scaling_factor(c_problem_t* c_problem) {
    auto x_nominal = assign_or_one(c_problem->x_nominal, c_problem->x_size);
    auto u_nominal = assign_or_one(c_problem->u_nominal, c_problem->u_size);
    auto p_nominal = assign_or_one(c_problem->p_nominal, c_problem->p_size);

    auto obj_nominal = (c_problem->obj_nominal ? *c_problem->obj_nominal : 1.0);
    auto f_nominal = assign_or_one(c_problem->f_nominal, c_problem->x_size);
    auto g_nominal = assign_or_one(c_problem->g_nominal, c_problem->g_size);
    auto r_nominal = assign_or_one(c_problem->r_nominal, c_problem->r_size);
}
*/

} // namespace C
