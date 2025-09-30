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

#ifndef MOO_C_PROBLEM_H
#define MOO_C_PROBLEM_H

#include <nlp/instances/gdop/problem.h>
#include <nlp/instances/gdop/gdop.h>

#include <interfaces/c/structures.h>

namespace C {

class FullSweep : public GDOP::FullSweep {
public:
    FullSweep(GDOP::FullSweepLayout&& layout_in,
              const GDOP::ProblemConstants& pc,
              c_problem_t* c_problem)
        : GDOP::FullSweep(std::move(layout_in), pc),
          c_callbacks(c_problem->callbacks),
          c_problem(c_problem) {};

    void callback_eval(const f64* xu_nlp, const f64* p) override;
    void callback_jac(const f64* xu_nlp, const f64* p) override;
    void callback_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, const f64* lambda) override;

private:
    inline f64* get_data_ij(int interval_i, int node_j) {
        if (!c_problem->data) { return nullptr; }
        int node_idx = 1 + pc.mesh->acc_nodes[interval_i][node_j];
        return c_problem->data + node_idx * c_problem->data_chunk_size;
    }

    c_callbacks_t* c_callbacks;
    c_problem_t* c_problem;
};

class BoundarySweep : public GDOP::BoundarySweep {
public:
    BoundarySweep(GDOP::BoundarySweepLayout&& layout_in,
                  const GDOP::ProblemConstants& pc,
                  c_problem_t* c_problem)
        : GDOP::BoundarySweep(std::move(layout_in), pc),
          c_callbacks(c_problem->callbacks),
          c_problem(c_problem) {};

    void callback_eval(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) override;
    void callback_jac(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) override;
    void callback_hes(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, const f64 mayer_factor, const f64* lambda) override;

private:
    inline f64* get_data_t0() {
        if (!c_problem->data) { return nullptr; }
        return c_problem->data;
    }

    inline f64* get_data_tf() {
        if (!c_problem->data) { return nullptr; }
        return c_problem->data + pc.mesh->node_count * c_problem->data_chunk_size;
    }

    c_callbacks_t* c_callbacks;
    c_problem_t* c_problem;
};

class Dynamics : public GDOP::Dynamics {
public:
    Dynamics(const GDOP::ProblemConstants& pc,
             c_problem_t* c_problem_,
             std::shared_ptr<Trajectory[]> raw_data);

    f64* get_data(f64 t);
    void eval(const f64* x, const f64* u, const f64* p, f64 t, f64* f, void* user_data) override;
    void jac(const f64* x, const f64* u, const f64* p, f64 t, f64* dfdx, void* user_data)  override;

private:
    c_callbacks_t* c_callbacks;
    c_problem_t* c_problem;
    std::unique_ptr<ControlTrajectory[]> raw_ctrl_data;
    FixedVector<f64> current_data;
    mutable f64 current_data_time;
};

class Problem : public GDOP::Problem {
public:
    static Problem create(c_problem_t* c_problem, const Mesh& mesh);

    c_callbacks_t* c_callbacks;
    c_problem_t* c_problem;

private:
    Problem(c_problem_t* c_problem, const Mesh& mesh, std::shared_ptr<Trajectory[]> raw_data);

    std::shared_ptr<Trajectory[]> raw_data;           // raw array of Trajectories, read from the user
    std::unique_ptr<Trajectory[]> interpolated_data;  // raw_data but interpolated onto a given mesh
    FixedVector<f64> flat_interpolated_input_data;    // flat representation of interpolated_data.u, owner field of c_problem->data

    void fill_runtime_parameters();
    void fill_data(const Mesh& mesh);
};

// === helpers ===

size_t accumulate_data_size(c_problem_t* c_problem, std::unique_ptr<ControlTrajectory[]>& raw_ctrl_data);

std::unique_ptr<ControlTrajectory[]> extract_ctrl_data(c_problem_t* c_problem, std::shared_ptr<Trajectory[]>& raw_data);

} // namespace C

#endif // MOO_C_PROBLEM_H
