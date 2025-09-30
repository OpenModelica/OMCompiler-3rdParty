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

#include <simulation/integrator/integrator.h>

namespace Simulation {

Integrator::Integrator(ODEFunction ode_fn,
                       std::vector<f64> grid,
                       f64* x_start_values,
                       int x_size,
                       void* user_data,
                       const f64* parameters,
                       int p_size,
                       const ControlTrajectory* controls,
                       JacobianFunction jac_fn,
                       Jacobian jac_pattern)
    : ode_func(ode_fn),
      jac_func(jac_fn),
      dense_output_grid(grid),
      output(nullptr),
      x_start_values(x_start_values),
      u(FixedVector<f64>(controls ? controls->u.size() : 0)),
      x_size(x_size),
      u_size(controls ? controls->u.size() : 0),
      p_size(p_size),
      user_data(user_data),
      internal_controls(controls),
      parameters(parameters),
      jac_pattern(jac_pattern),
      last_t(MINUS_INFINITY),
      sparse_jac(FixedVector<f64>(this->jac_pattern.nnz)) {}

/**
 * @brief Override current dense_output_grid and start simulating at new states x_start_values
 */
std::unique_ptr<Trajectory> Integrator::simulate(f64* x_start_values_, std::vector<f64> dense_output_grid_) {
    x_start_values = x_start_values_;
    dense_output_grid = dense_output_grid_;
    return simulate();
}

/**
 * @brief Override current dense_output_grid and start simulating at new states x_start_values
 */
std::unique_ptr<Trajectory> Integrator::simulate(f64* x_start_values_, f64 t0_, f64 tf_, int steps_) {
    x_start_values = x_start_values_;
    dense_output_grid = std::vector<f64>(steps_ + 1);
    f64 dt = (tf_ - t0_) / steps_;
    for (int i = 0; i <= steps_; i++) dense_output_grid[i] = t0_ + i * dt;
    return simulate();
}

std::unique_ptr<Trajectory> Integrator::simulate() {
    if (internal_controls) set_controls(dense_output_grid[0]);

    output = std::make_unique<Trajectory>();
    output->t.reserve(dense_output_grid.size());
    output->x.resize(x_size);
    output->u.resize(u_size);

    for (int x_idx = 0; x_idx < x_size; x_idx++) { output->x[x_idx].reserve(dense_output_grid.size()); }
    for (int u_idx = 0; u_idx < u_size; u_idx++) { output->u[u_idx].reserve(dense_output_grid.size()); }

    int err = internal_simulate();

    if (err < 0) {
        return nullptr;
    }

    output->p.resize(p_size);

    for (int p_idx = 0;  p_idx < p_size; p_idx++) {
        output->p[p_idx] = parameters[p_idx];
    }

    return std::move(output);
}

void Integrator::set_controls(f64 t) {
    if (internal_controls && t != last_t) {
        internal_controls->interpolate_at(t, u.raw());
    }

    last_t = t;
}

void Integrator::get_dense_jacobian(f64 t, f64* x, f64* out) {
    if (!jac_func) return;

    set_controls(t);

    if (jac_pattern.jfmt == JacobianFormat::DENSE) {
        jac_func(x, u.raw(), parameters, t, out, user_data);
    }
    else if (jac_pattern.jfmt == JacobianFormat::COO) {
        sparse_jac.fill_zero();
        jac_func(x, u.raw(), parameters, t, sparse_jac.raw(), user_data);

        for (int nz = 0; nz < jac_pattern.nnz; nz++) out[jac_pattern.i_row[nz] * x_size + jac_pattern.j_col[nz]] = sparse_jac[nz];

    }
    else if (jac_pattern.jfmt == JacobianFormat::CSC) {
        sparse_jac.fill_zero();
        jac_func(x, u.raw(), parameters, t, sparse_jac.raw(), user_data);

        for (int col = 0; col < x_size; col++) {
            for (int nz = jac_pattern.j_col[col]; nz < jac_pattern.j_col[col + 1]; nz++) {
                out[jac_pattern.i_row[nz] * x_size + col] = sparse_jac[nz];
            }
        }
    }
}

void Integrator::get_ode(f64 t, f64* x, f64* out) {
    set_controls(t);
    ode_func(x, u.raw(), parameters, t, out, user_data);
}

} // namespace Simulation
