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

#ifndef MOO_INTEGRATOR_H
#define MOO_INTEGRATOR_H

#include <base/util.h>
#include <base/export.h>
#include <base/fixed_vector.h>
#include <base/trajectory.h>

#include <simulation/integrator/integrator_util.h>

namespace Simulation {

class MOO_EXPORT Integrator {
public:
    Integrator(ODEFunction ode_fn,
               std::vector<f64> dense_output_grid,
               f64* x_start_values,
               int x_size,
               void* user_data,
               const f64* parameters,
               int p_size,
               const ControlTrajectory* controls,
               JacobianFunction jac_fn,
               Jacobian jac_pattern);

    virtual ~Integrator() = default;
    Integrator(Integrator&&) noexcept = default;

    std::unique_ptr<Trajectory> simulate();
    std::unique_ptr<Trajectory> simulate(f64* x_start_values_, std::vector<f64> dense_output_grid_);
    std::unique_ptr<Trajectory> simulate(f64* x_start_values_, f64 t0_, f64 tf_, int steps_);

    void get_ode(f64 t, f64* x, f64* out);
    void get_dense_jacobian(f64 t, f64* x, f64* out);
    // f64* get_sparse_jacobian(f64 t, f64* x);

    ODEFunction ode_func;
    JacobianFunction jac_func;

    std::vector<f64> dense_output_grid;

    std::unique_ptr<Trajectory> output;

    f64* x_start_values;

    void set_controls(f64 t);

    FixedVector<f64> u;

    int x_size;
    int u_size;
    int p_size;

private:
    virtual int internal_simulate() = 0;

    void* user_data;

    const ControlTrajectory* internal_controls;
    const f64* parameters;

    Jacobian jac_pattern;

    f64 last_t;

    FixedVector<f64> sparse_jac;
};

} // namespace Simulation

#endif // MOO_INTEGRATOR_H
