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

#ifndef MOO_SIMULATION_BUILDER_H
#define MOO_SIMULATION_BUILDER_H

#include <vector>

#include <base/export.h>
#include <simulation/integrator/integrator.h>

namespace Simulation {

template <typename BuilderImpl, typename IntegratorImpl>
class MOO_EXPORT IntegratorBuilder {
public:
    virtual ~IntegratorBuilder() = default;

    BuilderImpl& states(int x_size_, f64* x_start_values_ = nullptr) {
        x_size = x_size_;
        x_start_values = x_start_values_;
        return static_cast<BuilderImpl&>(*this);
    }

    BuilderImpl& ode(ODEFunction ode_func_) {
        ode_func = ode_func_;
        return static_cast<BuilderImpl&>(*this);
    }

    BuilderImpl& interval(f64 t0_, f64 tf_, int steps_) {
        if (steps_ <= 0 || tf_ <= t0_) Log::error("Invalid interval");

        t0 = t0_;
        tf = tf_;
        num_steps = steps_;

        dense_output_grid = std::vector<f64>(num_steps + 1);
        f64 dt = (tf - t0) / num_steps;
        for (int i = 0; i <= num_steps; i++) dense_output_grid[i] = t0 + i * dt;

        return static_cast<BuilderImpl&>(*this);
    }

    BuilderImpl& grid(const std::vector<f64> dense_output_grid_) {
        dense_output_grid = dense_output_grid_;
        return static_cast<BuilderImpl&>(*this);
    }

    BuilderImpl& userdata(void* user_data_) {
        user_data = user_data_;
        return static_cast<BuilderImpl&>(*this);
    }

    BuilderImpl& params(int p_size_, const f64* parameters_ = nullptr) {
        p_size = p_size_;
        parameters = parameters_;
        return static_cast<BuilderImpl&>(*this);
    }

    BuilderImpl& control(const ControlTrajectory* controls_) {
        controls = controls_;
        return static_cast<BuilderImpl&>(*this);
    }

    BuilderImpl& jacobian(JacobianFunction jac_func_,
                          Jacobian jac_pattern_) {
        jac_func = jac_func_;
        jac_pattern = jac_pattern_;
        return static_cast<BuilderImpl&>(*this);
    }

    virtual IntegratorImpl build() const = 0;

protected:
    ODEFunction ode_func = nullptr;
    f64* x_start_values = nullptr;
    int x_size = 0;

    std::vector<f64> dense_output_grid{};

    f64 t0 = 0.0;
    f64 tf = 0.0;
    int num_steps = 1;

    void* user_data = nullptr;
    const f64* parameters = nullptr;
    int p_size = 0;

    const ControlTrajectory* controls = nullptr;

    JacobianFunction jac_func = nullptr;
    Jacobian jac_pattern = Jacobian::dense();
};

} // namespace Simulation

#endif // MOO_SIMULATION_BUILDER_H
