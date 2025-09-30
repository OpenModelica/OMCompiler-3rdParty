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

#include <iomanip>
#include <cmath>

#include <simulation/radau/radau_builder.h>
#include <simulation/radau/test.h>
#include <base/log.h>

namespace Simulation {

void fcn(const f64* x, const f64* u, const f64* p, f64 t, f64* dxdt, void* user_data) {
    dxdt[0] = -x[0];
    dxdt[1] = x[0] - 2 * x[1] * p[1] + p[0] + u[0];
}

void jac(const f64* x, const f64* u, const f64* p, f64 t, f64* J, void* user_data) {
    J[0] = -1.0;
    J[1] = 1.0;
    J[2] = -2 * p[1];
}

int radau_wrapper_test() {
    FixedVector<f64> x_start = {1, -1};

    auto dummy_control = ControlTrajectory{{0, 1}, {{-1, 1}}};

    int nnz = 3;
    int rowidx[] = { 0, 1, 1 };
    int colptr[] = { 0, 2, 3 };

    Jacobian jac_pattern = Jacobian::sparse(JacobianFormat::CSC, rowidx, colptr, nnz);

    f64 parameters[] = { -1.0, 1.0 };

    auto radau_integrator = RadauBuilder()
                                .ode(fcn)
                                .states(2, x_start.raw())
                                .control(&dummy_control)
                                .params(2, parameters)
                                .jacobian(jac, jac_pattern)
                                .interval(0, 1, 10)
                                .radau_scheme(RadauScheme::ADAPTIVE)
                                .radau_h0(1e-5)
                                .radau_tol(1e-12, 1e-12)
                                .build();

    auto out = radau_integrator.simulate();

    out->print_table();

    return radau_integrator.return_code;
}

} // namespace Simulation
