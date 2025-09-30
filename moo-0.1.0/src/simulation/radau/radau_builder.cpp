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

#include <simulation/radau/radau_builder.h>

namespace Simulation {

RadauBuilder& RadauBuilder::radau_scheme(RadauScheme radau_scheme_) {
    scheme = radau_scheme_;
    return *this;
}

RadauBuilder& RadauBuilder::radau_h0(f64 h_init_) {
    h_init = h_init_;
    return *this;
}

RadauBuilder& RadauBuilder::radau_tol(f64 atol_, f64 rtol_) {
    atol = atol_;
    rtol = rtol_;
    return *this;
}

RadauBuilder& RadauBuilder::radau_max_it(int max_it_) {
    max_it = max_it_;
    return *this;
}

RadauIntegrator RadauBuilder::build() const {
    f64 h0 = h_init != 0.0 ? h_init : (dense_output_grid.back() - dense_output_grid[0]) / (2 * dense_output_grid.size());

    return RadauIntegrator(
        /* every Integrator */
        ode_func, dense_output_grid,
        x_start_values, x_size,
        user_data, parameters, p_size,
        controls,
        jac_func, jac_pattern,
        /* Radau specific */
        scheme,
        h0,
        atol,
        rtol,
        max_it
    );
}

} // namespace Simulation
