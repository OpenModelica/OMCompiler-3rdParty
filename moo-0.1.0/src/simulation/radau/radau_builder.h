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

#ifndef MOO_RADAU_BUILDER_H
#define MOO_RADAU_BUILDER_H

#include <vector>

#include <base/export.h>
#include <simulation/integrator/builder.h>
#include <simulation/radau/radau_integrator.h>

namespace Simulation {

class MOO_EXPORT RadauBuilder : public IntegratorBuilder<RadauBuilder, RadauIntegrator>{
public:
    RadauBuilder() : IntegratorBuilder() {}

    RadauBuilder& radau_scheme(RadauScheme radau_scheme_);
    RadauBuilder& radau_h0(f64 h_init_);
    RadauBuilder& radau_tol(f64 atol_, f64 rtol_);
    RadauBuilder& radau_max_it(int max_it_);

    RadauIntegrator build() const override;

private:
    RadauScheme scheme = RadauScheme::ADAPTIVE;
    f64 h_init = 1e-6;
    f64 atol = 1e-10;
    f64 rtol = 1e-10;
    int max_it = 100000;
};

} // namespace Simulation

#endif // MOO_RADAU_BUILDER_H
