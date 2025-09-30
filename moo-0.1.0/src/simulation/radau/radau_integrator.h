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

#ifndef MOO_RADAU_WRAPPER_H
#define MOO_RADAU_WRAPPER_H

#include <base/export.h>
#include <base/log.h>
#include <simulation/integrator/integrator.h>

extern "C" {
    void radau_(
        int* n,
        void (*fcn)(int*, f64*, f64*, f64*, f64*, int*),
        f64* x,
        f64* y,
        f64* xend,
        f64* h,
        f64* rtol,
        f64* atol,
        int* itol,
        void (*jac)(int*, f64*, f64*, f64*, int*, f64*, int*),
        int* ijac,
        int* mljac,
        int* mujac,
        void (*mas)(int*, f64*, int*, f64*, int*),
        int* imas,
        int* mlmas,
        int* mumas,
        void (*solout)(int*, f64*, f64*, f64*, f64*, int*, int*, f64*, int*, int*),
        int* iout,
        f64* work,
        int* lwork,
        int* iwork,
        int* liwork,
        f64* rpar,
        int* ipar,
        int* idid
    );

    f64 contra_(
        int* i,
        f64* s,
        f64* cont,
        int* lrc
    );
}

namespace Simulation {

enum RadauScheme {
    ADAPTIVE = 0,
    ONE = 1,
    FIVE = 5,
    NINE = 9,
    THIRTEEN = 13
};

class MOO_EXPORT RadauIntegrator : public Integrator {

friend class RadauBuilder;

public:
    RadauIntegrator(RadauIntegrator&&) noexcept = default;

    int internal_simulate() override;

    RadauScheme scheme;
    f64 h_init;
    f64 atol;
    f64 rtol;
    int max_it;

    int lwork;
    int liwork;

    std::vector<f64> work;
    std::vector<int> iwork;

    int ijac;

    /* args we dont change at all */
    int itol = 0;
    int mljac = 0;
    int mujac = 0;
    int imas = 0;
    int mlmas = 0;
    int mumas = 0;
    int return_code = 0;
    size_t dense_output_index = 0;
    int iout = 1;
    f64 rpar = 0;
    int ipar = 0;

private:
    RadauIntegrator(ODEFunction ode_fn,
                    std::vector<f64> dense_output_grid,
                    f64* x_start_values,
                    int x_size,
                    void* user_data,
                    const f64* parameters,
                    int p_size ,
                    const ControlTrajectory* controls,
                    JacobianFunction jac_fn,
                    Jacobian jac_pattern,
                    RadauScheme scheme,
                    f64 h_init,
                    f64 atol,
                    f64 rtol,
                    int max_it);
};

} // namespace Simulation

#endif // MOO_RADAU_WRAPPER_H
