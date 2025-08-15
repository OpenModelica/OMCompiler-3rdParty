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

#include "wrapper.h"

// TODO: make this more simple: just provide tolerances time ControlTrajectory etc. => Obtain simulated Trajectory (optionally at predefined points)

void radau_solver(
    int* n,
    void (*fcn)(int*, double*, double*, double*),
    double* x,
    double* y,
    double* xend,
    double* h,
    double* rtol,
    double* atol,
    int* itol,
    void (*jac)(int*, double*, double*, int*, int*, double*, double*),
    int* ijac,
    int* mljac,
    int* mujac,
    void (*mas)(int*, int*, double*),
    int* imas,
    int* mlmas,
    int* mumas,
    void (*solout)(int*, double*, double*, double*, double*),
    int* iout,
    double* work,
    int* lwork,
    int* iwork,
    int* liwork,
    double* rpar,
    int* ipar,
    int* idid,
    Integrator integrator)
{
    switch (integrator)
    {
        case RADAU_5:
            iwork[10] = 3; /* min_m */
            iwork[11] = 3; /* max_m */
            iwork[12] = 3; /* start_m */
            break;

        case RADAU_9:
            iwork[10] = 5; /* min_m */
            iwork[11] = 5; /* max_m */
            iwork[12] = 5; /* start_m */
            break;

        case RADAU_13:
            iwork[10] = 7; /* min_m */
            iwork[11] = 7; /* max_m */
            iwork[12] = 7; /* start_m */
            break;

        case RADAU_ADAPTIVE:
        default:
            iwork[10] = 1; /* min_m */
            iwork[11] = 7; /* max_m */
            iwork[12] = 5; /* start_m */
            break;
    }

    radau_(
        n, fcn, x, y, xend, h, rtol, atol, itol,
        jac, ijac, mljac, mumas, mas, imas, mlmas, mumas,
        solout, iout, work, lwork, iwork, liwork, rpar, ipar, idid
    );
}
