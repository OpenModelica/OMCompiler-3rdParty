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

enum Integrator {
    RADAU_ADAPTIVE = 0,
    RADAU_5 = 5,
    RADAU_9 = 9,
    RADAU_13 = 13
};

extern "C" {
    void radau_(
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
        int* idid
    );
}

MOO_EXPORT void radau_solver(
    int* n,           // N: size of the problem
    void (*fcn)(int*, double*, double*, double*), // FCN: function for dy/dx
    double* x,        // X: initial value of the independent variable
    double* y,        // Y: solution array
    double* xend,     // XEND: final value of the independent variable
    double* h,        // H: initial step size
    double* rtol,     // RTOL: relative tolerance
    double* atol,     // ATOL: absolute tolerance
    int* itol,        // ITOL: tolerance type
    void (*jac)(int*, double*, double*, int*, int*, double*, double*), // JAC: function for the Jacobian
    int* ijac,        // IJAC: jacobian type
    int* mljac,       // MLJAC: lower bandwidth of the Jacobian
    int* mujac,       // MUJAC: upper bandwidth of the Jacobian
    void (*mas)(int*, int*, double*), // MAS: function for the mass matrix
    int* imas,        // IMAS: mass matrix type
    int* mlmas,       // MLMAS: lower bandwidth of the mass matrix
    int* mumas,       // MUMAS: upper bandwidth of the mass matrix
    void (*solout)(int*, double*, double*, double*, double*), // SOLOUT: output function
    int* iout,        // IOUT: output type
    double* work,     // WORK: double precision work array
    int* lwork,       // LWORK: size of WORK
    int* iwork,       // IWORK: integer work array
    int* liwork,      // LIWORK: size of IWORK
    double* rpar,     // RPAR: user-defined double precision parameters
    int* ipar,        // IPAR: user-defined integer parameters
    int* idid,        // IDID: success indicator
    Integrator integrator
);

#endif // MOO_RADAU_WRAPPER_H
