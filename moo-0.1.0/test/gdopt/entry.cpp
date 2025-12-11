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

#include <simulation/radau/test.h>
#include <generated.h>
#include <sanity_check.h>
#include <full_functionality_check.h>

#include <base/log.h>

int main(int argc, char** argv) {
    // simulation test
    Simulation::radau_wrapper_test();

    // test for future generated code
    main_generated(argc, argv);

    // sanity check, testing time optimization for a minimal DAE example
    main_sanity_check(argc, argv);

    // full functionality test, including all types of constraints + variables + derivative test
    main_full_functionality(argc, argv);
}
