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

#ifndef MOO_C_GDOPT_MAIN
#define MOO_C_GDOPT_MAIN

#ifdef __cplusplus
extern "C" {
#endif

#include <base/export.h>
#include <interfaces/c/structures.h>

MOO_EXPORT int main_gdopt(int argc, char** argv, c_problem_t* c_problem);

#ifdef __cplusplus
}
#endif

#endif // MOO_C_GDOPT_MAIN