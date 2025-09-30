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

#include "orchestrator.h"

namespace GDOP {

void MeshRefinementOrchestrator::optimize() {
    strategies->reset(gdop);

    auto initial_guess = strategies->get_initial_guess(gdop);

    gdop.set_scaling_factory(strategies->scaling_factory);

    for(;;) {
        gdop.set_initial_guess(std::move(initial_guess));

        solver.optimize();

        // === mesh refinement ===

        // 1. detect intervals and degrees (new vectors)
        auto mesh_update = strategies->detect(gdop.get_mesh(), *gdop.get_optimal_solution());

        if (!mesh_update) { break; }

        // 2. create refined Mesh
        auto refined_mesh = Mesh::create_from_mesh_update(std::move(mesh_update));

        // 3. interpolate (x*, lambda*, z*) to new mesh -> new initial guess
        initial_guess = strategies->get_refined_initial_guess(gdop.get_mesh(), *refined_mesh, *gdop.get_optimal_solution());
        solver.solver_settings.set(NLP::Option::WarmStart, true);

        //initial_guess->costates->to_csv("costates_interp.csv", false);
        //initial_guess->lower_costates->to_csv("lower_costates_interp.csv", false);
        //initial_guess->upper_costates->to_csv("upper_costates_interp.csv", false);

        // 4. update gdop with new mesh
        gdop.update(refined_mesh);
    }

    strategies->verify(gdop, *gdop.get_optimal_solution());
    strategies->emit(*gdop.get_optimal_solution());

    //gdop.get_optimal_solution()->costates->to_csv("costates_final.csv", false);
    //gdop.get_optimal_solution()->lower_costates->to_csv("lower_costates_final.csv", false);
    //gdop.get_optimal_solution()->upper_costates->to_csv("upper_costates_final.csv", false);
}

} // namespace GDOP


/*
TODO:
    const auto& optimum = *gdop.get_optimal_solution();
    const auto& controls = optimum.primals->copy_extract_controls();

    strategies->simulate_step_activate(controls, FixedVector<f64>(optimum.primals->p));
    auto res = strategies->simulate_step(optimum.primals->extract_initial_states().raw(), optimum.primals->t[0], optimum.primals->t.back() / 2);
    res->print_table();
    auto newx = res->extract_final_states();
    res = strategies->simulate_step(newx.raw(), optimum.primals->t.back() / 2 , optimum.primals->t.back());
    res->print_table();

    strategies->simulate_step_reset();
*/