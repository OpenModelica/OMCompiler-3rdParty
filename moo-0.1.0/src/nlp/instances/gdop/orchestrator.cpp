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
#include <base/timing.h>

namespace GDOP {

Orchestrator::Orchestrator(GDOP& gdop,
                           std::unique_ptr<Strategies> strategies,
                           NLP::NLPSolver& solver)
    : gdop(gdop),
      strategies(std::move(strategies)),
      solver(solver) {}

MeshRefinementOrchestrator::MeshRefinementOrchestrator(GDOP& gdop,
                                                       std::unique_ptr<Strategies> strategies,
                                                       NLP::NLPSolver& solver)
    : Orchestrator(gdop, std::move(strategies), solver) {}

void MeshRefinementOrchestrator::optimize() {
    { ScopedTimer optimize{"MeshRefinementOrchestrator::optimize"};

    strategies->reset(gdop);

    auto initial_guess = strategies->get_initial_guess(gdop);

    gdop.set_scaling_factory(strategies->scaling_factory);

    for(;;) {
        gdop.set_initial_guess(std::move(initial_guess));

        solver.optimize();

        // append metadata to mesh refinement history
        history.add_block(gdop, solver);

        // === mesh refinement ===

        // 1. detect intervals and degrees (new vectors)
        auto mesh_update = strategies->detect(gdop.get_mesh(), *gdop.get_optimal_solution());

        if (!mesh_update) { break; }

        // 2. create refined Mesh
        auto refined_mesh = gdop.get_mesh().create_from_mesh_update(std::move(mesh_update));

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
    strategies->emit(*gdop.get_optimal_solution()); }

    history.finalize().print();

    //gdop.get_optimal_solution()->costates->to_csv("costates_final.csv", false);
    //gdop.get_optimal_solution()->lower_costates->to_csv("lower_costates_final.csv", false);
    //gdop.get_optimal_solution()->upper_costates->to_csv("upper_costates_final.csv", false);
}

MeshRefinementHistoryBlock::MeshRefinementHistoryBlock(const GDOP& gdop, const NLP::NLPSolver& solver)
    : objective(gdop.get_objective_value()),
      intervals(gdop.get_mesh().intervals),
      nlp_solver_iters(solver.get_iterations()),
      nlp_solver_total_nano(solver.get_total_time()),
      nlp_solver_self_nano(solver.get_solver_time()),
      nlp_solver_callback_nano(solver.get_callback_time()) {}

void MeshRefinementHistory::add_block(const GDOP& gdop, const NLP::NLPSolver& solver) {
    blocks.push_back(MeshRefinementHistoryBlock(gdop, solver));
}

MeshRefinementHistory& MeshRefinementHistory::finalize() {
    strategy_timings_nano = Timing::accumulate_blocks("Strategies::", "IpoptSolver::optimize");
    return *this;
}

void MeshRefinementHistory::print() {
    FixedTableFormat<8> ftf(
        {9, 20, 9, 10, 14, 15, 18, 15},
        {Align::Right, Align::Right, Align::Right, Align::Right,
         Align::Right, Align::Right, Align::Right, Align::Right});

    Log::start_module(ftf, "Mesh Refinement History");
    Log::row(ftf, "Iteration", "Objective", "Intervals", "Iterations",
                  "NLP Total [ms]", "NLP Solver [ms]", "NLP Callbacks [ms]", "Strategies [ms]");
    Log::dashes(ftf);

    f64 total_solver_total_nano = 0.0;
    f64 total_solver_self_nano = 0.0;
    f64 total_solver_callback_nano = 0.0;
    f64 total_strategy_nano = 0.0;
    int total_iters = 0;

    for (int i = 0; i < int_size(blocks); i++) {
        const auto& e = blocks[i];
        Log::row(ftf,
            fmt::format("{:>9}", i),
            fmt::format("{:>20.10e}", e.objective),
            fmt::format("{:>9}", e.intervals),
            fmt::format("{:>10}", e.nlp_solver_iters),
            fmt::format("{:>14.3f}", Timing::nano_to_ms(e.nlp_solver_total_nano)),
            fmt::format("{:>15.3f}", Timing::nano_to_ms(e.nlp_solver_self_nano)),
            fmt::format("{:>18.3f}", Timing::nano_to_ms(e.nlp_solver_callback_nano)),
            fmt::format("{:>15.3f}", Timing::nano_to_ms(strategy_timings_nano[i])));

        total_iters += e.nlp_solver_iters;
        total_solver_total_nano += e.nlp_solver_total_nano;
        total_solver_self_nano += e.nlp_solver_self_nano;
        total_solver_callback_nano += e.nlp_solver_callback_nano;
        total_strategy_nano += strategy_timings_nano[i];
    }

    Log::dashes(ftf);

    Log::row(ftf,
    fmt::format("{:>9}", "Overall"),
    fmt::format("{:>20.10e}", blocks.back().objective),
    fmt::format("{:>9}", blocks.back().intervals),
    fmt::format("{:>10}", total_iters),
    fmt::format("{:>14.3f}", Timing::nano_to_ms(total_solver_total_nano)),
    fmt::format("{:>15.3f}", Timing::nano_to_ms(total_solver_self_nano)),
    fmt::format("{:>18.3f}", Timing::nano_to_ms(total_solver_callback_nano)),
    fmt::format("{:>15.3f}", Timing::nano_to_ms(total_strategy_nano)));

    Log::dashes_ln(ftf);
}

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

} // namespace GDOP
