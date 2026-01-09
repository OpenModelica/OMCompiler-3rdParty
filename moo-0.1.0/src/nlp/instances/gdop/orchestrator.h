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

#ifndef MOO_GDOP_ORCHESTRATOR_H
#define MOO_GDOP_ORCHESTRATOR_H

#include <base/export.h>
#include <nlp/nlp_solver.h>
#include "gdop.h"

namespace GDOP {

class MOO_EXPORT Orchestrator {
public:
    GDOP& gdop;
    std::unique_ptr<Strategies> strategies;
    NLP::NLPSolver& solver;

    Orchestrator(GDOP& gdop,
                 std::unique_ptr<Strategies> strategies,
                 NLP::NLPSolver& solver);

    virtual ~Orchestrator() = default;

    virtual void optimize() = 0;
};

class MOO_EXPORT MeshRefinementHistoryBlock {
    friend class MeshRefinementHistory;

    f64 objective;
    int intervals;
    int nlp_solver_iters;
    f64 nlp_solver_total_nano;
    f64 nlp_solver_self_nano;
    f64 nlp_solver_callback_nano;

    MeshRefinementHistoryBlock(const GDOP& gdop, const NLP::NLPSolver& solver);
};

class MOO_EXPORT MeshRefinementHistory {
public:
    std::vector<MeshRefinementHistoryBlock> blocks;
    std::vector<f64> strategy_timings_nano;

    MeshRefinementHistory() = default;

    void add_block(const GDOP& gdop, const NLP::NLPSolver& solver);
    MeshRefinementHistory& finalize();
    void print();
};

class MOO_EXPORT MeshRefinementOrchestrator : public Orchestrator {
public:
    MeshRefinementOrchestrator(GDOP& gdop,
                               std::unique_ptr<Strategies> strategies,
                               NLP::NLPSolver& solver);

    void optimize() override;

private:
    MeshRefinementHistory history{};
};

} // namespace GDOP

#endif // MOO_GDOP_ORCHESTRATOR_H
