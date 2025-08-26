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
                 NLP::NLPSolver& solver)
    : gdop(gdop),
      strategies(std::move(strategies)),
      solver(solver) {}

    virtual ~Orchestrator() = default;

    virtual void optimize() = 0;
};

class MOO_EXPORT MeshRefinementOrchestrator : public Orchestrator {
public:
    MeshRefinementOrchestrator(GDOP& gdop,
                               std::unique_ptr<Strategies> strategies,
                               NLP::NLPSolver& solver)
    : Orchestrator(gdop, std::move(strategies), solver) {}

    void optimize() override;
};

} // namespace GDOP

#endif // MOO_GDOP_ORCHESTRATOR_H
