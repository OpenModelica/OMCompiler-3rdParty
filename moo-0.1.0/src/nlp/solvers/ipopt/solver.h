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

#ifndef MOO_IPOPT_SOLVER_H
#define MOO_IPOPT_SOLVER_H

#include <memory>

#include <nlp/nlp_solver.h>
#include <base/export.h>
#include <base/timing.h>

namespace IpoptSolver {

struct IpoptSolverData;
class IpoptTimingNode;

class MOO_EXPORT IpoptSolver : public NLP::NLPSolver {
    friend IpoptTimingNode;

public:
    IpoptSolver(NLP::NLP& nlp, NLP::NLPSolverSettings& solver_settings);

    virtual ~IpoptSolver();

    void optimize() override;
    void init_application();
    void set_settings();

    int get_iterations() const override;
    f64 get_total_time() const override;
    f64 get_solver_time() const override;
    f64 get_callback_time() const override;

private:
   IpoptSolverData* ipdata;
};

class IpoptTimingNode : public TimingNode {
    IpoptSolver* ipopt_solver;
    int total_iterations = 0;

    // TODO: maybe more here in the future?

public:
    IpoptTimingNode(std::string n, TimingNode* p = nullptr, IpoptSolver* ipopt_solver = nullptr);

    void finalize() override;
};

} // namespace IpoptSolver

#endif // MOO_IPOPT_SOLVER_H
