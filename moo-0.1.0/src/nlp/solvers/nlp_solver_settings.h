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

#ifndef MOO_NLP_SOLVER_SETTINGS_H
#define MOO_NLP_SOLVER_SETTINGS_H

#include <string>
#include <unordered_map>
#include <variant>
#include <optional>

#include <base/util.h>
#include <base/export.h>

namespace NLP {

enum class HessianOption {
    Exact,
    LBFGS,
    Const
};

enum class LinearSolverOption {
    MUMPS,
    MA27,
    MA57,
    MA77,
    MA86,
    MA97
};

enum class NLPSolverOption {
    Ipopt
};

enum class Option {
/* HessianOption         */    Hessian,
/* f64                   */    Tolerance,
/* int                   */    Iterations,
/* f64                   */    CPUTime,
/* LinearSolverOption    */    LinearSolver,
/* NLPSolverOption       */    NLPSolver,
/* bool                  */    IpoptDerivativeTest,
/* bool                  */    WarmStart,
/* bool                  */    QP,
};

using OptionValue = std::variant<std::string, f64, int, bool, HessianOption, LinearSolverOption, NLPSolverOption>;

class MOO_EXPORT NLPSolverSettings {
public:
    NLPSolverSettings(int argc, char** argv);

    void print() const;

    void set(Option option, const OptionValue& value);
    const OptionValue& get(Option option) const;

    bool option_is_true(Option option) const;
    bool option_matches(Option option, const std::string& str) const;

    template<typename T>
    T get_or_default(Option option) const;

private:
    std::unordered_map<Option, OptionValue> settings;
};

// Option enum to string
std::string to_string(Option option);

// string to Option enum
std::optional<Option> option_from_string(const std::string& name);

extern const std::unordered_map<Option, OptionValue> default_settings;

} // namespace NLP

#endif // MOO_NLP_SOLVER_SETTINGS_H
