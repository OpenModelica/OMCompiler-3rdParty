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
#include <sstream>
#include <iomanip>
#include <cstdlib>

#include <base/log.h>

#include "nlp_solver_settings.h"

// TODO: refactor this into big struct or what is best if the size expands?

namespace NLP {

// === Default values ===

const std::unordered_map<Option, OptionValue> default_settings = {
    {Option::Hessian,             HessianOption::Exact},
    {Option::Tolerance,           1e-10},
    {Option::Iterations,          5000},
    {Option::CPUTime,             3600.0},
    {Option::LinearSolver,        LinearSolverOption::MUMPS},
    {Option::NLPSolver,           NLPSolverOption::Ipopt},
    {Option::IpoptDerivativeTest, false},
    {Option::WarmStart,           false},
    {Option::QP,                  false},
};

// === Option string conversions ===

std::string to_string(Option option) {
    switch (option) {
        case Option::Hessian:             return "Hessian";
        case Option::Tolerance:           return "Tolerance";
        case Option::Iterations:          return "Iterations";
        case Option::CPUTime:             return "CPUTime";
        case Option::LinearSolver:        return "LinearSolver";
        case Option::NLPSolver:           return "NLPSolver";
        case Option::IpoptDerivativeTest: return "IpoptDerivativeTest";
        case Option::WarmStart:           return "WarmStart";
        case Option::QP:                  return "QP";
        default:                          return "Unknown";
    }
}

std::optional<Option> option_from_string(const std::string& name) {
    static const std::unordered_map<std::string, Option> map = {
        {"Hessian", Option::Hessian},
        {"Tolerance", Option::Tolerance},
        {"Iterations", Option::Iterations},
        {"CPUTime", Option::CPUTime},
        {"LinearSolver", Option::LinearSolver},
        {"NLPSolver", Option::NLPSolver},
        {"IpoptDerivativeTest", Option::IpoptDerivativeTest},
        {"WarmStart", Option::WarmStart},
        {"QP", Option::QP},
    };
    auto it = map.find(name);
    if (it != map.end()) return it->second;
    return std::nullopt;
}

// === Option Class ===

NLPSolverSettings::NLPSolverSettings(int argc, char** argv) {
    // start with defaults
    settings = default_settings;

    // parse CLI arguments of form: --Key=Value
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg.rfind("--", 0) == 0) {
            size_t eq = arg.find('=');
            if (eq != std::string::npos) {
                std::string key = arg.substr(2, eq - 2);
                std::string val = arg.substr(eq + 1);
                auto maybe_option = option_from_string(key);
                if (maybe_option) {
                    // Try to parse into correct type based on default
                    const auto& default_val = default_settings.at(*maybe_option);
                    if (std::holds_alternative<f64>(default_val)) {
                        settings[*maybe_option] = std::stod(val);
                    } else if (std::holds_alternative<int>(default_val)) {
                        settings[*maybe_option] = std::stoi(val);
                    } else if (std::holds_alternative<bool>(default_val)) {
                        settings[*maybe_option] = (val == "true" || val == "1");
                    } else {
                        settings[*maybe_option] = val;
                    }
                }
            }
        }
    }
}

void NLPSolverSettings::print() const {
    FixedTableFormat<2> ftf = {
        {25,            15},
        {Align::Center, Align::Center}
    };

    Log::start_module(ftf, "NLP Solver options");
    Log::row(ftf, "Option", "Value");
    Log::dashes(ftf);

    for (const auto& [option, value] : settings) {
        std::string val_str = std::visit([](const auto& v) -> std::string {
            using T = std::decay_t<decltype(v)>;

            if constexpr (std::is_same_v<T, bool>) {
                return v ? "true" : "false";
            } else if constexpr (std::is_same_v<T, std::string>) {
                return v;
            } else if constexpr (std::is_same_v<T, f64>) {
                return std::to_string(v);
            } else if constexpr (std::is_same_v<T, int>) {
                return std::to_string(v);
            } else if constexpr (std::is_same_v<T, HessianOption>) {
                switch (v) {
                    case HessianOption::Exact: return "Exact";
                    case HessianOption::LBFGS: return "LBFGS";
                    case HessianOption::Const: return "Const";
                    default: return "<invalid HessianOption>";
                }
            } else if constexpr (std::is_same_v<T, LinearSolverOption>) {
                switch (v) {
                    case LinearSolverOption::MUMPS: return "MUMPS";
                    case LinearSolverOption::MA27:  return "MA27";
                    case LinearSolverOption::MA57:  return "MA57";
                    case LinearSolverOption::MA77:  return "MA77";
                    case LinearSolverOption::MA86:  return "MA86";
                    case LinearSolverOption::MA97:  return "MA97";
                    default: return "<invalid LinearSolverOption>";
                }
            } else if constexpr (std::is_same_v<T, NLPSolverOption>) {
                switch (v) {
                    case NLPSolverOption::Ipopt: return "Ipopt";
                    default: return "<invalid NLPSolverOption>";
                }
            } else {
                return "<unknown>";
            }
        }, value);

        Log::row(ftf, to_string(option), val_str);
    }

    Log::dashes_ln(ftf);
}

void NLPSolverSettings::set(Option option, const OptionValue& value) {
    settings[option] = value;
}

const OptionValue& NLPSolverSettings::get(Option option) const {
    auto it = settings.find(option);
    if (it != settings.end()) return it->second;
    return default_settings.at(option);
}

template<typename T>
T NLPSolverSettings::get_or_default(Option option) const {
    auto it = settings.find(option);
    if (it != settings.end()) {
        if (auto val = std::get_if<T>(&it->second))
            return *val;
    }
    return std::get<T>(default_settings.at(option));
}

bool NLPSolverSettings::option_is_true(Option option) const {
    auto val = get_or_default<bool>(option);
    return val;
}

bool NLPSolverSettings::option_matches(Option option, const std::string& str) const {
    if (auto val = std::get_if<std::string>(&get(option)))
        return *val == str;
    return false;
}

// Explicit instantiations for template
template int NLPSolverSettings::get_or_default<int>(Option) const;
template f64 NLPSolverSettings::get_or_default<f64>(Option) const;
template bool NLPSolverSettings::get_or_default<bool>(Option) const;
template std::string NLPSolverSettings::get_or_default<std::string>(Option) const;
template HessianOption NLPSolverSettings::get_or_default<HessianOption>(Option) const;
template LinearSolverOption NLPSolverSettings::get_or_default<LinearSolverOption>(Option) const;
template NLPSolverOption NLPSolverSettings::get_or_default<NLPSolverOption>(Option) const;

} // namespace NLP
