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

#include <IpIpoptData.hpp>
#include <IpDenseVector.hpp>
#include <IpSmartPtr.hpp>
#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <IpJournalist.hpp>

#include <base/nlp_structs.h>
#include <base/log.h>

#include "adapter.h"
#include "solver.h"

namespace IpoptSolver {

std::string vformat_no_newline(const char* format, va_list args) {
    va_list args_copy;
    va_copy(args_copy, args);
    int len = vsnprintf(nullptr, 0, format, args_copy);
    va_end(args_copy);

    if (len < 0) {
        return "";
    }

    std::string str(len + 1, '\0');
    vsnprintf(str.data(), len + 1, format, args);
    str.pop_back();

    if (!str.empty() && str.back() == '\n') {
        str.pop_back();
    }

    return str;
}

class LoggerJournal : public Ipopt::Journal {
public:
    LoggerJournal(const std::string& name, Ipopt::EJournalLevel default_level)
    : Ipopt::Journal(name, Ipopt::EJournalLevel::J_ITERSUMMARY) {}

protected:
    void PrintImpl(Ipopt::EJournalCategory category, Ipopt::EJournalLevel level, const char* str) override
    {
        if (str == nullptr) return;
        Log::info(str);
    }

    void PrintfImpl(Ipopt::EJournalCategory category, Ipopt::EJournalLevel level, const char* pformat, va_list ap) override
    {
        std::string formatted_message = vformat_no_newline(pformat, ap);

        if (formatted_message.empty()) return;

        switch (level) {
            case Ipopt::J_ERROR:
                Log::error(formatted_message);
                break;
            case Ipopt::J_STRONGWARNING:
            case Ipopt::J_WARNING:
                Log::warning(formatted_message);
                break;
            case Ipopt::J_SUMMARY:
            case Ipopt::J_ITERSUMMARY:
            case Ipopt::J_DETAILED:
            case Ipopt::J_MOREDETAILED:
            case Ipopt::J_VECTOR:
            case Ipopt::J_MOREVECTOR:
            case Ipopt::J_MATRIX:
            case Ipopt::J_MOREMATRIX:
            case Ipopt::J_ALL:
                Log::info(formatted_message);
                break;
            default:
                Log::info(formatted_message);
                break;
        }
    }

    void FlushBufferImpl() override {}
};

struct IpoptSolverData {
    Ipopt::SmartPtr<IpoptAdapter> adapter;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app;
    Ipopt::SmartPtr<Ipopt::Journal> logger_journal;

    IpoptSolverData(NLP::NLP& nlp)
        : adapter(new IpoptAdapter(nlp)),
          app(IpoptApplicationFactory()),
          logger_journal(new LoggerJournal("LoggerImpl", Ipopt::J_DETAILED))
    {
        app->Jnlst()->DeleteAllJournals();
        app->Jnlst()->AddJournal(logger_journal);
    }
};

IpoptSolver::IpoptSolver(NLP::NLP& nlp, NLP::NLPSolverSettings& solver_settings)
    : NLPSolver(nlp, solver_settings), ipdata{new IpoptSolverData{nlp}}
{
    init_application();
}

IpoptSolver::~IpoptSolver() {
    delete ipdata;
}

// simple wrapper to adapter
void IpoptSolver::optimize() {
    set_settings();

    Ipopt::ApplicationReturnStatus status = ipdata->app->OptimizeTNLP(ipdata->adapter);

    switch (status) {
        case Ipopt::Solve_Succeeded:
            Log::success("[Ipopt Interface] Optimization succeeded.");
            break;

        case Ipopt::Solved_To_Acceptable_Level:
            Log::success("[Ipopt Interface] Optimization succeeded (acceptable).");
            break;

        case Ipopt::Infeasible_Problem_Detected:
            Log::error("[Ipopt Interface] Infeasible problem detected.");
            break;

        case Ipopt::Search_Direction_Becomes_Too_Small:
            Log::warning("[Ipopt Interface] Search direction became too small.");
            break;

        case Ipopt::Diverging_Iterates:
            Log::error("[Ipopt Interface] Diverging iterates.");
            break;

        case Ipopt::User_Requested_Stop:
            Log::warning("[Ipopt Interface] Optimization stopped by user request.");
            break;

        case Ipopt::Feasible_Point_Found:
            Log::error("[Ipopt Interface] Feasible point found.");
            break;

        case Ipopt::Maximum_Iterations_Exceeded:
            Log::warning("[Ipopt Interface] Maximum iterations exceeded.");
            break;

        case Ipopt::Restoration_Failed:
            Log::error("[Ipopt Interface] Restoration failed.");
            break;

        case Ipopt::Error_In_Step_Computation:
            Log::error("[Ipopt Interface] Error in step computation.");
            break;

        case Ipopt::Maximum_CpuTime_Exceeded:
            Log::warning("[Ipopt Interface] Maximum CPU time exceeded.");
            break;

        case Ipopt::Maximum_WallTime_Exceeded:
            Log::warning("[Ipopt Interface] Maximum wall time exceeded.");
            break;

        case Ipopt::Not_Enough_Degrees_Of_Freedom:
            Log::error("[Ipopt Interface] Not enough degrees of freedom.");
            break;

        case Ipopt::Invalid_Problem_Definition:
            Log::error("[Ipopt Interface] Invalid problem definition.");
            break;

        case Ipopt::Invalid_Option:
            Log::error("[Ipopt Interface] Invalid option.");
            break;

        case Ipopt::Invalid_Number_Detected:
            Log::error("[Ipopt Interface] Invalid number detected.");
            break;

        case Ipopt::Unrecoverable_Exception:
            Log::error("[Ipopt Interface] Unrecoverable exception occurred.");
            break;

        case Ipopt::NonIpopt_Exception_Thrown:
            Log::error("[Ipopt Interface] Non-Ipopt exception thrown.");
            break;

        case Ipopt::Insufficient_Memory:
            Log::error("[Ipopt Interface] Insufficient memory.");
            break;

        case Ipopt::Internal_Error:
            Log::error("[Ipopt Interface] Internal error.");
            break;

        default:
            Log::error("[Ipopt Interface] Unknown return status: {}", static_cast<int>(status));
            break;
    }
}

void IpoptSolver::init_application() {
    Ipopt::ApplicationReturnStatus status = ipdata->app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        Log::error("[Ipopt Interface] Error during application initialization.");
        abort();
    }
}

void IpoptSolver::set_settings() {
    // --- termination ---
    ipdata->app->Options()->SetIntegerValue("max_iter", solver_settings.get_or_default<int>(NLP::Option::Iterations));
    ipdata->app->Options()->SetNumericValue("max_cpu_time", solver_settings.get_or_default<f64>(NLP::Option::CPUTime));

    // --- tolerances, step sizes ---
    f64 tol = solver_settings.get_or_default<f64>(NLP::Option::Tolerance);
    ipdata->app->Options()->SetNumericValue("tol", tol);
    ipdata->app->Options()->SetNumericValue("acceptable_tol", tol * 1e3);
    ipdata->app->Options()->SetNumericValue("bound_push", 1e-2);
    ipdata->app->Options()->SetNumericValue("bound_frac", 1e-2);
    ipdata->app->Options()->SetNumericValue("alpha_red_factor", 0.5);

    // --- strategy settings ---
    ipdata->app->Options()->SetStringValue("mu_strategy", "adaptive");
    ipdata->app->Options()->SetStringValue("adaptive_mu_globalization", "kkt-error");
    ipdata->app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    ipdata->app->Options()->SetStringValue("fixed_variable_treatment", "make_parameter");

    // --- hessian options ---
    NLP::HessianOption hess_opt = solver_settings.get_or_default<NLP::HessianOption>(NLP::Option::Hessian);
    switch (hess_opt) {
        case NLP::HessianOption::LBFGS:
            ipdata->app->Options()->SetStringValue("hessian_approximation", "limited-memory");
            break;
        case NLP::HessianOption::Const:
            ipdata->app->Options()->SetStringValue("hessian_constant", "yes");
            break;
        case NLP::HessianOption::Exact:
            ipdata->app->Options()->SetStringValue("hessian_approximation", "exact");
            break;
    }

    // --- warm start ---
    if (solver_settings.option_is_true(NLP::Option::WarmStart)) {
        ipdata->app->Options()->SetStringValue("warm_start_init_point", "yes");
        ipdata->app->Options()->SetStringValue("mu_strategy", "monotone");
        ipdata->app->Options()->SetNumericValue("mu_init", 1e-14);
        ipdata->app->Options()->SetNumericValue("warm_start_bound_push", 1e-8);
        ipdata->app->Options()->SetNumericValue("warm_start_bound_frac", 1e-8);
        ipdata->app->Options()->SetNumericValue("warm_start_slack_bound_push", 1e-8);
        ipdata->app->Options()->SetNumericValue("warm_start_slack_bound_frac", 1e-8);
        ipdata->app->Options()->SetNumericValue("warm_start_mult_bound_push", 1e-8);
    }

    // --- linear solver ---
    NLP::LinearSolverOption linear_solver = solver_settings.get_or_default<NLP::LinearSolverOption>(NLP::Option::LinearSolver);
    switch (linear_solver) {
        case NLP::LinearSolverOption::MUMPS: ipdata->app->Options()->SetStringValue("linear_solver", "mumps"); break;
        case NLP::LinearSolverOption::MA27:  ipdata->app->Options()->SetStringValue("linear_solver", "ma27");  break;
        case NLP::LinearSolverOption::MA57:  ipdata->app->Options()->SetStringValue("linear_solver", "ma57");  break;
        case NLP::LinearSolverOption::MA77:  ipdata->app->Options()->SetStringValue("linear_solver", "ma77");  break;
        case NLP::LinearSolverOption::MA86:  ipdata->app->Options()->SetStringValue("linear_solver", "ma86");  break;
        case NLP::LinearSolverOption::MA97:  ipdata->app->Options()->SetStringValue("linear_solver", "ma97");  break;
    }

    // --- constant derivatives (assumed false for now) ---
    ipdata->app->Options()->SetStringValue("grad_f_constant", "no");
    ipdata->app->Options()->SetStringValue("jac_c_constant", "no");
    ipdata->app->Options()->SetStringValue("jac_d_constant", "no");

    // --- info ---
    ipdata->app->Options()->SetStringValue("timing_statistics", "yes");

    // --- derivative test (optional) ---
    if (solver_settings.option_is_true(NLP::Option::IpoptDerivativeTest)) {
        ipdata->app->Options()->SetStringValue("derivative_test", "second-order");
        ipdata->app->Options()->SetNumericValue("derivative_test_tol", 1e-2);
        ipdata->app->Options()->SetNumericValue("point_perturbation_radius", 0);
    }
}

} // namespace IpoptSolver
