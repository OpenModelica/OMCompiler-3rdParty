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

#include "strategies.h"
#include "gdop.h"

// TODO: add doxygen everywhere

namespace GDOP {

// ==================== no-op strategies ====================

// no simulation available
std::unique_ptr<Trajectory> NoSimulation::operator()(const ControlTrajectory& controls, const FixedVector<f64>& parameters,
                                                     int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) {
    Log::warning("No Simulation strategy set: returning nullptr.");
    return nullptr;
}

// no simulation step available
void NoSimulationStep::activate(const ControlTrajectory& controls, const FixedVector<f64>& parameters) {}

void NoSimulationStep::reset() {}

std::unique_ptr<Trajectory> NoSimulationStep::operator()(f64* x_start_values, f64 start_time, f64 stop_time) {
    Log::warning("No SimulationStep strategy set: returning nullptr.");
    return nullptr;
}

// no mesh refinement available
void NoMeshRefinement::reset(const GDOP& gdop) {}

std::unique_ptr<MeshUpdate> NoMeshRefinement::operator()(const Mesh& mesh, const PrimalDualTrajectory& trajectory) {
    Log::warning("No MeshRefinement strategy set: returning nullptr.");
    return nullptr;
}

// no emitter
int NoEmitter::operator()(const PrimalDualTrajectory& trajectory) {
    Log::warning("No Emitter strategy set: returning -1.");
    return -1;
}

// no verifier
bool NoVerifier::operator()(const GDOP& gdop, const PrimalDualTrajectory& trajectory) {
    Log::warning("No Verifier strategy set: returning false.");
    return false;
}

// no scaling
std::shared_ptr<NLP::Scaling> NoScalingFactory::operator()(const GDOP& gdop) {
    Log::warning("No ScalingFactory strategy set: fallback to NoScalingFactory.");
    return std::make_shared<NLP::NoScaling>(NLP::NoScaling());
};

// ==================== non no-op strategies ====================

// default initialization (is not really proper, but an implementation)
std::unique_ptr<PrimalDualTrajectory> ConstantInitialization::operator()(const GDOP& gdop) {
    Log::warning("No Initialization strategy set: fallback to ConstantInitialization.");

    const auto& problem = gdop.get_problem();

    const int x_size = problem.pc->x_size;
    const int u_size = problem.pc->u_size;
    const int p_size = problem.pc->p_size;

    // time vector with start and end times
    std::vector<f64> t = { 0.0, gdop.get_mesh().tf };

    std::vector<std::vector<f64>> x_guess(x_size);
    std::vector<std::vector<f64>> u_guess(u_size);
    std::vector<f64>              p_guess(p_size, 0.0);  // TODO: PARAMETERS implement

    InterpolationMethod interpolation = InterpolationMethod::LINEAR;

    // fill state guesses
    for (int x = 0; x < x_size; x++) {
        auto x0_opt = problem.pc->x0_fixed[x];
        auto xf_opt = problem.pc->xf_fixed[x];

        if (x0_opt && xf_opt) {
            // initial and final fixed: linear interpolation
            x_guess[x] = { *x0_opt, *xf_opt };
        } else if (x0_opt) {
            // initial fixed: constant at initial
            x_guess[x] = { *x0_opt, *x0_opt };
        } else if (xf_opt) {
            // final fixed: constant at final
            x_guess[x] = { *xf_opt, *xf_opt };
        } else {
            // nothing fixed: use midpoint of bounds or zero if unbounded
            f64 val = 0.0;
            if (problem.pc->x_bounds[x].has_lower() && problem.pc->x_bounds[x].has_upper()) {
                val = 0.5 * (problem.pc->x_bounds[x].lb + problem.pc->x_bounds[x].ub);
            }
            x_guess[x] = { val, val };
        }
    }

    // fill control guesses: use midpoint of bounds or zero
    for (int u = 0; u < u_size; u++) {
        f64 val = 0.0;
        if (problem.pc->u_bounds[u].has_lower() && problem.pc->u_bounds[u].has_upper()) {
            val = 0.5 * (problem.pc->u_bounds[u].lb + problem.pc->u_bounds[u].ub);
        }
        u_guess[u] = { val, val };
    }

    // fill parameter guesses: use midpoint of bounds or zero
    for (int p = 0; p < p_size; p++) {
        f64 val = 0.0;
        if (problem.pc->p_bounds[p].has_lower() && problem.pc->p_bounds[p].has_upper()) {
            val = 0.5 * (problem.pc->p_bounds[p].lb + problem.pc->p_bounds[p].ub);
        }
        p_guess[p] = val;
    }

    return std::make_unique<PrimalDualTrajectory>(std::make_unique<Trajectory>(t, x_guess, u_guess, p_guess, interpolation));
}

RadauIntegratorSimulation::RadauIntegratorSimulation(Dynamics& dynamics) : dynamics(dynamics) {}

std::unique_ptr<Trajectory> RadauIntegratorSimulation::operator()(
    const ControlTrajectory& controls,
    const FixedVector<f64>& parameters,
    int num_steps,
    f64 start_time,
    f64 stop_time,
    f64* x_start_values)
{
    dynamics.allocate();

    auto ode_f_fn = [this](const f64* x, const f64* u, const f64* p, f64 t, f64* f, void* user_data) {
        return this->dynamics.eval(x, u, p, t, f, user_data);
    };

    auto ode_jac_fn = [this](const f64* x, const f64* u, const f64* p, f64 t, f64* dfdx, void* user_data) {
        return this->dynamics.jac(x, u, p, t, dfdx, user_data);
    };

    auto res = ::Simulation::RadauBuilder()
                .interval(start_time, stop_time, num_steps)
                .states(dynamics.pc.x_size, x_start_values)
                .control(&controls)
                .params(parameters.int_size(), parameters.raw())
                .ode(ode_f_fn)
                .jacobian(ode_jac_fn, dynamics.jac_pattern)
                .radau_scheme(::Simulation::RadauScheme::ADAPTIVE)
                .radau_tol(1e-10, 1e-10)
                .radau_max_it(5e6)
                .build()
                .simulate();

    dynamics.free();

    return res;
}

RadauIntegratorSimulationStep::RadauIntegratorSimulationStep(Dynamics& dynamics) : dynamics(dynamics) {}

void RadauIntegratorSimulationStep::activate(
    const ControlTrajectory& controls,
    const FixedVector<f64>& parameters)
{
    auto ode_f_fn = [this](const f64* x, const f64* u, const f64* p, f64 t, f64* f, void* user_data) {
        return this->dynamics.eval(x, u, p, t, f, user_data);
    };

    auto ode_jac_fn = [this](const f64* x, const f64* u, const f64* p, f64 t, f64* dfdx, void* user_data) {
        return this->dynamics.jac(x, u, p, t, dfdx, user_data);
    };

    auto builder = ::Simulation::RadauBuilder()
                        .states(dynamics.pc.x_size)
                        .control(&controls)
                        .params(parameters.int_size(), parameters.raw())
                        .ode(ode_f_fn)
                        .jacobian(ode_jac_fn, dynamics.jac_pattern)
                        .radau_scheme(::Simulation::RadauScheme::ADAPTIVE)
                        .radau_tol(1e-10, 1e-10)
                        .radau_max_it(5e6);

    integrator = std::make_unique<::Simulation::RadauIntegrator>(builder.build());

    dynamics.allocate();
}

void RadauIntegratorSimulationStep::reset()
{
    integrator = nullptr; // delete captured integrator
    dynamics.free();      // free allocated structures from user callbacks
}

std::unique_ptr<Trajectory> RadauIntegratorSimulationStep::operator()(
    f64* x_start_values,
    f64 start_time,
    f64 stop_time)
{
    if (!integrator) {
        Log::error("RadauIntegrator has not been allocated: returning nullptr.");
        return nullptr;
    }

    return integrator->simulate(x_start_values, start_time, stop_time, 1);
}

std::vector<f64> LinearInterpolation::operator()(
    const Mesh& old_mesh,
    const Mesh& new_mesh,
    const std::vector<f64>& values)
{
    std::vector<f64> old_t = old_mesh.get_flat_t();
    std::vector<f64> new_t = new_mesh.get_flat_t();

    return interpolate_linear_single(old_t, values, new_t);
}

// TODO: reduce overhead in interpolation + copies, etc. if necessary, make this without aux allocation of new_t and old_t
// interpolate trajectory to new mesh with simple linear interpolation
InterpolationRefinedInitialization::InterpolationRefinedInitialization(std::shared_ptr<Interpolation> interpolation_,
                                                                       bool interpolate_primals_,
                                                                       bool interpolate_costates_constraints_,
                                                                       bool interpolate_costates_bounds_)
    : interpolation(interpolation_),
      interpolate_primals(interpolate_primals_),
      interpolate_costates_constraints(interpolate_costates_constraints_),
      interpolate_costates_bounds(interpolate_costates_bounds_) {}

// TODO: refactor this. Can we unify the interpolations even further, now we also need to interpolate controls better at callbacks...
// maybe we should do just 1 grant linear interpolator, and 1 polynomial interpolator class, which implement all the cases? think about this
std::unique_ptr<PrimalDualTrajectory> InterpolationRefinedInitialization::operator()(const Mesh& old_mesh,
                                                                                     const Mesh& new_mesh,
                                                                                     const PrimalDualTrajectory& trajectory)
{
    std::unique_ptr<Trajectory> new_primals              = nullptr;
    std::unique_ptr<CostateTrajectory> new_costates      = nullptr;
    std::unique_ptr<Trajectory> new_costate_bounds_lower = nullptr;
    std::unique_ptr<Trajectory> new_costate_bounds_upper = nullptr;

    // === interpolation of primal variables onto new_mesh ===
    if (interpolate_primals) {
        auto const& old_primals = trajectory.primals;
        new_primals = std::make_unique<Trajectory>();

        // copy time grid
        new_primals->t = new_mesh.get_flat_t();

        // interpolate states
        new_primals->x.resize(old_primals->x.size());
        for (size_t x_index = 0; x_index < old_primals->x.size(); x_index++) {
            new_primals->x[x_index] = (*interpolation)(old_mesh, new_mesh, old_primals->x[x_index]);
        }

        // interpolate controls
        new_primals->u.resize(old_primals->u.size());
        for (size_t u_index = 0; u_index < old_primals->u.size(); u_index++) {
            new_primals->u[u_index] = (*interpolation)(old_mesh, new_mesh, old_primals->u[u_index]);
        }

        // copy parameters
        new_primals->p.resize(old_primals->p.size());
        for (size_t p_index = 0; p_index < old_primals->p.size(); p_index++) {
            new_primals->p[p_index] = old_primals->p[p_index];
        }

        new_primals->inducing_mesh = new_mesh.shared_from_this();
    }

    // === interpolation of duals / costates onto new_mesh ===
    if (interpolate_costates_constraints) {
        auto const& old_costates = trajectory.costates;
        new_costates = std::make_unique<CostateTrajectory>();

        // copy time grid
        new_costates->t = new_mesh.get_flat_t();

        // interpolate lambda_f
        new_costates->costates_f.resize(old_costates->costates_f.size());
        for (size_t f_index = 0; f_index < old_costates->costates_f.size(); f_index++) {
            new_costates->costates_f[f_index] = (*interpolation)(old_mesh, new_mesh, old_costates->costates_f[f_index]);
        }

        // interpolate lambda_g
        new_costates->costates_g.resize(old_costates->costates_g.size());
        for (size_t g_index = 0; g_index < old_costates->costates_g.size(); g_index++) {
            new_costates->costates_g[g_index] = (*interpolation)(old_mesh, new_mesh, old_costates->costates_g[g_index]);
        }

        // copy lambda_r
        new_costates->costates_r.resize(old_costates->costates_r.size());
        for (size_t r_index = 0; r_index < old_costates->costates_r.size(); r_index++) {
            new_costates->costates_r[r_index] = old_costates->costates_r[r_index];
        }

        new_costates->inducing_mesh = new_mesh.shared_from_this();
    }

    // === interpolation of primal variables onto new_mesh ===
    if (interpolate_costates_bounds) {
        new_costate_bounds_lower = std::make_unique<Trajectory>();
        new_costate_bounds_upper = std::make_unique<Trajectory>();

        std::vector<const Trajectory*> old_duals = { trajectory.lower_costates.get(), trajectory.upper_costates.get() };
        std::vector<Trajectory*> new_duals = { new_costate_bounds_lower.get(), new_costate_bounds_upper.get() };

        for (short bound_idx = 0; bound_idx < 2; bound_idx++) {
            const Trajectory& old_dual = *old_duals[bound_idx];
            Trajectory& new_dual = *new_duals[bound_idx];

            // copy time grid
            new_dual.t = new_mesh.get_flat_t();

            // interpolate state-bound duals
            new_dual.x.resize(old_dual.x.size());
            for (size_t x_index = 0; x_index < old_dual.x.size(); x_index++) {
                new_dual.x[x_index] = (*interpolation)(old_mesh, new_mesh, old_dual.x[x_index]);
            }

            // interpolate control-bound duals
            new_dual.u.resize(old_dual.u.size());
            for (size_t u_index = 0; u_index < old_dual.u.size(); u_index++) {
                new_dual.u[u_index] = (*interpolation)(old_mesh, new_mesh, old_dual.u[u_index]);
            }

            new_dual.p = old_dual.p;

            new_dual.inducing_mesh = new_mesh.shared_from_this();
        }

    }

    return std::make_unique<PrimalDualTrajectory>(std::move(new_primals),
                                                  std::move(new_costates),
                                                  std::move(new_costate_bounds_lower),
                                                  std::move(new_costate_bounds_upper));
}

// proper simulation-based initialization strategy
SimulationInitialization::SimulationInitialization(std::shared_ptr<Initialization> initialization,
                                                   std::shared_ptr<Simulation> simulation)
  : initialization(initialization), simulation(simulation) {}

std::unique_ptr<PrimalDualTrajectory> SimulationInitialization::operator()(const GDOP& gdop) {
    auto simple_guess         = (*initialization)(gdop);                                             // call simple, e.g. constant guess
    auto& simple_guess_primal = simple_guess->primals;
    auto extracted_controls   = simple_guess_primal->copy_extract_controls();                        // extract controls from the guess
    auto extracted_parameters = FixedVector<f64>(simple_guess_primal->p);
    auto exctracted_x0        = simple_guess_primal->extract_initial_states();                       // extract x(t_0) from the guess
    auto simulated_guess      = (*simulation)(extracted_controls, extracted_parameters,              // perform simulation using the controls and gdop config
                                              gdop.get_mesh().node_count, 0.0,
                                              gdop.get_mesh().tf, exctracted_x0.raw());
    auto interpolated_sim     = simulated_guess->interpolate_onto_mesh(gdop.get_mesh()); // interpolate simulation to current mesh + collocation
    return std::make_unique<PrimalDualTrajectory>(std::make_unique<Trajectory>(interpolated_sim));
}

// csv emit
CSVEmitter::CSVEmitter(std::string filename, bool write_header) : filename(filename), write_header(write_header) {}

int CSVEmitter::operator()(const PrimalDualTrajectory& trajectory) { return trajectory.primals->to_csv(filename, write_header); }

// print emit
int PrintEmitter::operator()(const PrimalDualTrajectory& trajectory) { trajectory.primals->print_table(); return 0; }

// simulation-based verification
SimulationVerifier::SimulationVerifier(std::shared_ptr<Simulation> simulation,
                                       Linalg::Norm norm,
                                       FixedVector<f64>&& tolerances)
    : simulation(simulation), norm(norm), tolerances(std::move(tolerances)) {}

bool SimulationVerifier::operator()(const GDOP& gdop, const PrimalDualTrajectory& trajectory) {
    auto& trajectory_primal   = trajectory.primals;
    auto extracted_controls   = trajectory_primal->copy_extract_controls();   // extract controls from the trajectory
    auto extracted_parameters = FixedVector<f64>(trajectory_primal->p);       // extract parameters from the trajectory

    extracted_controls.interpolation = InterpolationMethod::POLYNOMIAL;
    auto exctracted_x0      = trajectory_primal->extract_initial_states();    // extract x(t_0) from the trajectory

    // perform simulation using the controls, gdop config and a high number of nodes
    int  high_node_count    = 1 * gdop.get_mesh().node_count;

    auto simulation_result  = (*simulation)(extracted_controls, extracted_parameters, high_node_count,
                                            0.0, gdop.get_mesh().tf, exctracted_x0.raw());

    // result of high resolution simulation is interpolated onto lower resolution mesh
    auto interpolated_opt   = trajectory_primal->interpolate_polynomial_onto_grid(simulation_result->t);

    // calculate errors for each state in given norm (between provided and simulated states)
    auto errors             = interpolated_opt.state_errors(*simulation_result, norm);

    bool is_valid = true;

    FixedTableFormat<4> ftf = {{7,             13,            11,            4},
                               {Align::Center, Align::Center, Align::Center, Align::Center}};

    Log::start_module(ftf, "Simulation-Based Verification");

    Log::row(ftf, "State", fmt::format("Error [{}]", Linalg::norm_to_string(norm)), "Tolerance", "Pass");
    Log::dashes(ftf);

    for (size_t x_idx = 0; x_idx < trajectory_primal->x.size(); x_idx++) {
        f64 tol = tolerances[x_idx];
        f64 err = errors[x_idx];
        bool pass = (err <= tol);

        Log::row(ftf,
            fmt::format("x[{}]", x_idx),
            fmt::format("{:.3e}", err),
            fmt::format("{:.3e}", tol),
            pass ? "PASS" : "FAIL");

        if (!pass) {
            is_valid = false;
        }
    }

    Log::dashes(ftf);

    if (is_valid) {
        Log::success("All state errors within tolerances.");
    } else {
        Log::warning("One or more state errors exceeded tolerances.");
    }

    Log::dashes_ln(ftf);

    return is_valid;
}

// interpolate trajectory to new mesh with collocation scheme - polynomial interpolation
std::vector<f64> PolynomialInterpolation::operator()(const Mesh& old_mesh,
                                                     const Mesh& new_mesh,
                                                     const std::vector<f64>& values) {
    const int grid_size = new_mesh.node_count + 1;
    std::vector<f64> new_values(grid_size);

    // === t = 0.0 ===
    new_values[0] = values[0];

    // === t = t_{i,j} ===
    int global_grid_index = 1;
    int current_old_interval = 0;
    int offset = 0;

    for (int i = 0; i < new_mesh.intervals; i++) {
        int scheme_new = new_mesh.nodes[i];

        for (int j = 0; j < scheme_new; j++) {
            f64 t_query = new_mesh.t[i][j];

            while (current_old_interval + 1 < old_mesh.intervals &&
                old_mesh.grid[current_old_interval + 1] < t_query) {
                offset += old_mesh.nodes[current_old_interval];
                current_old_interval++;
            }

            const int stride      = 1;
            const int old_p_order = old_mesh.nodes[current_old_interval];
            const f64* values_i   = values.data() + offset;

            f64 t_start = old_mesh.grid[current_old_interval];
            f64 t_end   = old_mesh.grid[current_old_interval + 1];

            new_values[global_grid_index] = fLGR::interpolate(
                old_p_order, true, values_i, stride,
                t_start, t_end, t_query
            );

            global_grid_index++;
        }
    }

    assert(global_grid_index == grid_size);

    return new_values;
}

L2BoundaryNorm::L2BoundaryNorm(int max_phase_one_iterations, int max_phase_two_iterations, f64 phase_two_level)
    : max_phase_one_iterations(max_phase_one_iterations),
      max_phase_two_iterations(max_phase_two_iterations),
      phase_two_level(phase_two_level) {}

// L2BoundaryNorm reset for next call
void L2BoundaryNorm::reset(const GDOP& gdop) {
    phase_one_iteration = 0;
    phase_two_iteration = 0;

    // on-interval
    mesh_size_zero = gdop.get_mesh().intervals;

    // corner
    CTOL_1 = FixedVector<f64>(gdop.get_problem().pc->u_size);
    CTOL_2 = FixedVector<f64>(gdop.get_problem().pc->u_size);
    for (size_t i = 0; i < CTOL_1.size(); i++) { CTOL_1[i] = 0.1; }
    for (size_t i = 0; i < CTOL_2.size(); i++) { CTOL_2[i] = 0.1; }
}

// L2BoundaryNorm mesh refinement algorithm
std::unique_ptr<MeshUpdate> L2BoundaryNorm::operator()(const Mesh& mesh, const PrimalDualTrajectory& trajectory) {
    auto& trajectory_primal = trajectory.primals;

    // constant degree for all intervals
    const int p = mesh.nodes[0];

    bool terminated = true;
    FixedVector<f64> new_grid;

    // L2BoundaryNorm
    if (phase_one_iteration < max_phase_one_iterations) {
        // phase I - full bisection case
        terminated = false;
        new_grid = FixedVector<f64>(2 * mesh.grid.size() - 1);
        for (int i = 0; i < mesh.grid.int_size() - 1; i++) {
            new_grid[2 * i] = mesh.grid[i];
            new_grid[2 * i + 1] = 0.5 * (mesh.grid[i] + mesh.grid[i + 1]);
        }
        new_grid[2 * mesh.grid.size() - 2] = mesh.tf;
        phase_one_iteration++;
    }
    else if (phase_two_iteration < max_phase_two_iterations) {
        // phase II - non-smoothness detection
        std::set<f64> set_new_grid;

        // p', p'' in `Enhancing Collocation-Based Dynamic Optimization through Adaptive Mesh Refinement` (Algorithm 3.2)
        FixedVector<f64> p_1(fLGR::get_max_scheme() + 1);
        FixedVector<f64> p_2(fLGR::get_max_scheme() + 1);

        // p'_i(t_{i+1}) and p''_i(t_{i+1}) for boundary condition
        bool has_last_boundary = false;
        f64 p_boundary_1_last_end;
        f64 p_boundary_2_last_end;
        f64 p_boundary_1_this_end;
        f64 p_boundary_2_this_end;

        for (size_t u_idx = 0; u_idx < trajectory_primal->u.size(); u_idx++) {
            auto const& u_vec = trajectory_primal->u[u_idx];

            // compute range of u
            auto [min_it, max_it] = std::minmax_element(u_vec.begin(), u_vec.end());
            f64 u_range = *max_it - *min_it;

            // TODO: what about nominals? What about u == const. ?
            f64 TOL_1 = u_range * pow(10, -phase_two_level) / mesh_size_zero;
            f64 TOL_2 = TOL_1 / 2;

            for (int i = 0; i < mesh.grid.int_size() - 1; i++) {
                int start_index = mesh.acc_nodes[i][0];

                // apply D to the controls on interval i
                fLGR::diff_matrix_multiply(mesh.nodes[i], &u_vec[start_index], p_1.raw());   // p'  = D^(1) * u_hat
                fLGR::diff_matrix_multiply(mesh.nodes[i], p_1.raw(), p_2.raw());             // p'' = D^(2) * u_hat = D^(1) * p'

                // remember last value for next boundary condition
                p_boundary_1_this_end = p_1.back();
                p_boundary_2_this_end = p_2.back();

                // retrieve pointers at index 1 (ignore index 0 as its not needed for the condition) - note index 0 untouched!
                f64 *q_1 = &p_1[1];
                f64 *q_2 = &p_2[1];

                // element-wise square (inplace)
                Linalg::square(mesh.nodes[i], q_1); // q'  = ((p_1')^2,  (p_2')^2,  ..., (p_m')^2)
                Linalg::square(mesh.nodes[i], q_2); // q'' = ((p_1'')^2, (p_2'')^2, ..., (p_m'')^2)

                // calculate the L_2 norm with exact quadrature
                f64 norm_p_1 = sqrt(fLGR::integrate(mesh.nodes[i], q_1)); // norm_p_1 := sqrt(b^T * q')
                f64 norm_p_2 = sqrt(fLGR::integrate(mesh.nodes[i], q_2)); // norm_p_2 := sqrt(b^T * q'')

                // === on-interval condition ===
                if (norm_p_1 > TOL_1 || norm_p_2 > TOL_2) {
                    terminated = false;

                    // left interval
                    if (i > 0) {
                        set_new_grid.insert(0.5 * (mesh.grid[i - 1] + mesh.grid[i]));
                    }

                    // this / center interval
                    set_new_grid.insert(0.5 * (mesh.grid[i] + mesh.grid[i + 1]));

                    // right interval
                    if (i < mesh.grid.int_size() - 2) {
                        set_new_grid.insert(0.5 * (mesh.grid[i + 1] + mesh.grid[i + 2]));
                    }
                }

                // === boundary condition: Compare end of i-1 with start of i ===
                if (has_last_boundary) {
                    f64 p1_i1 = p_1[0];  // this interval, start i - index was untouched in on-interval computation
                    f64 p2_i1 = p_2[0];

                    // boundary condition plus-1-error
                    f64 ERR_1 = std::abs(p1_i1 - p_boundary_1_last_end) / (1.0 + std::min(std::abs(p1_i1), std::abs(p_boundary_1_last_end)));
                    f64 ERR_2 = std::abs(p2_i1 - p_boundary_2_last_end) / (1.0 + std::min(std::abs(p2_i1), std::abs(p_boundary_2_last_end)));

                    if ((i > 0) && (ERR_1 > CTOL_1[u_idx] || ERR_2 > CTOL_2[u_idx])) {
                        terminated = false;

                        // insert this / center midpoint + left / previous midpoint
                        set_new_grid.insert(0.5 * (mesh.grid[i] + mesh.grid[i + 1]));
                        set_new_grid.insert(0.5 * (mesh.grid[i - 1] + mesh.grid[i]));
                    }
                }

                // save last p_k for next round (last = end of this interval)
                p_boundary_1_last_end = p_boundary_1_this_end;
                p_boundary_2_last_end = p_boundary_2_this_end;
                has_last_boundary = true;

                set_new_grid.insert(mesh.grid[i]);
            }
        }

        set_new_grid.insert(mesh.grid.back());

        // create new_grid if not terminated
        if (!terminated) {
            new_grid = FixedVector<f64>(set_new_grid.begin(), set_new_grid.end());
        }

        phase_two_iteration++;
    }

    if (terminated) {
        return nullptr;
    }

    // set constant polynomial degree
    FixedVector<int> new_nodes(new_grid.size() - 1);
    for (int i = 0; i < new_nodes.int_size(); i++) {
        new_nodes[i] = p;
    }

    return std::make_unique<MeshUpdate>(std::move(new_grid),
                                        std::move(new_nodes));
}

// default strategy collection
Strategies Strategies::default_strategies() {
    Strategies strategies;
    strategies.initialization          = std::make_shared<ConstantInitialization>();
    strategies.simulation              = std::make_shared<NoSimulation>();
    strategies.simulation_step         = std::make_shared<NoSimulationStep>();
    strategies.mesh_refinement         = std::make_shared<NoMeshRefinement>();
    strategies.interpolation           = std::make_shared<LinearInterpolation>();
    strategies.emitter                 = std::make_shared<NoEmitter>();
    strategies.verifier                = std::make_shared<NoVerifier>();
    strategies.scaling_factory         = std::make_shared<NoScalingFactory>();
    strategies.refined_initialization  = std::make_shared<InterpolationRefinedInitialization>(
                                            InterpolationRefinedInitialization(strategies.interpolation, true, true, true));
    return strategies;
};

} // namespace GDOP
