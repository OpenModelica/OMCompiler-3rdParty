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

#ifndef MOO_TRAJECTORY_H
#define MOO_TRAJECTORY_H

#include <map>

#include <base/log.h>
#include <base/mesh.h>
#include <base/linalg.h>
#include <base/fLGR.h>
#include <base/export.h>

enum class InterpolationMethod {
    LINEAR = 0,
    POLYNOMIAL = 1
};

struct DynamicField {
    std::string name;
    std::vector<std::vector<f64>> data;
};

struct StaticField {
    std::string name;
    std::vector<f64> data;
};

struct MOO_EXPORT ControlTrajectory {
    std::vector<f64> t;                       // time grid, monotonic increasing
    std::vector<std::vector<f64>> u;          // u[i][j] = value of i-th control at t[j]
    mutable InterpolationMethod interpolation = InterpolationMethod::LINEAR;

    // optional mesh observer (nullptr if not set)
    std::shared_ptr<const Mesh> inducing_mesh = nullptr;

    // for repeated interpolation, cache last index in this->t or last mesh interval
    mutable size_t last_t_index = 0;
    mutable int last_mesh_interval = 0;

    void interpolate_at(f64 t_query, f64* interpolation_values) const;
    void interpolate_at_linear(f64 t_query, f64* interpolation_values) const;
    void interpolate_at_polynomial(f64 t_query, f64* interpolation_values) const;

    // dumps
    void print_table() const;

    // I/O
    int to_csv(const std::string& filename, bool write_header = true) const;
    static ControlTrajectory from_csv(const std::string& filename);
};

// given some data trajectories t, x(t), u(t), p
struct MOO_EXPORT Trajectory {
    // x[i][j] = x_i(t_j) = x_i at time t[j] 
    std::vector<f64> t;
    std::vector<std::vector<f64>> x;
    std::vector<std::vector<f64>> u;
    std::vector<f64> p;
    mutable InterpolationMethod interpolation = InterpolationMethod::LINEAR;

    // optional mesh observer (nullptr if not set)
    std::shared_ptr<const Mesh> inducing_mesh = nullptr;

    Trajectory() = default;

    Trajectory(std::vector<f64> t, std::vector<std::vector<f64>> x, std::vector<std::vector<f64>> u,
               std::vector<f64> p, InterpolationMethod interpolation = InterpolationMethod::LINEAR,
               const Mesh* inducing_mesh = nullptr)
        : t(t), x(x), u(u), p(p), interpolation(interpolation), inducing_mesh(inducing_mesh ? inducing_mesh->shared_from_this() : nullptr) {}

    Trajectory(const Trajectory& other)
        : t(other.t),
          x(other.x),
          u(other.u),
          p(other.p),
          interpolation(other.interpolation),
          inducing_mesh(other.inducing_mesh) {}

    // create new trajectories based on mesh & collocation
    Trajectory interpolate_onto_mesh(const Mesh& mesh) const;
    Trajectory interpolate_onto_mesh_linear(const Mesh& mesh) const;
    Trajectory interpolate_onto_mesh_polynomial(const Mesh& mesh) const;

    Trajectory interpolate_polynomial_onto_grid(const std::vector<f64>& time_grid) const;

    // extract + copy information from the trajectory
    ControlTrajectory copy_extract_controls() const;
    FixedVector<f64> extract_initial_states() const;
    FixedVector<f64> extract_final_states() const;

    // compare with other trajectory
    FixedVector<f64> state_errors(const Trajectory& other, Linalg::Norm norm) const;
    FixedVector<f64> state_errors_inf_norm(const Trajectory& other) const;
    FixedVector<f64> state_errors_2_norm(const Trajectory& other) const;
    FixedVector<f64> state_errors_1_norm(const Trajectory& other) const;

    // dumps
    void print();
    void print_table();

    // I/O
    int to_csv(const std::string& filename, bool write_header = true) const;
    static Trajectory from_csv(const std::string& filename);
};

// dual trajectory for [costates_f, costates_g]_{ij} constraints, costates_r constraints
struct MOO_EXPORT CostateTrajectory {
    // costates_f[i][j] = costates_f_i(t_j) = costates_f_i at time t[j] 
    std::vector<f64> t;
    std::vector<std::vector<f64>> costates_f;
    std::vector<std::vector<f64>> costates_g;
    std::vector<f64> costates_r;
    mutable InterpolationMethod interpolation;

    // optional mesh observer (nullptr if not set)
    std::shared_ptr<const Mesh> inducing_mesh = nullptr;

    CostateTrajectory() = default;

    CostateTrajectory(std::vector<f64> t, std::vector<std::vector<f64>> costates_f, std::vector<std::vector<f64>> costates_g,
                      std::vector<f64> costates_r, InterpolationMethod interpolation = InterpolationMethod::LINEAR,
                      const Mesh* inducing_mesh = nullptr)
        : t(t),
          costates_f(costates_f),
          costates_g(costates_g),
          costates_r(costates_r),
          interpolation(interpolation),
          inducing_mesh(inducing_mesh ? inducing_mesh->shared_from_this() : nullptr) {}

    CostateTrajectory(const CostateTrajectory& other)
        : t(other.t),
          costates_f(other.costates_f),
          costates_g(other.costates_g),
          costates_r(other.costates_r),
          interpolation(other.interpolation),
          inducing_mesh(other.inducing_mesh) {}

    // create new trajectories based on mesh & collocation
    CostateTrajectory interpolate_onto_mesh(const Mesh& mesh) const;
    CostateTrajectory interpolate_onto_mesh_linear(const Mesh& mesh) const;
    CostateTrajectory interpolate_onto_mesh_polynomial(const Mesh& mesh) const;

    CostateTrajectory interpolate_polynomial_onto_grid(const std::vector<f64>& time_grid) const;

    // dumps
    void print();
    void print_table();

    // I/O
    int to_csv(const std::string& filename, bool write_header = true) const;
    static CostateTrajectory from_csv(const std::string& filename);
};

struct MOO_EXPORT PrimalDualTrajectory {
    std::unique_ptr<Trajectory> primals;          // unscaled primal variables: x
    std::unique_ptr<CostateTrajectory> costates;  // transformed constraint multipliers / costates: \hat{lambda_g}
    std::unique_ptr<Trajectory> lower_costates;   // transformed (lb) variable bound multipliers / costates: \hat{x_L}
    std::unique_ptr<Trajectory> upper_costates;   // transformed (ub) variable bound multipliers / costates: \hat{x_U}

    // full constructor
    PrimalDualTrajectory(std::unique_ptr<Trajectory> primals_,
                         std::unique_ptr<CostateTrajectory> costates_,
                         std::unique_ptr<Trajectory> lower_costates_,
                         std::unique_ptr<Trajectory> upper_costates_)
        : primals(std::move(primals_)),
          costates(std::move(costates_)),
          lower_costates(std::move(lower_costates_)),
          upper_costates(std::move(upper_costates_)) {}

    // primals + costates constructor
    PrimalDualTrajectory(std::unique_ptr<Trajectory> primals_,
                         std::unique_ptr<CostateTrajectory> costates_)
        : primals(std::move(primals_)),
          costates(std::move(costates_)),
          lower_costates(nullptr),
          upper_costates(nullptr) {}

    // primals only constructor
    PrimalDualTrajectory(std::unique_ptr<Trajectory> primals_)
        : primals(std::move(primals_)),
          costates(nullptr),
          lower_costates(nullptr),
          upper_costates(nullptr) {}
};

// === shared helpers for Trajectory / CostateTrajectory / ControlTrajectory ===

std::vector<f64> interpolate_polynomial_onto_grid_single(
    const Mesh& mesh,
    const std::vector<f64>& values,
    const std::vector<f64>& time_grid);

std::vector<std::vector<f64>> interpolate_polynomial_onto_grid_multiple(
    const Mesh& mesh,
    const std::vector<std::vector<f64>>& values,
    const std::vector<f64>& time_grid);

std::vector<f64> interpolate_linear_single(
    const std::vector<f64>& t,
    const std::vector<f64>& values,
    const std::vector<f64>& new_t);

std::vector<std::vector<f64>> interpolate_linear_multiple(
    const std::vector<f64>& t,
    const std::vector<std::vector<f64>>& values,
    const std::vector<f64>& new_t);

bool check_time_compatibility(
    const std::vector<f64>& t_vec,
    const std::vector<std::vector<std::vector<f64>>>& fields_to_check,
    const Mesh& mesh);

int write_csv(const std::string& filename,
              const std::vector<f64>& t,
              const std::vector<DynamicField>& dynamics,
              const std::vector<StaticField>& statics,
              bool write_header = true);

void read_csv(
    const std::string& filename,
    std::vector<f64>& t,
    std::map<std::string, std::vector<std::vector<f64>>>& fields,
    std::map<std::string, std::vector<f64>>& static_fields);

void print_trajectory(
    const std::vector<f64>& t,
    const std::vector<std::pair<std::string, std::vector<std::vector<f64>>>>& fields,
    const std::string& static_name,
    const std::vector<f64>& static_field);

void print_trajectory_table(
    const std::vector<f64>& t,
    const std::vector<std::pair<std::string, std::vector<std::vector<f64>>>>& fields,
    const std::string& static_name,
    const std::vector<f64>& static_field,
    const std::string title = "Trajectory Table");

#endif // MOO_TRAJECTORY_H
