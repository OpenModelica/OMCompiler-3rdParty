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

#include "trajectory.h"

Trajectory Trajectory::interpolate_onto_mesh(const Mesh& mesh) const {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            return interpolate_onto_mesh_linear(mesh);
        case InterpolationMethod::POLYNOMIAL:
            return interpolate_onto_mesh_polynomial(mesh);
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

Trajectory Trajectory::interpolate_onto_mesh_linear(const Mesh& mesh) const {
    Trajectory new_traj;

    std::vector<f64> new_t = mesh.get_flat_t();

    new_traj.t = new_t;
    new_traj.x = interpolate_linear_multiple(t, x, new_t);
    new_traj.u = interpolate_linear_multiple(t, u, new_t);
    new_traj.p = p;
    new_traj.inducing_mesh = mesh.shared_from_this();

    return new_traj;
}

Trajectory Trajectory::interpolate_onto_mesh_polynomial(const Mesh& mesh) const {
    std::vector<f64> mesh_grid = mesh.get_flat_t();
    auto new_traj = interpolate_polynomial_onto_grid(mesh_grid);
    new_traj.inducing_mesh = mesh.shared_from_this();

    return new_traj;
}

Trajectory Trajectory::interpolate_polynomial_onto_grid(const std::vector<f64>& time_grid) const {
    Trajectory new_traj;
    new_traj.t = time_grid;
    new_traj.x = interpolate_polynomial_onto_grid_multiple(*inducing_mesh, x, time_grid);
    new_traj.u = interpolate_polynomial_onto_grid_multiple(*inducing_mesh, u, time_grid);
    new_traj.p = p;
    new_traj.inducing_mesh = nullptr;

    return new_traj;
}

void Trajectory::print() {
    print_trajectory(t, {
        {"x", x},
        {"u", u}
    }, "p", p);
}

int Trajectory::to_csv(const std::string& filename) const {
    return write_trajectory_csv(filename, t, {
        {"x", x},
        {"u", u}
    }, "p", p);
}

FixedVector<f64> Trajectory::extract_initial_states() const {
    FixedVector<f64> x0(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        x0[i] = x[i][0];
    }

    return x0;
}

FixedVector<f64> Trajectory::state_errors_inf_norm(const Trajectory& other) const {
    FixedVector<f64> max_abs_errors(x.size());

    for (size_t x_idx = 0; x_idx < x.size(); x_idx++) {
        const auto& x_traj_1 = x[x_idx];
        const auto& x_traj_2 = other.x[x_idx];

        if (x_traj_1.size() != x_traj_2.size()) {
            LOG_ERROR("State trajectory length mismatch in Trajectory::state_errors_inf_norm.");
            return max_abs_errors;
        }

        f64* max_err = &max_abs_errors[x_idx];
        for (size_t t_idx = 0; t_idx < x_traj_1.size(); t_idx++) {
            f64 diff = std::abs(x_traj_1[t_idx] - x_traj_2[t_idx]);
            if (diff > *max_err) {
                *max_err = diff;
            }
        }
    }

    return max_abs_errors;
}

FixedVector<f64> Trajectory::state_errors_2_norm(const Trajectory& other) const {
    FixedVector<f64> norm_errors(x.size());

    for (size_t x_idx = 0; x_idx < x.size(); x_idx++) {
        const auto& x_traj_1 = x[x_idx];
        const auto& x_traj_2 = other.x[x_idx];

        if (x_traj_1.size() != x_traj_2.size()) {
            LOG_ERROR("State trajectory length mismatch in Trajectory::state_errors_2_norm.");
            return norm_errors;
        }

        f64 sum_sq_diff = 0.0;
        for (size_t t_idx = 0; t_idx < x_traj_1.size(); t_idx++) {
            f64 diff = x_traj_1[t_idx] - x_traj_2[t_idx];
            sum_sq_diff += diff * diff;
        }
        norm_errors[x_idx] = std::sqrt(sum_sq_diff);
    }

    return norm_errors;
}

FixedVector<f64> Trajectory::state_errors_1_norm(const Trajectory& other) const {
    FixedVector<f64> norm_errors(x.size());

    for (size_t x_idx = 0; x_idx < x.size(); x_idx++) {
        const auto& x_traj_1 = x[x_idx];
        const auto& x_traj_2 = other.x[x_idx];

        if (x_traj_1.size() != x_traj_2.size()) {
            LOG_ERROR("State trajectory length mismatch in Trajectory::state_errors_1_norm.");
            return norm_errors;
        }

        f64 sum_abs_diff = 0.0;
        for (size_t t_idx = 0; t_idx < x_traj_1.size(); t_idx++) {
            sum_abs_diff += std::abs(x_traj_1[t_idx] - x_traj_2[t_idx]);
        }
        norm_errors[x_idx] = sum_abs_diff;
    }

    return norm_errors;
}

FixedVector<f64> Trajectory::state_errors(const Trajectory& other, Linalg::Norm norm) const {
    if (this->t.size() != other.t.size()) {
        throw std::runtime_error("Time vector size mismatch in Trajectory::state_errors.");
    }

    switch (norm) {
        case Linalg::Norm::NORM_INF:
            return state_errors_inf_norm(other);
        case Linalg::Norm::NORM_2:
            return state_errors_2_norm(other);
        case Linalg::Norm::NORM_1:
            return state_errors_1_norm(other);
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

// === Control Trajectory ===

ControlTrajectory Trajectory::copy_extract_controls() const {
    ControlTrajectory controls_copy;
    controls_copy.t = t;
    controls_copy.u = u;
    controls_copy.interpolation = interpolation;
    controls_copy.last_t_index = 0;
    controls_copy.last_mesh_interval = 0;
    controls_copy.inducing_mesh = inducing_mesh;

    return controls_copy;
}

void ControlTrajectory::interpolate_at_linear(f64 t_query, f64* interpolation_values) const {
    const size_t t_len = t.size();

    // out of bounds cases
    if (t_query <= t.front()) {
        for (size_t u_idx = 0; u_idx < u.size(); u_idx++) {
            interpolation_values[u_idx] = u[u_idx][0];
        }
        last_t_index = 0;
        return;
    }
    if (t_query >= t.back()) {
        for (size_t u_idx = 0; u_idx < u.size(); u_idx++) {
            interpolation_values[u_idx] = u[u_idx].back();
        }
        last_t_index = t_len - 2; // safe last segment
        return;
    }

    // search for the correct interval using the last_t_index
    size_t i = last_t_index;
    while (i + 1 < t_len && t_query > t[i + 1]) {
        i++;
    }
    while (i > 0 && t_query < t[i]) {
        i--;
    }

    // interval [t[i], t[i+1]]
    f64 t1 = t[i];
    f64 t2 = t[i + 1];
    f64 alpha = (t_query - t1) / (t2 - t1);

    for (size_t k = 0; k < u.size(); ++k) {
        f64 u1 = u[k][i];
        f64 u2 = u[k][i + 1];
        interpolation_values[k] = u1 + alpha * (u2 - u1);
    }

    last_t_index = i; // keep last_index for next call
    return;
}

void ControlTrajectory::interpolate_at_polynomial(f64 t_query, f64* interpolation_values) const {
    // out-of-bounds cases
    if (t_query <= t.front()) {
        for (size_t u_idx = 0; u_idx < u.size(); u_idx++) {
            interpolation_values[u_idx] = u[u_idx][0];
        }
        last_mesh_interval = 0;
        return;
    }
    if (t_query >= t.back()) {
        for (size_t u_idx = 0; u_idx < u.size(); u_idx++) {
            interpolation_values[u_idx] = u[u_idx].back();
        }
        last_mesh_interval = inducing_mesh->intervals - 1;
        return;
    }

    // search for the correct mesh interval using the last_mesh_interval
    int i = last_mesh_interval;
    while (i + 1 < inducing_mesh->intervals && t_query > inducing_mesh->grid[i + 1]) {
        i++;
    }
    while (i > 0 && t_query < inducing_mesh->grid[i]) {
        i--;
    }

    const int mesh_p_order = inducing_mesh->nodes[i];
    const f64 t_start      = inducing_mesh->grid[i];
    const f64 t_end        = inducing_mesh->grid[i + 1];
    const int offset       = inducing_mesh->acc_nodes[i][0];

    for (size_t u_idx = 0; u_idx < u.size(); u_idx++) {
        const f64* values_i = u[u_idx].data() + offset;
        interpolation_values[u_idx] = fLGR::interpolate(
            mesh_p_order, true, values_i, 1,
            t_start, t_end, t_query
        );
    }

    last_mesh_interval = i;
}

void ControlTrajectory::interpolate_at(f64 t_query, f64* interpolation_values) const {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            interpolate_at_linear(t_query, interpolation_values);
            return;
        case InterpolationMethod::POLYNOMIAL:
            if (!inducing_mesh) {
                LOG_WARNING("ControlTrajectory in interpolate_at() not induced from a Mesh. Can't perform polynomial interpolation: fallback to interpolate_at_linear().");
                interpolate_at_linear(t_query, interpolation_values);
            } else {
                interpolate_at_polynomial(t_query, interpolation_values);
            }
            return;
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

// === Dual Trajectory ===

CostateTrajectory CostateTrajectory::interpolate_onto_mesh(const Mesh& mesh) const {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            return interpolate_onto_mesh_linear(mesh);
        case InterpolationMethod::POLYNOMIAL:
            if (!inducing_mesh) {
                LOG_WARNING("CostateTrajectory in interpolate_onto_mesh() not induced from a Mesh. Can't perform polynomial interpolation: fallback to interpolate_onto_mesh_linear().");
                return interpolate_onto_mesh_linear(mesh);
            } else {
                return interpolate_onto_mesh_polynomial(mesh);
            }
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

CostateTrajectory CostateTrajectory::interpolate_onto_mesh_linear(const Mesh& mesh) const {
    CostateTrajectory new_dual;
    std::vector<f64> new_t = mesh.get_flat_t();

    new_dual.t = new_t;
    new_dual.costates_f = interpolate_linear_multiple(t, costates_f, new_t);
    new_dual.costates_g = interpolate_linear_multiple(t, costates_g, new_t);
    new_dual.costates_r = costates_r;
    new_dual.inducing_mesh = mesh.shared_from_this();

    return new_dual;
}

CostateTrajectory CostateTrajectory::interpolate_onto_mesh_polynomial(const Mesh& mesh) const {
    std::vector<f64> new_t = mesh.get_flat_t();

    auto new_dual = interpolate_polynomial_onto_grid(new_t);
    new_dual.inducing_mesh = mesh.shared_from_this();

    return new_dual;
}

CostateTrajectory CostateTrajectory::interpolate_polynomial_onto_grid(const std::vector<f64>& time_grid) const {
    CostateTrajectory new_traj;
    new_traj.t = time_grid;
    new_traj.costates_f = interpolate_polynomial_onto_grid_multiple(*inducing_mesh, costates_f, time_grid);
    new_traj.costates_g = interpolate_polynomial_onto_grid_multiple(*inducing_mesh, costates_g, time_grid);
    new_traj.costates_r = costates_r;
    new_traj.inducing_mesh = nullptr;

    return new_traj;
}

void CostateTrajectory::print() {
    print_trajectory(t, {
        {"costates_f", costates_f},
        {"costates_g", costates_g}
    }, "costates_r", costates_r);
}

int CostateTrajectory::to_csv(const std::string& filename) const {
    return write_trajectory_csv(filename, t, {
        {"costates_f", costates_f},
        {"costates_g", costates_g}
    }, "costates_r", costates_r);
}

// === helpers for Dual and standard Trajectory ===

bool check_time_compatibility(
    const std::vector<f64>& t_vec,
    const std::vector<std::vector<std::vector<f64>>>& fields_to_check,
    const Mesh& mesh)
{
    const f64 tol = 1e-12;

    int expected_size = mesh.node_count + 1;
    if (static_cast<int>(t_vec.size()) != expected_size) {
        LOG_WARNING("Time array is not compatible with given Mesh.");
        return false;
    }

    int time_idx = 1;
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            if (std::abs(t_vec[time_idx++] - mesh.t[i][j]) > tol) {
                LOG_WARNING("Time array is not compatible with given Mesh.");
                return false;
            }
        }
    }

    for (const auto& vect : fields_to_check) {
        for (const auto& field : vect) {
            if (static_cast<int>(field.size()) != static_cast<int>(t_vec.size())) {
                LOG_WARNING("Time array is not compatible with given Mesh.");
                return false;
            }
        }
    }

    return true;
}

// === default vector + time_grid interpolations ===

// interpolate trajectory to new mesh with collocation scheme - polynomial interpolation
std::vector<f64> interpolate_polynomial_onto_grid_single(const Mesh& mesh,
                                                         const std::vector<f64>& values,
                                                         const std::vector<f64>& time_grid) {
    std::vector<f64> new_values(time_grid.size());
    new_values[0] = values[0]; // t = 0

    int offset = 0;
    int curr_mesh_interval = 0;

    for (size_t k = 0; k < time_grid.size(); k++) {
        f64 t_query = time_grid[k];

        while (curr_mesh_interval + 1 < mesh.intervals &&
            mesh.grid[curr_mesh_interval + 1] < t_query) {
            offset += mesh.nodes[curr_mesh_interval];
            curr_mesh_interval++;
        }

        const int stride       = 1;
        const int mesh_p_order = mesh.nodes[curr_mesh_interval];
        const f64* values_i    = values.data() + offset;

        f64 t_start = mesh.grid[curr_mesh_interval];
        f64 t_end   = mesh.grid[curr_mesh_interval + 1];

        new_values[k] = fLGR::interpolate(
            mesh_p_order, true, values_i, stride,
            t_start, t_end, t_query
        );
    }

    return new_values;
}

std::vector<std::vector<f64>> interpolate_polynomial_onto_grid_multiple(const Mesh& mesh,
                                                                        const std::vector<std::vector<f64>>& values,
                                                                        const std::vector<f64>& time_grid)
{
    std::vector<std::vector<f64>> out(values.size());
    for (size_t i = 0; i < values.size(); i++) {
        out[i] = interpolate_polynomial_onto_grid_single(mesh, values[i], time_grid);
    }
    return out;
}

std::vector<f64> interpolate_linear_single(
    const std::vector<f64>& old_t,
    const std::vector<f64>& values,
    const std::vector<f64>& new_t)
{
    std::vector<f64> out_values(new_t.size());

    for (int i = 0; i < int(new_t.size()); i++) {
        f64 t_new = new_t[i];
        auto it = std::lower_bound(old_t.begin(), old_t.end(), t_new);

        if (it == old_t.begin()) {
            out_values[i] = values[0];
        } else if (it == old_t.end()) {
            out_values[i] = values.back();
        } else {
            int idx = std::distance(old_t.begin(), it);
            f64 t1 = old_t[idx - 1];
            f64 t2 = old_t[idx];
            f64 y1 = values[idx - 1];
            f64 y2 = values[idx];
            out_values[i] = y1 + (t_new - t1) * (y2 - y1) / (t2 - t1);
        }
    }

    return out_values;
}

std::vector<std::vector<f64>> interpolate_linear_multiple(
    const std::vector<f64>& old_t,
    const std::vector<std::vector<f64>>& values,
    const std::vector<f64>& new_t)
{
    int fields = int(values.size());
    std::vector<std::vector<f64>> out_values(fields);
    for (int field = 0; field < fields; field++) {
        out_values[field] = interpolate_linear_single(old_t, values[field], new_t);
    }
    return out_values;
}

// === prints ===

void print_trajectory(
    const std::vector<f64>& t,
    const std::vector<std::pair<std::string, std::vector<std::vector<f64>>>>& fields,
    const std::string& static_name,
    const std::vector<f64>& static_field)
{
    auto print_vector = [](const std::string& name, const std::vector<f64>& vec) {
        std::cout << name << " = [";
        for (size_t i = 0; i < vec.size(); ++i) {
            std::cout << vec[i];
            if (i + 1 < vec.size()) std::cout << ", ";
        }
        std::cout << "]\n";
    };

    auto print_matrix = [](const std::string& name, const std::vector<std::vector<f64>>& mat) {
        std::cout << name << " = [\n";
        for (const auto& row : mat) {
            std::cout << "  [";
            for (size_t j = 0; j < row.size(); ++j) {
                std::cout << row[j];
                if (j + 1 < row.size()) std::cout << ", ";
            }
            std::cout << "],\n";
        }
        std::cout << "]\n";
    };

    print_vector("t", t);
    for (const auto& [name, mat] : fields) {
        print_matrix(name, mat);
    }
    print_vector(static_name, static_field);
}

int write_trajectory_csv(
    const std::string& filename,
    const std::vector<f64>& t,
    const std::vector<std::pair<std::string, std::vector<std::vector<f64>>>>& fields,
    const std::string& static_name,
    const std::vector<f64>& static_field) // static_field should ideally have only one value or a set of static values
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "[Warning] Failed to open file for writing: " << filename << "\n";
        return -1;
    }

    // header
    file << "time";
    for (const auto& [name, mat] : fields) {
        for (size_t i = 0; i < mat.size(); ++i) {
            file << "," << name << "[" << i << "]";
        }
    }

    // static field header(s)
    for (size_t i = 0; i < static_field.size(); ++i) {
        file << "," << static_name;
        if (static_field.size() > 1) {
            file << "[" << i << "]";
        }
    }
    file << "\n";

    file << std::setprecision(16);

    for (size_t k = 0; k < t.size(); ++k) {
        file << t[k];
        for (const auto& [_, mat] : fields) {
            for (const auto& series : mat) {
                if (k < series.size()) {
                    file << "," << series[k];
                } else {
                    file << ",";
                }
            }
        }

        for (size_t i = 0; i < static_field.size(); ++i) {
            file << "," << static_field[i];
        }
        file << "\n";
    }

    file.close();
    return 0;
}
