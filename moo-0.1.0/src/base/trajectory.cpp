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

#include <algorithm>

#include "trajectory.h"

Trajectory Trajectory::interpolate_onto_mesh(const Mesh& mesh) const {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            return interpolate_onto_mesh_linear(mesh);
        case InterpolationMethod::POLYNOMIAL:
            return interpolate_onto_mesh_polynomial(mesh);
        default:
            Log::error("Unknown interpolation method!");
            abort();
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

void Trajectory::print_table() {
    print_trajectory_table(t, {
        {"x", x},
        {"u", u}
    }, "p", p,
    "Trajectory Table");
}

int Trajectory::to_csv(const std::string& filename, bool write_header) const {
    std::vector<DynamicField> dynamics;
    if (!x.empty()) dynamics.push_back({"x", x});
    if (!u.empty()) dynamics.push_back({"u", u});

    std::vector<StaticField> statics;
    if (!p.empty()) statics.push_back({"p", p});

    return write_csv(filename, t, dynamics, statics, write_header);
}

Trajectory Trajectory::from_csv(const std::string& filename) {
    std::vector<f64> t_new;
    std::map<std::string, std::vector<std::vector<f64>>> dynamics;
    std::map<std::string, std::vector<f64>> statics;

    read_csv(filename, t_new, dynamics, statics);

    return Trajectory(
        t_new,
        dynamics.count("x") ? dynamics["x"] : std::vector<std::vector<f64>>{},
        dynamics.count("u") ? dynamics["u"] : std::vector<std::vector<f64>>{},
        statics.count("p") ? statics["p"] : std::vector<f64>{},
        InterpolationMethod::LINEAR
    );
}

FixedVector<f64> Trajectory::extract_initial_states() const {
    FixedVector<f64> x0(x.size());

    for (size_t x_idx = 0; x_idx < x.size(); x_idx++) {
        x0[x_idx] = x[x_idx][0];
    }

    return x0;
}

FixedVector<f64> Trajectory::extract_final_states() const {
    FixedVector<f64> xf(x.size());

    for (size_t x_idx = 0; x_idx < x.size(); x_idx++) {
        xf[x_idx] = x[x_idx].back();
    }

    return xf;
}

FixedVector<f64> Trajectory::state_errors_inf_norm(const Trajectory& other) const {
    FixedVector<f64> max_abs_errors(x.size());

    for (size_t x_idx = 0; x_idx < x.size(); x_idx++) {
        const auto& x_traj_1 = x[x_idx];
        const auto& x_traj_2 = other.x[x_idx];

        if (x_traj_1.size() != x_traj_2.size()) {
            Log::error("State trajectory length mismatch in Trajectory::state_errors_inf_norm.");
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
            Log::error("State trajectory length mismatch in Trajectory::state_errors_2_norm.");
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
            Log::error("State trajectory length mismatch in Trajectory::state_errors_1_norm.");
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
        Log::error("Time vector size mismatch in Trajectory::state_errors!");
        abort();
    }

    switch (norm) {
        case Linalg::Norm::NORM_INF:
            return state_errors_inf_norm(other);
        case Linalg::Norm::NORM_2:
            return state_errors_2_norm(other);
        case Linalg::Norm::NORM_1:
            return state_errors_1_norm(other);
        default:
            Log::error("Unknown interpolation method!");
            abort();
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

void ControlTrajectory::print_table() const {
    print_trajectory_table(t, { {"u", u} }, "", {}, "ControlTrajectory Table");
}

int ControlTrajectory::to_csv(const std::string& filename, bool write_header) const {
    std::vector<DynamicField> dynamics;
    if (!u.empty()) dynamics.push_back({"u", u});
    std::vector<StaticField> statics;
    return write_csv(filename, t, dynamics, statics, write_header);
}

ControlTrajectory ControlTrajectory::from_csv(const std::string& filename) {
    std::vector<f64> t_new;
    std::map<std::string, std::vector<std::vector<f64>>> dynamics;
    std::map<std::string, std::vector<f64>> statics;

    read_csv(filename, t_new, dynamics, statics);

    ControlTrajectory traj;
    traj.t = t_new;
    traj.u = dynamics.count("u") ? dynamics["u"] : std::vector<std::vector<f64>>{};
    traj.interpolation = InterpolationMethod::LINEAR;
    return traj;
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
    size_t t_idx = last_t_index;
    while (t_idx + 1 < t_len && t_query > t[t_idx + 1]) {
        t_idx++;
    }
    while (t_idx > 0 && t_query < t[t_idx]) {
        t_idx--;
    }

    // interval [t[t_idx], t[t_idx+1]]
    f64 t1 = t[t_idx];
    f64 t2 = t[t_idx + 1];
    f64 alpha = (t_query - t1) / (t2 - t1);

    for (size_t u_idx = 0; u_idx < u.size(); u_idx++) {
        f64 u1 = u[u_idx][t_idx];
        f64 u2 = u[u_idx][t_idx + 1];
        interpolation_values[u_idx] = u1 + alpha * (u2 - u1);
    }

    last_t_index = t_idx; // keep last_index for next call
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
    int m_interval = last_mesh_interval;
    while (m_interval + 1 < inducing_mesh->intervals && t_query > inducing_mesh->grid[m_interval + 1]) {
        m_interval++;
    }
    while (m_interval > 0 && t_query < inducing_mesh->grid[m_interval]) {
        m_interval--;
    }

    const int mesh_p_order = inducing_mesh->nodes[m_interval];
    const f64 t_start      = inducing_mesh->grid[m_interval];
    const f64 t_end        = inducing_mesh->grid[m_interval + 1];
    const int offset       = inducing_mesh->acc_nodes[m_interval][0];

    for (size_t u_idx = 0; u_idx < u.size(); u_idx++) {
        const f64* values_i = u[u_idx].data() + offset;
        interpolation_values[u_idx] = fLGR::interpolate(
            mesh_p_order, true, values_i, 1,
            t_start, t_end, t_query
        );
    }

    last_mesh_interval = m_interval;
}

void ControlTrajectory::interpolate_at(f64 t_query, f64* interpolation_values) const {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            interpolate_at_linear(t_query, interpolation_values);
            return;
        case InterpolationMethod::POLYNOMIAL:
            if (!inducing_mesh) {
                Log::warning("ControlTrajectory in interpolate_at() not induced from a Mesh. Can't perform polynomial interpolation: fallback to interpolate_at_linear().");
                interpolate_at_linear(t_query, interpolation_values);
            } else {
                interpolate_at_polynomial(t_query, interpolation_values);
            }
            return;
        default:
            Log::error("Unknown interpolation method!");
            abort();
    }
}

// === Dual Trajectory ===

int CostateTrajectory::to_csv(const std::string& filename, bool write_header) const {
    std::vector<DynamicField> dynamics;
    if (!costates_f.empty()) dynamics.push_back({"λ_f", costates_f});
    if (!costates_g.empty()) dynamics.push_back({"λ_g", costates_g});

    std::vector<StaticField> statics;
    if (!costates_r.empty()) statics.push_back({"λ_r", costates_r});

    return write_csv(filename, t, dynamics, statics, write_header);
}

CostateTrajectory CostateTrajectory::from_csv(const std::string& filename) {
    std::vector<f64> t_new;
    std::map<std::string, std::vector<std::vector<f64>>> dynamics;
    std::map<std::string, std::vector<f64>> statics;

    read_csv(filename, t_new, dynamics, statics);

    return CostateTrajectory(
        t_new,
        dynamics.count("λ_f") ? dynamics["λ_f"] : std::vector<std::vector<f64>>{},
        dynamics.count("λ_g") ? dynamics["λ_g"] : std::vector<std::vector<f64>>{},
        statics.count("λ_r") ? statics["λ_r"] : std::vector<f64>{},
        InterpolationMethod::LINEAR
    );
}

CostateTrajectory CostateTrajectory::interpolate_onto_mesh(const Mesh& mesh) const {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            return interpolate_onto_mesh_linear(mesh);
        case InterpolationMethod::POLYNOMIAL:
            if (!inducing_mesh) {
                Log::warning("CostateTrajectory in interpolate_onto_mesh() not induced from a Mesh. Can't perform polynomial interpolation: fallback to interpolate_onto_mesh_linear().");
                return interpolate_onto_mesh_linear(mesh);
            } else {
                return interpolate_onto_mesh_polynomial(mesh);
            }
        default:
            Log::error("Unknown interpolation method!");
            abort();
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
        {"λ_f", costates_f},
        {"λ_g", costates_g}
    }, "λ_r", costates_r);
}

void CostateTrajectory::print_table() {
    print_trajectory_table(t, {
        {"λ_f", costates_f},
        {"λ_g", costates_g}
    }, "λ_r", costates_r,
    "CostateTrajectory Table");
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
        Log::warning("Time array is not compatible with given Mesh.");
        return false;
    }

    int t_idx = 1;
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            if (std::abs(t_vec[t_idx++] - mesh.t[i][j]) > tol) {
                Log::warning("Time array is not compatible with given Mesh.");
                return false;
            }
        }
    }

    for (const auto& vect : fields_to_check) {
        for (const auto& field : vect) {
            if (static_cast<int>(field.size()) != static_cast<int>(t_vec.size())) {
                Log::warning("Time array is not compatible with given Mesh.");
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
    for (size_t var_idx = 0; var_idx < values.size(); var_idx++) {
        out[var_idx] = interpolate_polynomial_onto_grid_single(mesh, values[var_idx], time_grid);
    }
    return out;
}

std::vector<f64> interpolate_linear_single(
    const std::vector<f64>& old_t,
    const std::vector<f64>& values,
    const std::vector<f64>& new_t)
{
    std::vector<f64> out_values(new_t.size());

    for (int t_idx = 0; t_idx < int(new_t.size()); t_idx++) {
        f64 t_new = new_t[t_idx];
        auto it = std::lower_bound(old_t.begin(), old_t.end(), t_new);

        if (it == old_t.begin()) {
            out_values[t_idx] = values[0];
        } else if (it == old_t.end()) {
            out_values[t_idx] = values.back();
        } else {
            int idx = std::distance(old_t.begin(), it);
            f64 t1 = old_t[idx - 1];
            f64 t2 = old_t[idx];
            f64 y1 = values[idx - 1];
            f64 y2 = values[idx];
            out_values[t_idx] = y1 + (t_new - t1) * (y2 - y1) / (t2 - t1);
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
        std::string s = name + " = [";
        for (size_t i = 0; i < vec.size(); ++i) {
            s += fmt::format("{}", vec[i]);
            if (i + 1 < vec.size()) s += ", ";
        }
        s += "]";
        Log::info(s);
    };

    auto print_matrix = [](const std::string& name, const std::vector<std::vector<f64>>& mat) {
        std::string s = name + " = [\n";
        for (const auto& row : mat) {
            s += "  [";
            for (size_t j = 0; j < row.size(); ++j) {
                s += fmt::format("{}", row[j]);
                if (j + 1 < row.size()) s += ", ";
            }
            s += "],\n";
        }
        s += "]";
        Log::info(s);
    };

    print_vector("time", t);
    for (const auto& [name, mat] : fields) {
        print_matrix(name, mat);
    }
    print_vector(static_name, static_field);
}

void print_trajectory_table(
    const std::vector<f64>& t,
    const std::vector<std::pair<std::string, std::vector<std::vector<f64>>>>& fields,
    const std::string& static_name,
    const std::vector<f64>& static_field,
    const std::string title)
{
    size_t N = t.size();

    std::vector<std::string> col_names;
    col_names.push_back("time");
    for (const auto& [name, mat] : fields) {
        for (size_t var_idx = 0; var_idx < mat.size(); var_idx++) {
            col_names.push_back(fmt::format("{}[{}]", name, var_idx + 1));
        }
    }

    std::vector<int> widths(col_names.size(), 8);
    for (size_t col = 0; col < col_names.size(); col++) {
        widths[col] = std::max(widths[col], static_cast<int>(col_names[col].size()));
    }
    for (size_t t_idx = 0; t_idx < N; t_idx++) {
        widths[0] = std::max(widths[0], static_cast<int>(fmt::format("{:.6e}", t[t_idx]).size()));
        size_t col = 1;
        for (const auto& [_, mat] : fields) {
            for (size_t var_idx = 0; var_idx < mat.size(); var_idx++, col++) {
                widths[col] = std::max(widths[col], static_cast<int>(fmt::format("{:.6e}", mat[var_idx][t_idx]).size()));
            }
        }
    }

    std::vector<Align> aligns(col_names.size(), Align::Right);
    TableFormat tf(widths, aligns);

    Log::start_module(tf, title);
    Log::row(tf, col_names);
    Log::dashes(tf);

    for (size_t t_idx = 0; t_idx < N; t_idx++) {
        std::vector<std::string> row;
        row.push_back(fmt::format("{:.6e}", t[t_idx]));
        for (const auto& [_, mat] : fields) {
            for (size_t var_idx = 0; var_idx < mat.size(); var_idx++) {
                row.push_back(fmt::format("{:.6e}", mat[var_idx][t_idx]));
            }
        }
        Log::row(tf, row);
    }
    Log::dashes_ln(tf);

    if (static_name != "") {
        std::string header1 = static_name;
        std::string header2 = "Value";

        int w1 = static_cast<int>(header1.size());
        int w2 = static_cast<int>(header2.size());

        for (size_t var_idx = 0; var_idx < static_field.size(); var_idx++) {
            w1 = std::max(w1, static_cast<int>(fmt::format("{}[{}]", static_name, var_idx + 1).size()));
            w2 = std::max(w2, static_cast<int>(fmt::format("{:.6e}", static_field[var_idx]).size()));
        }

        TableFormat tfp({w1, w2}, {Align::Right, Align::Right});

        Log::row(tfp, {header1, header2});
        Log::dashes(tfp);
        for (size_t var_idx = 0; var_idx < static_field.size(); var_idx++) {
            Log::row(tfp, {fmt::format("{}[{}]", static_name, var_idx + 1), fmt::format("{:.6e}", static_field[var_idx])});
        }
        Log::dashes_ln(tfp);
    }
}

// === Shared I/O Helpers ===

/**
 * @brief write CSV with schema header.
 */
int write_csv(
    const std::string& filename,
    const std::vector<f64>& t,
    const std::vector<DynamicField>& dynamics,
    const std::vector<StaticField>& statics,
    bool write_header)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        Log::error("Could not open file {}", filename);
        abort();
    }

    file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);

    if (write_header) {
        // --- sections / header ---
        size_t col_index = 0;
        file << "# time: " << col_index << "\n";
        col_index++;

        for (const auto& dyn : dynamics) {
            size_t dim = dyn.data.empty() ? 0 : dyn.data.size();
            file << "# dynamic: \"" << dyn.name << "\", [";
            for (size_t var_idx = 0; var_idx < dim; var_idx++) {
                file << (col_index + var_idx);
                if (var_idx + 1 < dim) file << ",";
            }
            file << "]\n";
            col_index += dim;
        }

        for (const auto& st : statics) {
            size_t dim = st.data.size();
            file << "# static: \"" << st.name << "\", [";
            for (size_t var_idx = 0; var_idx < dim; var_idx++) {
                file << (col_index + var_idx);
                if (var_idx + 1 < dim) file << ",";
            }
            file << "]\n";
            col_index += dim;
        }
    }

    // --- column names ---
    file << "time";
    for (const auto& dyn : dynamics) {
        size_t dim = dyn.data.empty() ? 0 : dyn.data.size();
        for (size_t var_idx = 0; var_idx < dim; var_idx++) file << "," << dyn.name << "_" << var_idx;
    }
    for (const auto& st : statics) {
        for (size_t var_idx = 0; var_idx < st.data.size(); var_idx++) file << "," << st.name << "_" << var_idx;
    }
    file << "\n";

    // --- data ---
    for (size_t t_idx = 0; t_idx < t.size(); t_idx++) {
        file << t[t_idx];

        for (const auto& dyn : dynamics) {
            for (size_t var_idx = 0; var_idx < dyn.data.size(); var_idx++) {
                file << "," << dyn.data[var_idx][t_idx];
            }
        }

        if (t_idx == 0 || !write_header) {
            for (const auto& st : statics) {
                for (f64 val : st.data) file << "," << val;
            }
        }

        file << "\n";
    }

    return 0;
}

std::vector<size_t> parse_schema_indices(const std::string& colblock_str) {
    std::vector<size_t> indices;
    std::string colblock = colblock_str;
    colblock.erase(std::remove(colblock.begin(), colblock.end(), ' '), colblock.end());
    size_t l = colblock.find('[');
    size_t r = colblock.find(']');
    if (l != std::string::npos && r != std::string::npos) {
        std::string inside = colblock.substr(l + 1, r - l - 1);
        std::stringstream strsi(inside);
        std::string idx_str;
        while (std::getline(strsi, idx_str, ',')) {
            try {
                indices.push_back(std::stoul(idx_str));
            } catch (const std::exception& e) {
                Log::error("Invalid index in schema: {}", idx_str.c_str());
                return {};
            }
        }
    }
    return indices;
}

void read_csv(
    const std::string& filename,
    std::vector<f64>& t,
    std::map<std::string, std::vector<std::vector<f64>>>& dynamics,
    std::map<std::string, std::vector<f64>>& statics)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        Log::error("Could not open file {}", filename);
        abort();
    }

    std::string line;
    std::map<std::string, std::vector<size_t>> schema_dyn_indices;
    std::map<std::string, std::vector<size_t>> schema_stat_indices;
    size_t schema_time_index = MAX_SIZE;
    bool has_schema = false;
    
    // --- parse schema headers if present ---
    while (std::getline(file, line) && line[0] == '#') {
        has_schema = true;
        std::stringstream strs(line.substr(1));
        std::string tag, field_name;

        strs >> tag;
        
        if (tag == "time:") {
            strs >> schema_time_index;
        } else {
            std::string temp_name;
            strs >> temp_name;

            size_t name_start = line.find('"') + 1;
            size_t name_end = line.find('"', name_start);
            field_name = line.substr(name_start, name_end - name_start);

            size_t indices_start = line.find('[', name_end);
            size_t indices_end = line.find(']', indices_start);
            std::string indices_str = line.substr(indices_start, indices_end - indices_start + 1);
            
            if (tag == "dynamic:") {
                schema_dyn_indices[field_name] = parse_schema_indices(indices_str);
            } else if (tag == "static:") {
                schema_stat_indices[field_name] = parse_schema_indices(indices_str);
            }
        }
    }
    
    if (line.empty() && !file.eof()) {
        Log::error("File contains schema but no header or data rows.");
        return;
    }
    
    std::map<std::string, std::vector<size_t>> final_dyn_indices;
    std::map<std::string, std::vector<size_t>> final_stat_indices;
    size_t final_time_index = std::numeric_limits<size_t>::max();

    if (has_schema) {
        final_time_index = schema_time_index;
        final_dyn_indices = schema_dyn_indices;
        final_stat_indices = schema_stat_indices;
    } else {
        Log::warning("No schema header found. Defaulting to dynamic fields based on column names.");
        std::stringstream strs_header(line);
        std::string cell;
        size_t col_idx = 0;
        while (std::getline(strs_header, cell, ',')) {
            cell.erase(std::remove_if(cell.begin(), cell.end(), isspace), cell.end());
            if (cell == "time") {
                final_time_index = col_idx;
            } else if (!cell.empty()) {
                size_t underscore = cell.rfind('_');
                if (underscore != std::string::npos) {
                    std::string base_name = cell.substr(0, underscore);
                    final_dyn_indices[base_name].push_back(col_idx);
                }
            }
            col_idx++;
        }
    }

    bool stat_read = false;

    for (auto& kv : final_dyn_indices) {
        dynamics[kv.first].resize(kv.second.size());
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream strs_data(line);
        std::vector<std::string> cells;
        std::string cell;
        while (std::getline(strs_data, cell, ',')) cells.push_back(cell);

        if (cells.empty()) continue;

        if (final_time_index != std::numeric_limits<size_t>::max() && final_time_index < cells.size() && !cells[final_time_index].empty()) {
            try {
                t.push_back(std::stod(cells[final_time_index]));
            } catch (const std::exception& e) {
                t.push_back(std::numeric_limits<f64>::quiet_NaN());
            }
        } else {
            t.push_back(std::numeric_limits<f64>::quiet_NaN());
        }

        for (auto& kv : final_dyn_indices) {
            size_t var_idx = 0;
            for (size_t idx : kv.second) {
                auto& row_data = dynamics[kv.first][var_idx];

                if (idx < cells.size() && !cells[idx].empty()) {
                    try {
                        row_data.push_back(std::stod(cells[idx]));
                    } catch (const std::exception& e) {
                        row_data.push_back(std::numeric_limits<f64>::quiet_NaN());
                    }
                } else {
                    row_data.push_back(std::numeric_limits<f64>::quiet_NaN());
                }
                var_idx++;
            }
        }

        if (has_schema && !stat_read) {
            for (auto& kv : final_stat_indices) {
                for (size_t idx : kv.second) {
                    if (idx < cells.size() && !cells[idx].empty()) {
                        try {
                            statics[kv.first].push_back(std::stod(cells[idx]));
                        } catch (const std::exception& e) {
                            statics[kv.first].push_back(std::numeric_limits<f64>::quiet_NaN());
                        }
                    } else {
                         statics[kv.first].push_back(std::numeric_limits<f64>::quiet_NaN());
                    }
                }
            }
            stat_read = true;
        }
    }
}
