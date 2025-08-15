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
#include "mesh.h"

/** 
 * @brief Create a mesh with fixed (h, p) for final time tf
 *
 * @param intervals  Number of Intervals
 * @param tf         Final Time
 * @param stages     Number of Nodes for each Interval
 * @return Mesh      Mesh
 */
std::shared_ptr<const Mesh> Mesh::create_equidistant_fixed_stages(f64 tf, int intervals, int stages) {
    FixedVector<f64> grid(intervals + 1);
    FixedVector<f64> delta_t(intervals);
    FixedVector<int> nodes(intervals);
    FixedField<int, 2> acc_nodes(intervals, stages);
    FixedField<f64, 2> t(intervals, stages);

    f64 h = tf / intervals;
    for (int i = 0; i < intervals; i++) {
        grid[i] = i * h;
    }
    grid[intervals] = tf;

    for (int i = 0; i < intervals; i++) {
        delta_t[i] = h;
        nodes[i] = stages;
        for (int j = 0; j < stages; j++) {
            acc_nodes[i][j] = stages * i + j;
            t[i][j] = grid[i] + delta_t[i] * fLGR::get_c(stages, j);
        }
    }
    int node_count = stages * intervals;
    return std::shared_ptr<const Mesh>(
        new Mesh(
            intervals,
            tf,
            std::move(grid),
            std::move(delta_t),
            std::move(t),
            std::move(nodes),
            std::move(acc_nodes),
            node_count
        )
    );
}

std::shared_ptr<const Mesh> Mesh::create_from_mesh_update(std::unique_ptr<MeshUpdate> mesh_update) {
    FixedVector<f64> grid = std::move(mesh_update->new_grid);
    FixedVector<int> nodes = std::move(mesh_update->new_nodes_per_interval);

    int intervals = nodes.int_size();
    f64 tf = grid.back();

    FixedVector<f64> delta_t(intervals);
    for (int i = 0; i < intervals; i++) {
        delta_t[i] = grid[i + 1] - grid[i];
    }

    int node_count = std::accumulate(nodes.begin(), nodes.end(), 0);

    FixedField<f64, 2> t(intervals);
    FixedField<int, 2> acc_nodes(intervals);

    int global_index = 0;
    for (int i = 0; i < intervals; ++i) {
        int p = nodes[i];
        f64 h = delta_t[i];

        t[i] = FixedVector<f64>(p);
        acc_nodes[i] = FixedVector<int>(p);

        for (int j = 0; j < p; ++j) {
            t[i][j] = grid[i] + h * fLGR::get_c(p, j);
            acc_nodes[i][j] = global_index++;
        }
    }

    return std::shared_ptr<const Mesh>(
        new Mesh(
            intervals,
            tf,
            std::move(grid),
            std::move(delta_t),
            std::move(t),
            std::move(nodes),
            std::move(acc_nodes),
            node_count
        )
    );
}

std::vector<f64> Mesh::get_flat_t() const {
    std::vector<f64> flat_t;

    flat_t.reserve(node_count + 1);

    flat_t.push_back(0.0);

    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < nodes[i]; j++) {
            flat_t.push_back(t[i][j]);
        }
    }

    return flat_t;
}
