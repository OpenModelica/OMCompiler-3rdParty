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

Mesh::Mesh(int intervals,
           f64 t0,
           f64 tf,
           FixedVector<f64>&& grid,
           FixedVector<f64>&& delta_t,
           FixedField<f64, 2>&& t,
           FixedVector<int>&& nodes,
           FixedField<int, 2>&& acc_nodes,
           int node_count)
    : intervals(intervals),
      node_count(node_count),
      t0(t0),
      tf(tf),
      grid(std::move(grid)),
      delta_t(std::move(delta_t)),
      t(std::move(t)),
      nodes(std::move(nodes)),
      acc_nodes(std::move(acc_nodes)) {}

std::shared_ptr<Mesh> Mesh::create_same_type(
    int intervals,
    f64 t0,
    f64 tf,
    FixedVector<f64>&& grid,
    FixedVector<f64>&& delta_t,
    FixedField<f64, 2>&& t,
    FixedVector<int>&& nodes,
    FixedField<int, 2>&& acc_nodes,
    int node_count) const
{
    return std::shared_ptr<Mesh>(
        new Mesh(
            intervals,
            t0,
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

/** 
 * @brief Create a mesh with fixed (h, p) for final time tf
 *
 * @param intervals  Number of Intervals
 * @param tf         Final Time
 * @param stages     Number of Nodes for each Interval
 * @return Mesh      Mesh
 */
std::shared_ptr<Mesh> Mesh::create_equidistant_fixed_stages(f64 t0, f64 tf, int intervals, int stages, MeshType creation_type) {
    FixedVector<f64> grid(intervals + 1);
    FixedVector<f64> delta_t(intervals);
    FixedVector<int> nodes(intervals);
    FixedField<int, 2> acc_nodes(intervals, stages);
    FixedField<f64, 2> t(intervals, stages);

    f64 h = (tf - t0) / intervals;
    for (int i = 0; i < intervals; i++) {
        grid[i] = t0 + i * h;
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

    switch (creation_type) {
        case MeshType::Spectral:
        {
            FixedVector<f64> spectral_grid(grid.size());
            f64 invT = 1.0 / (tf - t0);

            for (size_t i = 0; i < grid.size(); i++) {
                spectral_grid[i] = invT * (grid[i] - t0);
            }

            FixedVector<f64> spectral_delta_t(intervals);
            for (int i = 0; i < intervals; i++) {
                spectral_delta_t[i] = spectral_grid[i + 1] - spectral_grid[i];
            }

            return std::shared_ptr<SpectralMesh>(
                new SpectralMesh(
                    intervals,
                    t0,
                    tf,
                    std::move(grid),
                    std::move(delta_t),
                    std::move(t),
                    std::move(nodes),
                    std::move(acc_nodes),
                    node_count,
                    std::move(spectral_grid),
                    std::move(spectral_delta_t)
                )
            );
        }
        case MeshType::Physical:
        default:
            return std::shared_ptr<Mesh>(
                new Mesh(
                    intervals,
                    t0,
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
}

std::shared_ptr<Mesh> Mesh::create_from_mesh_update(std::unique_ptr<MeshUpdate> mesh_update) const {
    FixedVector<f64> new_grid = std::move(mesh_update->new_grid);
    FixedVector<int> new_nodes = std::move(mesh_update->new_nodes_per_interval);

    int new_intervals = new_nodes.int_size();
    f64 new_t0 = new_grid[0];
    f64 new_tf = new_grid.back();

    FixedVector<f64> new_delta_t(new_intervals);
    for (int i = 0; i < new_intervals; i++) {
        new_delta_t[i] = new_grid[i + 1] - new_grid[i];
    }

    int new_node_count = std::accumulate(new_nodes.begin(), new_nodes.end(), 0);

    FixedField<f64, 2> new_t(new_intervals);
    FixedField<int, 2> new_acc_nodes(new_intervals);

    int global_index = 0;
    for (int i = 0; i < new_intervals; i++) {
        int p = new_nodes[i];
        f64 h = new_delta_t[i];

        new_t[i] = FixedVector<f64>(p);
        new_acc_nodes[i] = FixedVector<int>(p);

        for (int j = 0; j < p; j++) {
            new_t[i][j] = new_grid[i] + h * fLGR::get_c(p, j);
            new_acc_nodes[i][j] = global_index++;
        }
    }

    // let polymorphism handle the type
    return create_same_type(
        new_intervals,
        new_t0,
        new_tf,
        std::move(new_grid),
        std::move(new_delta_t),
        std::move(new_t),
        std::move(new_nodes),
        std::move(new_acc_nodes),
        new_node_count
    );
}

std::vector<f64> Mesh::get_flat_t() const {
    std::vector<f64> flat_t;

    flat_t.reserve(node_count + 1);
    flat_t.push_back(t0);

    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < nodes[i]; j++) {
            flat_t.push_back(t[i][j]);
        }
    }

    return flat_t;
}

SpectralMesh::SpectralMesh(int intervals,
                           f64 t0,
                           f64 tf,
                           FixedVector<f64>&& grid,
                           FixedVector<f64>&& delta_t,
                           FixedField<f64, 2>&& t,
                           FixedVector<int>&& nodes,
                           FixedField<int, 2>&& acc_nodes,
                           int node_count,
                           FixedVector<f64>&& spectral_grid,
                           FixedVector<f64>&& spectral_delta_t)
    : Mesh(intervals,
           t0,
           tf,
           std::move(grid),
           std::move(delta_t),
           std::move(t),
           std::move(nodes),
           std::move(acc_nodes),
           node_count),
      spectral_grid(std::move(spectral_grid)),
      spectral_delta_t(std::move(spectral_delta_t)) {}

// virtual (pseudo static) SpectralMesh factory
// is non static for override with derived SpectralMesh
std::shared_ptr<Mesh> SpectralMesh::create_same_type(
    int intervals,
    f64 t0,
    f64 tf,
    FixedVector<f64>&& grid,
    FixedVector<f64>&& delta_t,
    FixedField<f64, 2>&& t,
    FixedVector<int>&& nodes,
    FixedField<int, 2>&& acc_nodes,
    int node_count) const
{
    assert(intervals == grid.int_size() - 1);

    FixedVector<f64> new_spectral_grid(grid.size());
    f64 invT = 1.0 / (tf - t0);

    for (size_t i = 0; i < grid.size(); i++) {
        new_spectral_grid[i] = invT * (grid[i] - t0);
    }

    FixedVector<f64> new_spectral_delta_t(intervals);
    for (int i = 0; i < intervals; i++) {
        new_spectral_delta_t[i] = new_spectral_grid[i + 1] - new_spectral_grid[i];
    }

    assert(new_spectral_delta_t.size() == new_spectral_grid.size() - 1);

    return std::shared_ptr<SpectralMesh>(
        new SpectralMesh(
            intervals,
            t0,
            tf,
            std::move(grid),
            std::move(delta_t),
            std::move(t),
            std::move(nodes),
            std::move(acc_nodes),
            node_count,
            std::move(new_spectral_grid),
            std::move(new_spectral_delta_t)
        )
    );
}

// updates t0, tf, grid, delta_t and t from base Mesh object
void SpectralMesh::update_physical_from_spectral(f64 new_t0, f64 new_tf) {
    assert(grid.size() == spectral_grid.size());
    assert(delta_t.size() == spectral_delta_t.size());

    if (new_t0 == t0 && new_tf == tf) {
        return;
    }

    t0 = new_t0;
    tf = new_tf;

    f64 len = tf - t0;

    for (int i = 0; i < grid.int_size(); i++) {
        grid[i] = t0 + len * spectral_grid[i];
    }

    for (int i = 0; i < intervals; i++) {
        delta_t[i] = grid[i + 1] - grid[i];
    }

    for (int i = 0; i < intervals; i++) {
        int p = nodes[i];
        for (int j = 0; j < p; j++) {
            t[i][j] = grid[i] + delta_t[i] * fLGR::get_c(p, j);
        }
    }
}
