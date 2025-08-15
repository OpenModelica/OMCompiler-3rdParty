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

#ifndef MOO_MESH_H
#define MOO_MESH_H

#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>

#include "fixed_vector.h"
#include "fLGR.h"
#include "util.h"

struct MeshUpdate {
    FixedVector<f64> new_grid;
    FixedVector<int> new_nodes_per_interval;

    MeshUpdate(FixedVector<f64>&& new_grid,
               FixedVector<int>&& new_nodes_per_interval)
    : new_grid(std::move(new_grid)),
      new_nodes_per_interval(std::move(new_nodes_per_interval)) {}
};

class Mesh : public std::enable_shared_from_this<Mesh> {
public:
    int intervals;                // number of intervals
    int node_count;               // number of collocation nodes (sum over all intervals)
    f64 tf;                       // final time
    FixedVector<f64>   grid;      // grid base points
    FixedVector<f64>   delta_t;   // step size h for each interval
    FixedField<f64, 2> t;         // mesh points t_{i, j} = grid[i] + delta_t[i] * c[i][j] (c: collocation nodes for interval i)
    FixedVector<int>   nodes;     // number of collocation nodes p for each interval
    FixedField<int, 2> acc_nodes; // number of nodes to the left of index (i, j)

    static std::shared_ptr<const Mesh> create_equidistant_fixed_stages(f64 tf, int intervals, int p);

    static std::shared_ptr<const Mesh> create_from_mesh_update(std::unique_ptr<MeshUpdate> mesh_update);

    std::vector<f64> get_flat_t() const;


private:
    Mesh(int intervals,
         f64 tf,
         FixedVector<f64>&& grid,
         FixedVector<f64>&& delta_t,
         FixedField<f64, 2>&& t,
         FixedVector<int>&& nodes,
         FixedField<int, 2>&& acc_nodes,
         int node_count)
         : intervals(intervals),
           node_count(node_count),
           tf(tf),
           grid(std::move(grid)),
           delta_t(std::move(delta_t)),
           t(std::move(t)),
           nodes(std::move(nodes)),
           acc_nodes(std::move(acc_nodes)) {
    }
};

#endif  // MOO_MESH_H
