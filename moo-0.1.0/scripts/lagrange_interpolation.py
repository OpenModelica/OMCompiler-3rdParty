# SPDX-License-Identifier: LGPL-3.0-or-later
#
# This file is part of MOO - Modelica / Model Optimizer
# Copyright (C) 2025 University of Applied Sciences and Arts
# Bielefeld, Faculty of Engineering and Mathematics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np

t_nodes = np.array([0.0, 0.0038762756430420551, 0.016123724356957945, 0.025])
u_nodes = np.array([0.26372968347885795, 0.26390015970602487, 0.26444040300798038, 0.26483346836268973])
t_query = 1.6777216e-05

def lagrange_interpolate(x_nodes, y_nodes, x_query):
    n = len(x_nodes)
    result = 0.0
    for i in range(n):
        term = y_nodes[i]
        for j in range(n):
            if i != j:
                term *= (x_query - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
        result += term
    return result

def linear_interpolate(t_nodes, u_nodes, t_query):
    for i in range(len(t_nodes) - 1):
        if t_nodes[i] <= t_query <= t_nodes[i + 1]:
            t1, t2 = t_nodes[i], t_nodes[i + 1]
            u1, u2 = u_nodes[i], u_nodes[i + 1]
            alpha = (t_query - t1) / (t2 - t1)
            return u1 + alpha * (u2 - u1)
    if t_query < t_nodes[0]:
        return u_nodes[0]
    else:
        return u_nodes[-1]

u_interp = lagrange_interpolate(t_nodes, u_nodes, t_query)
print(u_interp)

u_interp = linear_interpolate(t_nodes, u_nodes, t_query)
print(u_interp)
