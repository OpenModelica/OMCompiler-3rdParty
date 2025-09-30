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

#include "fLGR.h"

/**
 * @brief Approximate the integral of a function over [0, 1] using collocation weights. @runtime: O(scheme).
 * 
 * @note  Provides exact integration for all polynomials of degree <= 2 * scheme - 2
 * 
 * This computes:
 * \f[
 * \int_0^1 f(t)\,dt \approx \sum_{k=1}^{m} b_k \cdot f(c_k)
 * \f]
 * where:
 * - \( f(c_k) \) are the input values at collocation nodes,
 * - \( b_k \) are the collocation weights for integration.
 *
 * @param scheme  Degree of the collocation scheme (number of collocation points).
 * @param values  Array of function values at the collocation nodes \( f(c_1), \dots, f(c_m) \), of length `scheme`.
 * @return        Approximation of the integral over the interval [0, 1].
 */
f64 fLGR::integrate(int scheme, const f64* values) {
    return Linalg::dot(scheme, values, get_b(scheme));
};

/**
 * @brief Apply the full collocation differentiation matrix D to an input vector. @runtime: O(scheme^2).
 *
 * This computes y := D * x where D is the (scheme + 1) × (scheme + 1) collocation differentiation matrix.
 *
 * @param scheme  Degree of the collocation scheme (number of collocation points in the interval).
 * @param in      Input vector of length (scheme + 1), e.g., values at t_0, ..., t_p.
 * @param out     Output vector of length (scheme + 1), contains the result of D * in.
 */
void fLGR::diff_matrix_multiply(int scheme, const f64* in, f64* out) {
    Linalg::matrix_vector(scheme + 1, 'R', get_D(scheme), in, out);
}

/**
 * @brief Applies the collocation differentiation matrix to a block of state vectors. @runtime: O(scheme^2 * x_size).
 *
 * This function computes one block of collocation constraints using the differentiation matrix
 * for a given collocation scheme. It performs a strided matrix-vector multiplication, where
 * both the input and output are flattened with strides to accommodate interleaved variables.
 *
 * where `x_new` is accessed with a stride (`x_stride`) to handle flattened and interleaved
 * state/control vectors (e.g. x + u sizes), and `out` is written using a stride (`out_stride`) to account for
 * the full size of one collocation constraint block (e.g. f + g sizes).
 *
 * @param scheme     Degree of the collocation scheme (number of collocation points in the interval).
 * @param x_size     Number of state variables (size of each x vector).
 * @param x_stride   Stride between consecutive vectors in `x_new` (i.e. size of x + u).
 * @param out_stride Stride between output rows (i.e. size of one full constraint block, size of f + g).
 * @param x_prev     Pointer to the previous node's state/control values (at t_i == t_{i-1, m_{i-1}}).
 * @param x_new      Pointer to strided new values at collocation points (from t_{i, 0} to t_{i, m_{i}}),
 * @param out        Pointer to the output buffer where results are accumulated (has stride out_stride).
 */
void fLGR::diff_matrix_multiply_block_strided(int scheme, int x_size, int x_stride, int out_stride,
                                              const f64* x_prev, const f64* x_new, f64* out) {
    for (int row = 1; row < scheme + 1; row++) {
        int out_row_base = (row - 1) * out_stride;

        for (int x_index = 0; x_index < x_size; x_index++) {
            int out_row = out_row_base + x_index;

            // col = 0 -> x_pref; D_{:, 0} * x_0
            out[out_row] += get_D(scheme, row, 0) * x_prev[x_index];

            // col > 0 -> x_new; D_{:, col} * x_col, col > 0
            for (int col = 1; col < scheme + 1; col++) {
                out[out_row] += get_D(scheme, row, col) * x_new[(col - 1) * x_stride + x_index];
            }
        }
    }
};

/**
 * @brief Interpolates a value at a given point T using barycentric interpolation. @runtime: O(scheme).
 *
 * This method performs barycentric interpolation based on a set of collocation nodes
 * and corresponding barycentric weights. It first rescales the target point T
 * to the nominal interval domain of the nodes and then computes the interpolated
 * value. It includes a check for exact matches with existing nodes to prevent
 * division by zero.
 *
 * @param scheme         The number of collocation points (nodes and weights).
 * @param contains_zero  A boolean indicating whether the `c0`/`w0` (containing zero)
 *                       or `c`/`w` (not containing zero) sets of nodes/weights should be used.
 * @param values         A pointer to an array of function values corresponding to the collocation nodes.
 * @param increment      The stride to use when accessing elements in the `values` array.
 * @param interval_start The start of the physical interval. (time of values[0], first)
 * @param interval_end   The end of the physical interval.   (time of values[-1], last)
 * @param point          The physical time at which to interpolate.
 * @return The interpolated value at point.
 */
f64 fLGR::interpolate(int scheme, bool contains_zero, const f64* values, int increment,
                      f64 interval_start, f64 interval_end, f64 point) {
    const auto& nodes   = contains_zero ? get_c0(scheme) : get_c(scheme);
    const auto& weights = contains_zero ? get_w0(scheme) : get_w(scheme);
    int node_count = scheme + static_cast<int>(contains_zero);

    // rescale T to the [0, 1] nominal interval domain
    f64 h          = interval_end - interval_start;
    f64 node_start = nodes[0];
    f64 node_end   = nodes[node_count - 1];

    // check if interval is empty or if just one point is given -> return constant polynomial
    if ((std::abs(node_end - node_start) < 1e-14) || (!contains_zero && scheme == 1)) {
        return values[0];
    }

    f64 point_hat = (point - interval_start) / h * (node_end - node_start) + node_start;

    // check for exact match with any node to avoid division by zero
    for (int j = 0; j < node_count; j++) {
        if (std::abs(point_hat - nodes[j]) < 1e-14) return values[j];
    }

    // compute the barycentric interpolant
    f64 numerator = 0.0;
    f64 denominator = 0.0;
    for (int j = 0; j < node_count; j++) {
        f64 temp = weights[j] / (point_hat - nodes[j]);
        numerator += temp * values[increment * j];
        denominator += temp;
    }

    return numerator / denominator;
};
