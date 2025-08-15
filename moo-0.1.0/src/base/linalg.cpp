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
#include <cstring>

#include "linalg.h"

namespace Linalg {

/**
 * @brief Perform inplace squaring (x_1^2, x_2^2, ..., x_size^2)
 *
 * @param size      Vector length
 * @param x         Vector ptr
 */
void square(int size, f64* x) {
    for (int i = 0; i < size; i++) {
        x[i] *= x[i];
    }
}


/** 
 * @brief Perform dot product: ret = transpose(x) * y
 *
 * @param size      Vector length
 * @param x         First Vector ptr
 * @param y         Second Vector ptr
 * @return f64   transpose(vec1) * vec2
 */
f64 dot(int size, const f64* x, const f64* y) {
    f64 ret = 0.0;
    for (int i = 0; i < size; i++) {
        ret += x[i] * y[i];
    }
    return ret;
}

/**
 * @brief Matrix-vector multiply: out = matrix * vector
 *
 * @param size    Size of the matrix
 * @param format  'R' = row-major, 'C' = column-major
 * @param matrix  Pointer to flat matrix (size x size)
 * @param vector  Pointer to input vector (size)
 * @param out     Pointer to output vector (size)
 */
void matrix_vector(int size, char format, const double* matrix,
                   const double* vector, double* out)
{
    std::memset(out, 0, size * sizeof(f64));

    if (format == 'C' || format == 'c') {
        // Column-major
        for (int col = 0; col < size; col++) {
            for (int row = 0; row < size; row++) {
                out[row] += matrix[row + col * size] * vector[col];
            }
        }
    } else {
        // Row-major
        for (int row = 0; row < size; row++) {
            for (int col = 0; col < size; col++) {
                out[row] += matrix[row * size + col] * vector[col];
            }
        }
    }
}

void diagmat_vec(const f64* D, bool invD, const f64* x, int size, f64* out) {
    int i;

    if (invD) {
        for (i = 0; i < size; i++) {
            out[i] = x[i] / D[i];
        }
    }
    else {
        for (i = 0; i < size; i++) {
            out[i] = D[i] * x[i];
        }
    }
}

void diagmat_vec_inplace(const f64* D, bool invD, f64* x, int size) {
    int i;

    if (invD) {
        for (i = 0; i < size; i++) {
            x[i] /= D[i];
        }
    }
    else {
        for (i = 0; i < size; i++) {
            x[i] *= D[i];
        }
    }
}

/** 
 * @brief Perform diagonal scaled axpy: out = D * (x + beta * y) or out = D^(-1) * (x + beta * y)
 *
 * @param size      Vector length
 * @param x         Vector x
 * @param y         Vector y
 * @param D         Diagonal Matrix D
 * @param beta      Scaling factor of beta
 * @param invD      Invert Matrix D
 * @param out       Vector to fill
 */
void dsaxpy(int size, const f64* x, const f64* y, const f64* D, f64 beta, bool invD, f64* out) {
    int i;

    if (invD) {
        for (i = 0; i < size; i++) {
            out[i] = (x[i] + beta * y[i]) / D[i];
        }
    }
    else {
        for (i = 0; i < size; i++) {
            out[i] = D[i] * (x[i] + beta * y[i]);
        }
    }
}

/** 
 * @brief Perform diagonal general matrix vector: out = D * x + beta * y or D^(-1) * x + beta * y
 *
 * @param size      Vector length
 * @param x         Vector x
 * @param y         Vector y
 * @param D         Diagonal Matrix D
 * @param beta      Scaling factor of beta
 * @param invD      Invert Matrix D
 * @param out       Vector to fill
 */
void dgmv(int size, const f64* x, const f64* y, const f64* D, f64 beta, bool invD, f64* out) {
    int i;

    if (invD) {
        for (i = 0; i < size; i++) {
            out[i] = x[i] / D[i] + beta * y[i];
        }
    }
    else {
        for (i = 0; i < size; i++) {
            out[i] =  D[i] * x[i] + beta * y[i];
        }
    }
}

} // namespace Linalg
