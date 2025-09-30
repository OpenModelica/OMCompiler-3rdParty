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

#ifndef MOO_UTIL_H
#define MOO_UTIL_H

#include <functional>
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <any>

/* simple typedef for the Number, for using f32 or something later */
typedef double f64;

/* max f64 == _DBL_MAX_  */
const f64 PLUS_INFINITY = std::numeric_limits<f64>::max();

/* min f64 == -_DBL_MAX_ */
const f64 MINUS_INFINITY = -std::numeric_limits<f64>::max();

/* max size_t */
const size_t MAX_SIZE = std::numeric_limits<size_t>::max();

template <typename T>
inline int int_size(const std::vector<T>& vec) {
    return static_cast<int>(vec.size());
}

template <typename T>
inline T sign(T value) { return (value > 0 ? 1.0 : -1.0); }

template <typename T>
inline T apply_threshold_floor(T value, T tol, T min_magnitude) {
    return (std::abs(value) < tol) ? sign(value) * min_magnitude : value;
}

#endif // MOO_UTIL_H
