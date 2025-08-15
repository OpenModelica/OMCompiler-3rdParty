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

const f64 PLUS_INFINITY = std::numeric_limits<f64>::infinity();
const f64 MINUS_INFINITY = -std::numeric_limits<f64>::infinity();

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

/* AutoFree manages a list of raw pointers and their corresponding free functions.
 * When an AutoFree object goes out of scope, it automatically calls each stored free function
 * on its associated pointer to release resources and avoid memory leaks.
 * You register pointers and their free functions using the attach() method.
 * Used for pure C-style mallocs / callocs in the MOO module */

/* @deprecated not in use
 class AutoFree {
public:
    ~AutoFree() {
        for (auto& [ptr, free_fn] : _to_free) {
            free_fn(ptr);
        }
    }

    template <typename T>
    void attach(T* ptr, void (*free_fn)(T*)) {
        if (ptr != nullptr)
            _to_free.emplace_back(ptr, [free_fn](void* p) { free_fn(static_cast<T*>(p)); });
    }

    template <typename T>
    void attach(std::initializer_list<T*> ptrs, void (*free_fn)(T*)) {
        for (T* ptr : ptrs) {
            if (ptr != nullptr)
                attach(ptr, free_fn);
        }
    }

private:
    std::vector<std::pair<void*, std::function<void(void*)>>> _to_free;
};
*/
#endif // MOO_UTIL_H
