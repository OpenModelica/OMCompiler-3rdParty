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

#ifndef MOO_FIXED_VECTOR_H
#define MOO_FIXED_VECTOR_H

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <memory>
#include <type_traits>
#include <iostream>

#include <base/util.h>
#include <base/export.h>

// helpers to distinguish iterator and range base field initialization
template<typename It, typename = void>
struct is_iterator : std::false_type {};

template<typename It>
struct is_iterator<It, std::void_t<
    typename std::iterator_traits<It>::iterator_category >> : std::true_type {};

template<typename It>
inline constexpr bool is_iterator_v = is_iterator<It>::value;

/**
 * @brief Templated, zero-initialized vector with fixed size.
 *
 * This class represents a fixed-size vector that is zero-initialized upon construction.
 * It provides basic functionality for accessing and manipulating elements, as well as
 * support for copying, nulling and moving vectors.
 *
 * @tparam T The type of elements stored in the vector.
 */
template<typename T>
class MOO_EXPORT FixedVector {
public:

    constexpr FixedVector() noexcept : _size{0}, _data{} {}

    // TODO: try using T[] instead of ptr<T[]>
    explicit FixedVector(std::size_t size) : _size{size}, _data{std::make_unique<T[]>(size)} {}

    constexpr FixedVector(const FixedVector &other) : FixedVector(other._size) {
        if constexpr (std::is_trivially_copyable_v<T>) {
            std::memcpy(_data.get(), other._data.get(), _size * sizeof(T));
        } else {
            for (std::size_t i = 0; i < _size; i++) {
                (*this)[i] = T(other[i]);
            }
        }
    }

    FixedVector(FixedVector &&other) noexcept = default;

    // assign based on iterator begin() and end(), guard only iterator
    template<typename It, std::enable_if_t<is_iterator_v<It>, int> = 0>
    constexpr FixedVector(It first, It last) : FixedVector(static_cast<std::size_t>(std::distance(first, last))) {
        assign(first, last);
    }

    // recursive constructor for nested FixedVector / FixedField
    template<typename... Args>
    FixedVector(std::size_t size, Args&&... args) : FixedVector(size) {
        for (std::size_t i = 0; i < size; i++) {
            _data[i] = T(std::forward<Args>(args)...);
        }
    }

    // access vector at index 0 <= index < vec.size()
    constexpr T& operator[](std::size_t index) {
        assert(index < _size); // this assert interacts very badly with the debugger (disable it when observing questionable behavior)

        return _data[index];
    }

    // access vector at index 0 <= index < vec.size()
    constexpr const T& operator[](std::size_t index) const {
        assert(index < _size);

        return _data[index];
    }

    // assign vector to other vector of equal size
    constexpr FixedVector& operator=(const FixedVector &other) {
        assert(_size == other._size);

        std::copy(other._data.get(), other._data.get() + _size, _data.get());

        return *this;
    }

    FixedVector& operator=(FixedVector&& other) noexcept = default;

    // fill entire vector with 0
    constexpr void fill_zero() {
        std::memset(_data.get(), 0, _size * sizeof(T));
    }

    // fills the vector with some data of the same len
    constexpr void assign(const T* data) {
        assert(data != nullptr);

        std::memcpy(_data.get(), data, _size * sizeof(T));
    }

    // fills the vector with some data with given len: vector[offset] = data[0], ..., vector[offset + len - 1] = data[len - 1]
    constexpr void assign(const T* data, std::size_t len, std::size_t offset = 0) {
        assert(data != nullptr);
        assert(len <= _size - offset);

        std::memcpy(_data.get() + offset, data, len * sizeof(T));
    }

    // assign the vector from a given iterator: vector[offset] = first, ..., vector[offset + len - 1] = last
    template<typename It>
    constexpr void assign(It first, It last, std::size_t offset = 0) {
        assert(static_cast<std::size_t>(std::distance(first, last)) <= _size - offset);

        std::copy(first, last, _data.get() + offset);
    }

    // write from vector -> data_buffer: data[i] = vector[offset + i], ..., data[len - 1] = vector[offset + len - 1],
    constexpr void write_to(T* data, std::size_t len, std::size_t offset = 0) const {
        assert(data != nullptr);
        assert(len <= _size - offset);

        std::memcpy(data, _data.get() + offset, len * sizeof(T));
    }

    // write from vector -> data_buffer: data[0] = vector[0], ..., data[_size - 1] = vector[_size - 1]
    constexpr void write_to(T* data) const {
        assert(data != nullptr);

        std::memcpy(data, _data.get(), _size * sizeof(T));
    }

    constexpr std::size_t size() const {
        return _size;
    }

    constexpr inline int int_size() const {
        return static_cast<int>(_size);
    }

    constexpr T& back() {
        assert(_size != 0);

        return _data[_size - 1];
    }

    constexpr const T& back() const {
        assert(_size != 0);

        return _data[_size - 1];
    }

    constexpr T* raw() {
        return _data.get();
    }

    constexpr const T* raw() const {
        return _data.get();
    }

    constexpr T* begin() noexcept {
        return _data.get();
    }

    constexpr T* end() noexcept {
        return _data.get() + _size;
    }

    constexpr const T* begin() const noexcept {
        return _data.get();
    }

    constexpr const T* end() const noexcept {
        return _data.get() + _size;
    }

    constexpr const bool empty() const noexcept {
        return (_size == 0);
    }

    void print() const {
        std::cout << "[";
        for (std::size_t i = 0; i < _size; i++) {
            std::cout << _data[i];
            if (i < _size - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
    }

private:
    std::size_t _size;
    std::unique_ptr<T[]> _data;
};

template<typename T, std::size_t Dim>
struct MOO_EXPORT FixedFieldRecursive {
    using type = FixedVector<typename FixedFieldRecursive<T, Dim - 1>::type>;
};

template<typename T>
struct MOO_EXPORT FixedFieldRecursive<T, 1> {
    using type = FixedVector<T>;
};

template<typename T, std::size_t Dim>
using FixedField = typename FixedFieldRecursive<T, Dim>::type;

#endif // MOO_FIXED_VECTOR_H
