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

#ifndef MOO_BLOCK_SPARSITY_H
#define MOO_BLOCK_SPARSITY_H

#include <set>

#include <base/fixed_vector.h>
#include <base/nlp_structs.h>
#include <base/export.h>
#include <base/log.h>


enum class BlockType {
    Exact,
    Offset,
    RowOffset,
};

struct MOO_EXPORT BlockSparsity {
    // use traditional C-like struct with type, because the polymorphism is very simple
    BlockType type;

    // data for all block types
    FixedField<int, 2> block;
    int nnz = 0;

    // for row-based offsets (e.g. block F in GDOP)
    FixedVector<int> row_offset_prev;
    FixedVector<int> row_size;

    // for constant offsets (e.g. block B in GDOP)
    int off_prev = 0;

    /* creates a dense lower triangular (with diagonal) block structure
       0 | 1 2 . s
       -----------
       1 | x      
       2 | x x    
       . | x x x  
       s | x x x x
    */
    static BlockSparsity create_lower_triangular(const int size, const BlockType block_type);

    /* creates a dense rectangular block structure Rows x Cols
       0 | 1 2 . C
       -----------
       1 | x x x x
       . | x x x x
       R | x x x x
    */
    static BlockSparsity create_rectangular(const int rows, const int cols, const BlockType block_type);

    /* creates a dense square block structure Size x Size
       0 | 1 . C
       -----------
       1 | x x x
       . | x x x
       R | x x x
    */
    inline static BlockSparsity create_square(const int size, const BlockType block_type) {
        return create_rectangular(size, size, block_type);
    }

    inline void insert(const int row, const int col, const int index) {
        block[row][col] = index;
        nnz++;
    }

    // mapping (row, col) -> index in some larger sparsity structure
    int access(const int row, const int col) const;

    // mapping (row, col, block_count) -> index in some larger sparsity structure
    int access(const int row, const int col, const int block_count) const;

    void print() const;
};

struct MOO_EXPORT DenseRectangularBlockSparsity {
    // for row-based offsets (e.g. block F in GDOP)
    FixedVector<int> row_offset_prev;
    int row_size; // == #cols here
    int nnz;

    /* creates a **fully** dense rectangular block structure Rows x Cols
       0 | 1 2 . C
       -----------
       1 | x x x x
       . | x x x x
       R | x x x x
    */
    static DenseRectangularBlockSparsity create(const int rows, const int cols);

    // mapping (row, col) -> index in some larger sparsity structure
    inline int access(const int row, const int col, const int block_count) const {
        return row_offset_prev[row] + row_size * block_count + col;
    }
};

struct MOO_EXPORT OrderedIndexSet {
    struct Compare {
        bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const;
    };

    std::set<std::pair<int, int>, Compare> set;

    // set to true if the block is on the diagonal (e.g. xu + xu, but not p with xu) => enforces row >= col
    bool is_diag_block = false;

    // standard Hessian insertion
    void insert_sparsity(const std::vector<HessianSparsity>& hes, int row_off, int col_off);

    // insertion for Jacobian, e.g. if because of product rule Jacobian terms must be included (see GDOP Block K)
    //                              these stem from the D * x - deltaT * f(x, u, p) -> partial w.r.t. tf / t0 and p
    void insert_sparsity(std::vector<int>& rows, const std::vector<JacobianSparsity>& jac, int row_off, int col_off);

    inline int size() const {
        return set.size();
    }

    inline void clear(bool is_diag) {
        is_diag_block = is_diag;
        set.clear();
    }
};

#endif // MOO_BLOCK_SPARSITY_H
