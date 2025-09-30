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
#include <memory>
#include <stdexcept>

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
    static BlockSparsity create_lower_triangular(const int size, const BlockType block_type) {
        BlockSparsity b;
        b.type = block_type;
        b.block = FixedField<int, 2>(size);

        for (int i = 0; i < size; i++) {
            b.block[i] = FixedVector<int>(i + 1);
        }
        switch (block_type) {
            case BlockType::Offset:
                b.off_prev = 0;
                break;
            case BlockType::RowOffset:
                b.row_offset_prev = FixedVector<int>(size);
                b.row_size = FixedVector<int>(size);
                break;
            default:
                break;
        }
        return b;
    }

    /* creates a dense rectangular block structure Rows x Cols
       0 | 1 2 . C
       -----------
       1 | x x x x
       . | x x x x
       R | x x x x
    */
    static BlockSparsity create_rectangular(const int rows, const int cols, const BlockType block_type) {
        BlockSparsity b;
        b.type = block_type;
        b.block = FixedField<int, 2>(rows, cols);

        switch (block_type) {
            case BlockType::Offset:
                b.off_prev = 0;
                break;
            case BlockType::RowOffset:
                b.row_offset_prev = FixedVector<int>(rows);
                b.row_size = FixedVector<int>(rows);
                break;
            default:
                break;
        }
        return b;
    }

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
    int access(const int row, const int col) const {
        switch (type) {
            case BlockType::Exact:
                return block[row][col];
            default:
                Log::error("Unknown BlockType in BlockSparsity::access().");
                abort();
        }
    }

    // mapping (row, col, block_count) -> index in some larger sparsity structure
    int access(const int row, const int col, const int block_count) const {
        switch (type) {
            // for A -> B: row, col, offset_prev - #nnz at end of blocktype (e.g. |A|), block_count (e.g. (i,j) = (0, 2) => 2)
            case BlockType::Offset:
                return off_prev + nnz * block_count + block[row][col];

            // for E -> F
            case BlockType::RowOffset:
                return row_offset_prev[row] + row_size[row] * block_count + block[row][col];

            default:
               return access(row, col);
        }
    }

    void print() const {
        for (auto const& v : block) {
            v.print();
        }
    }
};

struct MOO_EXPORT OrderedIndexSet {
    struct Compare {
        bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
            if (a.first != b.first) {
                return a.first < b.first;
            } else {
                return a.second < b.second;
            }
        }
    };

    std::set<std::pair<int, int>, Compare> set;

    void insert_sparsity(const std::vector<HessianSparsity>& hes, int row_off, int col_off) {
        for (auto coo : hes) {
            assert(coo.row + row_off >= coo.col + col_off); // Hessian must be lower triangular!
            set.insert({coo.row + row_off, coo.col + col_off});
        }
    }

    inline int size() const {
        return set.size();
    }

    void clear() {
        set.clear();
    }
};

#endif // MOO_BLOCK_SPARSITY_H
