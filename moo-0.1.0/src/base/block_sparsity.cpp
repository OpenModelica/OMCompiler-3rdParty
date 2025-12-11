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

#include <set>
#include <memory>
#include <stdexcept>

#include <base/fixed_vector.h>
#include <base/nlp_structs.h>
#include <base/export.h>
#include <base/log.h>
#include <base/block_sparsity.h>


/* creates a dense lower triangular (with diagonal) block structure
    0 | 1 2 . s
    -----------
    1 | x      
    2 | x x    
    . | x x x  
    s | x x x x
*/
BlockSparsity BlockSparsity::create_lower_triangular(const int size, const BlockType block_type) {
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
BlockSparsity BlockSparsity::create_rectangular(const int rows, const int cols, const BlockType block_type) {
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

// mapping (row, col) -> index in some larger sparsity structure
int BlockSparsity::access(const int row, const int col) const {
    switch (type) {
        case BlockType::Exact:
            return block[row][col];
        default:
            Log::error("Unknown BlockType in BlockSparsity::access().");
            abort();
    }
}

// mapping (row, col, block_count) -> index in some larger sparsity structure
int BlockSparsity::access(const int row, const int col, const int block_count) const {
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

void BlockSparsity::print() const {
    for (auto const& v : block) {
        v.print();
    }
}


/* creates a **fully** dense rectangular block structure Rows x Cols
    0 | 1 2 . C
    -----------
    1 | x x x x
    . | x x x x
    R | x x x x
*/
DenseRectangularBlockSparsity DenseRectangularBlockSparsity::create(const int rows, const int cols) {
    DenseRectangularBlockSparsity b;
    b.row_offset_prev = FixedVector<int>(rows);
    b.row_size = cols;
    b.nnz = rows * cols;
    return b;
}

bool OrderedIndexSet::Compare::operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
    if (a.first != b.first) {
        return a.first < b.first;
    } else {
        return a.second < b.second;
    }
}

// standard Hessian insertion
void OrderedIndexSet::insert_sparsity(const std::vector<HessianSparsity>& hes, int row_off, int col_off) {
    for (const auto& coo : hes) {
        assert((!is_diag_block || coo.row + row_off >= coo.col + col_off) && "Hessian must be lower triangular!");
        set.insert({coo.row + row_off, coo.col + col_off});
    }
}

// insertion for Jacobian, e.g. if because of product rule Jacobian terms must be included (see GDOP Block K)
//                              these stem from the D * x - deltaT * f(x, u, p) -> partial w.r.t. tf / t0 and p
void OrderedIndexSet::insert_sparsity(std::vector<int>& rows, const std::vector<JacobianSparsity>& jac, int row_off, int col_off) {
    for (const auto row : rows) {
        for (const auto& jac_elem : jac) {
            assert((!is_diag_block || row + row_off >= jac_elem.col + col_off) && "Hessian must be lower triangular!");
            set.insert({row + row_off, jac_elem.col + col_off});
        }
    }
}
