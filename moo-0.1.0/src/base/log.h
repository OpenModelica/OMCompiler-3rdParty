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

// TODO: add proper logging with SPD log, override Ipopt Journalist -> useful for 3rd party logging; add colors as runtime option?

#ifndef MOO_LOG_H
#define MOO_LOG_H

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/color.h>
#include <string>
#include <array>

// ---------------------------
// General Log Macros
// ---------------------------

#define LOG(...) \
    do { fmt::print(__VA_ARGS__); fmt::print("\n"); } while (0)

#define LOG_T(tabs, ...) \
    do { fmt::print("{0:{1}}", "", (tabs) * 4); fmt::print(__VA_ARGS__); } while (0)

#define LOG_PREFIX(c, ...) \
    do { fmt::print("{} ", c); fmt::print(__VA_ARGS__); } while (0)

#define LOG_SUCCESS(...) do { fmt::print("SUCCESS - "); fmt::println(__VA_ARGS__); } while (0)
#define LOG_WARNING(...) do { fmt::print("WARNING - "); fmt::println(__VA_ARGS__); } while (0)
#define LOG_ERROR(...) do { fmt::print("ERROR - "); fmt::println(__VA_ARGS__); } while (0)

#define LOG_START_MODULE(fmt, title) \
do { \
    log_start_module(title, (fmt).total_width()); \
} while (0)

#define LOG_ROW(ftf, ...)       do { log_row_fast(ftf, __VA_ARGS__); }         while (0)
#define LOG_HEADER(ftf, ...)    do { log_row_fast(ftf, __VA_ARGS__); }         while (0)
#define LOG_DASHES(ftf)         do { log_dashes_fast(ftf); }                   while (0)
#define LOG_DASHES_LN(ftf)      do { log_dashes_fast(ftf); fmt::println(""); } while (0)

enum class Align { Left, Right, Center };

constexpr const char* align_to_char(Align a) {
    switch (a) {
        case Align::Left:   return "<";
        case Align::Right:  return ">";
        case Align::Center: return "^";
    }
    return "^";  // fallback center
}

template <size_t N>
struct FixedTableFormat {
    std::array<int, N> col_widths;
    std::array<std::string, N> fmt_strings;

    constexpr FixedTableFormat(const std::array<int, N>& widths,
                                const std::array<Align, N>& aligns)
        : col_widths(widths) {
        for (size_t i = 0; i < N; ++i) {
            fmt_strings[i] = fmt::format("{{:{}{}}}", align_to_char(aligns[i]), col_widths[i]);
        }
    }

    int total_width() const {
        int total = -1;
        for (auto w : col_widths)
            total += w + 3;
        return total;
    }
};

inline void log_start_module(const std::string& title, int table_width) {
    const std::string base_sep = "===";
    int base_sep_len = 3;
    int min_pad = 1;

    int total_title_len = base_sep_len + min_pad + static_cast<int>(title.size()) + min_pad + base_sep_len;
    int width = std::max(total_title_len, table_width);
    int padding = width - total_title_len;

    int left_eqs  = padding / 2;
    int right_eqs = padding - left_eqs;

    fmt::print("\n");
    fmt::print("{:=<{}}", base_sep, base_sep_len + left_eqs);
    fmt::print(" {} ", title);
    fmt::print("{:=>{}}\n\n", base_sep, base_sep_len + right_eqs);
}

template <size_t N>
constexpr inline void log_dashes_fast(const FixedTableFormat<N>& ftf) {
    fmt::print("{:-<{}}\n", "", ftf.total_width());
}

template <size_t N, typename... Args>
constexpr inline void log_row_fast(const FixedTableFormat<N>& ftf, Args&&... args) {
    static_assert(sizeof...(Args) == N, "Number of columns must match format definition.");
    size_t i = 0;
    const char* sep = "";
    ((fmt::print("{}{}", sep, fmt::format(ftf.fmt_strings[i++], args)), sep = " | "), ...);
    fmt::print("\n");
}

#endif // MOO_LOG_H
