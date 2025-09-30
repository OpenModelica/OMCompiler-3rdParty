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

#ifndef MOO_LOG_H
#define MOO_LOG_H

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/color.h>
#include <cassert>
#include <string>
#include <vector>
#include <array>
#include <memory>

#include <base/export.h>

// ---------------------------
// logger API (pluggable)
// ---------------------------

enum class LogLevel { Info, Success, Warning, Error };

class Logger {
public:
    virtual ~Logger() = default;

    // implementations should print a newline
    virtual void log(LogLevel lvl, std::string msg) = 0;
};

// default stdout logger
class StdoutLogger : public Logger {
public:
    void log(LogLevel lvl, std::string msg) override;
};

namespace Log {

MOO_EXPORT void set_global_logger(std::unique_ptr<Logger>&& logger);
MOO_EXPORT Logger* global_logger();

} // namespace Log

// ---------------------------
// formatting helpers (produce std::string)
// ---------------------------

struct LogFormatter {
    // overload: pass already formatted string
    static std::string info_string(const std::string& s) { return s; }
    static std::string success_string(const std::string& s) { return s; }
    static std::string warning_string(const std::string& s) { return s; }
    static std::string error_string(const std::string& s) { return s; }

    // format variants
    template <typename... Args>
    static std::string info_string(const char* fmtstr, Args&&... args) {
        return fmt::format(fmtstr, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static std::string success_string(const char* fmtstr, Args&&... args) {
        return fmt::format(fmtstr, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static std::string warning_string(const char* fmtstr, Args&&... args) {
        return fmt::format(fmtstr, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static std::string error_string(const char* fmtstr, Args&&... args) {
        return fmt::format(fmtstr, std::forward<Args>(args)...);
    }
};

// helpers for indentation/prefix
template <typename... Args>
static inline std::string format_with_indent(int tabs, const char* fmtstr, Args&&... args) {
    std::string msg = fmt::format(fmtstr, std::forward<Args>(args)...);
    if (tabs <= 0) return msg;
    return std::string(static_cast<size_t>(tabs) * 4, ' ') + msg;
}
static inline std::string format_with_indent(int tabs, const std::string& s) {
    if (tabs <= 0) return s;
    return std::string(static_cast<size_t>(tabs) * 4, ' ') + s;
}

template <typename... Args>
static inline std::string format_with_prefix(char c, const char* fmtstr, Args&&... args) {
    std::string msg = fmt::format(fmtstr, std::forward<Args>(args)...);
    return fmt::format("{} {}", c, msg);
}
static inline std::string format_with_prefix(char c, const std::string& s) {
    return fmt::format("{} {}", c, s);
}

// ---------------------------
// table formatting utilities
// ---------------------------

enum class Align { Left, Right, Center };

static constexpr const char* align_to_char(Align a) {
    switch (a) {
        case Align::Left:   return "<";
        case Align::Right:  return ">";
        case Align::Center: return "^";
    }
    return "^";  // fallback center
}

template <size_t N>
struct MOO_EXPORT FixedTableFormat {
    std::array<int, N> col_widths;
    std::array<std::string, N> fmt_strings;

    constexpr FixedTableFormat(const std::array<int, N>& widths,
                               const std::array<Align, N>& aligns)
        : col_widths(widths) {
        for (size_t i = 0; i < N; i++) {
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

// format-only helpers (return std::string)
static inline std::string format_start_module(const std::string& title, int table_width) {
    const std::string base_sep = "===";
    int base_sep_len = 3;
    int min_pad = 1;

    int total_title_len = base_sep_len + min_pad + static_cast<int>(title.size()) + min_pad + base_sep_len;
    int width = std::max(total_title_len, table_width);
    int padding = width - total_title_len;

    int left_eqs  = padding / 2;
    int right_eqs = padding - left_eqs;

    std::string msg;
    msg += "\n";
    msg += fmt::format("{:=<{}}", base_sep, base_sep_len + left_eqs);
    msg += fmt::format(" {} ", title);
    msg += fmt::format("{:=>{}}\n", base_sep, base_sep_len + right_eqs);
    return msg;
}

template <size_t N>
static inline std::string format_dashes(const FixedTableFormat<N>& ftf) {
    return fmt::format("{:-<{}}", "", ftf.total_width());
}

template <size_t N, typename... Args>
static inline std::string format_row(const FixedTableFormat<N>& ftf, Args&&... args) {
    static_assert(sizeof...(Args) == N, "Number of columns must match format definition.");
    size_t i = 0;
    const char* sep = "";
    std::string row;
    ((row += fmt::format("{}{}", sep, fmt::format(ftf.fmt_strings[i++], args)), sep = " | "), ...);
    return row;
}

struct TableFormat {
    std::vector<int> widths;
    std::vector<Align> aligns;

    TableFormat(std::initializer_list<int> w, std::initializer_list<Align> a)
        : widths(w.begin(), w.end()), aligns(a.begin(), a.end())
    {
        if (widths.size() != aligns.size()) {
            throw std::invalid_argument("Widths and aligns must have the same size");
        }
    }

    TableFormat(const std::vector<int>& w, const std::vector<Align>& a)
        : widths(w), aligns(a)
    {
        if (widths.size() != aligns.size()) {
            throw std::invalid_argument("Widths and aligns must have the same size");
        }
    }

    std::vector<std::string> fmt_strings() const {
        std::vector<std::string> fs;
        for (size_t i = 0; i < widths.size(); i++)
            fs.push_back(fmt::format("{{:{}{}}}", align_to_char(aligns[i]), widths[i]));
        return fs;
    }

    int total_width() const {
        int total = -1;
        for (auto w : widths) total += w + 3;
        return total;
    }
};

static inline std::string format_row(const TableFormat& tf, const std::vector<std::string>& cols) {
    if (cols.size() != tf.widths.size()) {
        return "[LOG ERROR] Column count does not match table format!";
    }

    auto fs = tf.fmt_strings();
    const char* sep = "";
    std::string row;
    for (size_t i = 0; i < cols.size(); i++) {
        row += fmt::format("{}{}", sep, fmt::format(fs[i], cols[i]));
        sep = " | ";
    }
    return row;
}

static inline std::string format_dashes(const TableFormat& tf) {
    return fmt::format("{:-<{}}", "", tf.total_width());
}

// ---------------------------
// general logging
// ---------------------------

namespace Log {

// string overloads (already-formatted string)
inline void info(const char* s)    { global_logger()->log(LogLevel::Info, std::string(s)); }
inline void success(const char* s) { global_logger()->log(LogLevel::Success, std::string(s)); }
inline void warning(const char* s) { global_logger()->log(LogLevel::Warning, std::string(s)); }
inline void error(const char* s)   { global_logger()->log(LogLevel::Error, std::string(s)); }

inline void info(const std::string& s)    { global_logger()->log(LogLevel::Info, s); }
inline void success(const std::string& s) { global_logger()->log(LogLevel::Success, s); }
inline void warning(const std::string& s) { global_logger()->log(LogLevel::Warning, s); }
inline void error(const std::string& s)   { global_logger()->log(LogLevel::Error, s); }

// format overloads
template <typename... Args>
inline void info(const char* fmtstr, Args&&... args) {
    auto s = LogFormatter::info_string(fmtstr, std::forward<Args>(args)...);
    global_logger()->log(LogLevel::Info, s);
}
template <typename... Args>
inline void success(const char* fmtstr, Args&&... args) {
    auto s = LogFormatter::success_string(fmtstr, std::forward<Args>(args)...);
    global_logger()->log(LogLevel::Success, s);
}
template <typename... Args>
inline void warning(const char* fmtstr, Args&&... args) {
    auto s = LogFormatter::warning_string(fmtstr, std::forward<Args>(args)...);
    global_logger()->log(LogLevel::Warning, s);
}
template <typename... Args>
inline void error(const char* fmtstr, Args&&... args) {
    auto s = LogFormatter::error_string(fmtstr, std::forward<Args>(args)...);
    global_logger()->log(LogLevel::Error, s);
}

// indent then log
template <typename... Args>
inline void info_t(int tabs, const char* fmtstr, Args&&... args) {
    auto s = format_with_indent(tabs, fmtstr, std::forward<Args>(args)...);
    global_logger()->log(LogLevel::Info, s);
}
inline void info_t(int tabs, const std::string& s) {
    auto s2 = format_with_indent(tabs, s);
    global_logger()->log(LogLevel::Info, s2);
}

// prefix char then log
template <typename... Args>
inline void prefixed(char c, const char* fmtstr, Args&&... args) {
    auto s = format_with_prefix(c, fmtstr, std::forward<Args>(args)...);
    global_logger()->log(LogLevel::Info, s);
}
inline void prefixed(char c, const std::string& s) {
    auto s2 = format_with_prefix(c, s);
    global_logger()->log(LogLevel::Info, s2);
}

// ---------------------------
// Table logging
// ---------------------------

template <typename TableType>
inline void start_module(const TableType& table, const std::string& title) {
    auto s = format_start_module(title, table.total_width());
    global_logger()->log(LogLevel::Info, s);
}

template <typename TableType>
inline void dashes(const TableType& table) {
    auto s = format_dashes(table);
    global_logger()->log(LogLevel::Info, s);
}

template <typename TableType>
inline void dashes_ln(const TableType& table) {
    dashes(table);
    global_logger()->log(LogLevel::Info, "");
}

template <size_t N, typename... Args>
inline void row(const FixedTableFormat<N>& ftf, Args&&... args) {
    auto s = format_row(ftf, std::forward<Args>(args)...);
    global_logger()->log(LogLevel::Info, s);
}

inline void row(const TableFormat& tf, const std::vector<std::string>& cols) {
    auto s = format_row(tf, cols);
    global_logger()->log(LogLevel::Info, s);
}

} // namespace Log

#endif // MOO_LOG_H
