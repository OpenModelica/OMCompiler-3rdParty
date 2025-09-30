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

#include "log.h"

void StdoutLogger::log(LogLevel lvl, std::string msg)
{
    switch (lvl) {
        case LogLevel::Info:
            fmt::print("{}\n", msg);
            break;
        case LogLevel::Success:
            fmt::print("SUCCESS - {}\n", msg);
            break;
        case LogLevel::Warning:
            fmt::print("WARNING - {}\n", msg);
            break;
        case LogLevel::Error:
            fmt::print("ERROR - {}\n", msg);
            break;
    }
}

namespace Log {

static std::unique_ptr<Logger> _globl_logger = std::make_unique<StdoutLogger>();

Logger* global_logger()
{
    return _globl_logger.get();
}

void set_global_logger(std::unique_ptr<Logger>&& logger)
{
    _globl_logger = std::move(logger);
}

} // namespace Log
