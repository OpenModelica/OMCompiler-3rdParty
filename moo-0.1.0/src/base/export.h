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

#ifndef MOO_EXPORT_H
#define MOO_EXPORT_H

#if defined(_WIN32) || defined(__CYGWIN__)
  #if defined(MOO_DLL_EXPORT)
    #define MOO_EXPORT __declspec(dllexport)
  #else
    #define MOO_EXPORT __declspec(dllimport)
  #endif
#else
  #if __GNUC__ >= 4
    #define MOO_EXPORT __attribute__((visibility("default")))
  #else
    #define MOO_EXPORT
  #endif
#endif

#endif // MOO_EXPORT_H
