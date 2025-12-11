[![Build-Linux-x86](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-linux.yml/badge.svg)](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-linux.yml)
[![Build-Linux-ARM](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-linux-arm.yml/badge.svg)](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-linux-arm.yml)
[![Build-Windows](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-windows.yml/badge.svg)](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-windows.yml)
[![Build-macOS](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-macos.yml/badge.svg)](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-macos.yml)

# **MOO: Modelica / Model Optimizer - A Generic Framework for Dynamic Optimization**

This is **MOO: Modelica / Model Optimizer v.0.1.0**, a flexible and extensible
framework for solving optimization problems. **MOO** provides a generic
Nonlinear Programming (NLP) layer with built-in scaling support and a generic
NLP solver interface. While primarily designed for dynamic optimization in the
Modelica ecosystem (General Dynamic Optimization Problems - GDOPs, training of
Physics-enhanced Neural ODEs - PeN-ODEs), it is equally applicable to other
domains, e.g. model predictive control (MPC).

## Working with this repository

Clone with submodules:

```bash
git clone --recurse-submodules git@github.com:AMIT-HSBI/MOO.git
```

## Compilation

**MOO** uses CMake to compile.

### Dependencies

Install with your favorite package manager

#### Necessary
- LAPACK

  - Debian / Ubuntu: `apt install liblapack-dev`
  - Latest OpenBLAS build from source: add `-DUSE_SYSTEM_LAPACK=OFF -DDOWNLOAD_LAPACK=ON` to the CMake configure command

- gfortran (TODO: check if Flang build works, is this necessary?)

  - Debian / Ubuntu: `apt install gfortran`

#### Optional
- METIS

  - Debian / Ubuntu `apt install libmetis-dev`
  - if METIS is not available: add `-DMUMPS_HAS_METIS=OFF` to the CMake configure command

- HSL

  - add `-DIPOPT_HAS_HSL` to the CMake configure command (package should be found via PkgConfig)

### Configure

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
```

Possible configuration arguments:

- `MOO_WITH_RADAU`: Build with RADAU fortran code to perform simulations (default: `ON`)
- `MOO_WITH_GDOPT`: Build with GDOPT interface (default: `ON`, WIP)
- `MOO_WITH_C_INTERFACE`: Build MOO with the generic C interface (default: `ON`, WIP)

### Build

```bash
cmake --build build --parallel <Nr. of cores> --target all
```

### Test

Add `-DMOO_TESTS=ON` to the CMake configure command.
After building run the tests:

```bash
cmake --build build --target test
```

### Development

Use [act](https://github.com/nektos/act) to test the GitHub workflows locally:

```bash
act -P ubuntu-latest=catthehacker/ubuntu:full-latest --artifact-server-path $PWD/.artifacts
```

## License

The **Modelica/Model Optimizer (MOO)** is distributed under the **GNU Lesser
General Public License (LGPL) Version 3**. See the full license text for
detailed terms and conditions:
[LGPL-3.0](https://www.gnu.org/licenses/lgpl-3.0.html)


# Features and Functionality
## General Dynamic Optimization Problem (GDOP)


### Free Initial or Final Time
The most generic problem directly implemented in MOO is a General Dynamic Optimization Problem (GDOP) with free control variables $u(t)$, free static parameters
$p$, free final time $t_0$ and free final time $t_f$ of the form:

$$
\begin{aligned}
\min_{u(t), p, t_0, t_f}\quad & M(x_0, u_0, x_f, u_f, p, t_0, t_f)
+ \int_{t_0}^{t_f} L\bigl(x(t), u(t), p\bigr)\, dt \\[6pt]
\text{s.t.}\quad
& \frac{dx}{dt} = f\bigl(x(t), u(t), p\bigr), \quad t \in [t_0, t_f], \\[4pt]
& g^L \le g\bigl(x(t), u(t), p\bigr) \le g^U, \quad t \in [t_0, t_f], \\[4pt]
& r^L \le r\bigl(x_0, u_0, x_f, u_f, p, t_0, t_f\bigr) \le r^U, \\[4pt]
& x^L \le x(t) \le x^U, \quad
  u^L \le u(t) \le u^U, \quad
  t \in [t_0, t_f] \\
& p^L \le p \le p^U, \quad
  t_0^L \le t_0 \le t_0^U, \quad
  t_f^L \le t_f \le t_f^U.
\end{aligned}
$$

### Fixed Initial or Final Time

If the problem is given on a fixed time horizon $[t_0, t_f]$, the library can also solve problems of the form:

$$
\begin{aligned}
\min_{u(t), p}\quad & M(x_0, u_0, x_f, u_f, p)
+ \int_{t_0}^{t_f} L\bigl(x(t), u(t), p, t\bigr)\, dt \\[6pt]
\text{s.t.}\quad
& \frac{dx}{dt} = f\bigl(x(t), u(t), p, t\bigr), \quad t \in [t_0, t_f], \\[4pt]
& g^L \le g\bigl(x(t), u(t), p, t\bigr) \le g^U, \quad t \in [t_0, t_f], \\[4pt]
& r^L \le r\bigl(x_0, u_0, x_f, u_f, p\bigr) \le r^U, \\[4pt]
& x^L \le x(t) \le x^U, \quad
  u^L \le u(t) \le u^U, \quad
  t \in [t_0, t_f] \\
& p^L \le p \le p^U. \quad
\end{aligned}
$$

In this version, the user is also allowed to use the provided time variable $t$ in the functions $L, f, g$.
