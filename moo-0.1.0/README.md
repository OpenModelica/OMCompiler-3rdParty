[![Build-Linux-x86_64](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-linux.yml/badge.svg)](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-linux.yml)
[![Build-Linux-Arm](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-linux-arm.yml/badge.svg)](https://github.com/AMIT-HSBI/MOO/actions/workflows/build-linux-arm.yml)
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

Install with your favorit package manager

- metis

  - Debian / Ubuntu `apt install metis`

- LAPACK

  - Debian / Ubuntu `apt install libblas-dev liblapack-dev gfortran`

### Configure

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
```

Possible configuration arguments:

- `MOO_METIS_LIB`: Location of metis
- `MOO_LAPACK_LIB`: Location of LAPACK
- `MOO_GFORTRAN_LIB`: Location of GNU Fortran

### Build

```bash
cmake --build build --parallel <Nr. of cores> --target all
```

### Test

Add `-DMOO_TESTS=ON` to the CMake configuration step.
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
