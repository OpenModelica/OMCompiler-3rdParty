# Building latest Ipopt 3.14.19 + MUMPS 5.8.1 with CMake

Coin-OR Ipopt repository: [coin-or/Ipopt](https://github.com/coin-or/Ipopt.git)

This CMake version is forked from: [rjodon/coinor-ipopt-with-cmake](https://github.com/rjodon/coinor-ipopt-with-cmake.git)

This repository is still work in progress. You are welcome to open pull requests with fixes / updates or extensions, thank you!

---

## How to use

Create a build directory:
```bash
mkdir build && cd build
```

Run CMake configure (the default options should work out):
```bash
cmake .. <Options>
```

### Note:
You need a C, C++ and Fortran compiler on your system. Otherwise you can't build the necessary dependencies.

The default configuration requires **LAPACK** and **METIS** to be available on your system.
It will automatically build and install MUMPS and Ipopt.

However, you can also build LAPACK (latest OpenBLAS) from source by setting `-DDOWNLOAD_LAPACK=ON` and `-DUSE_SYSTEM_LAPACK=OFF`.

If the project cannot find your METIS installation, you can also provide the library `-DMETIS_LIB_PATH=path/to/metis.so` and include path `-DMETIS_INC_PATH=path/to/include` directly.

### Available CMake Options

You can generate the full table of all options with:
```bash
grep -R --include=CMakeLists.txt --include=\*.cmake "^[[:space:]]*option(" . \
  | sed -E 's/.*option[[:space:]]*\(([A-Za-z0-9_]+)[[:space:]]+"([^"]+)"[[:space:]]+([A-Z]+)\).*/|\1|\2|\3|/' \
  | column -t -s '|'
```

### Note:
Most options are untested for now!

### Build

Compile everything with
```bash
cmake --build . --parallel <Nr. of Cores>
```