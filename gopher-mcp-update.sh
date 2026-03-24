#!/bin/bash -ex

test $# = 1 # I merged a branch last time; hopefully main should work in the future

test -d gopher-mcp.git || git clone https://github.com/GopherSecurity/gopher-mcp.git gopher-mcp.git
(cd gopher-mcp.git && git fetch && git checkout $1 && git clean -fdx)

rm -rf gopher-mcp
mkdir gopher-mcp
cp -a gopher-mcp.git/* gopher-mcp

(cd gopher-mcp.git && cmake -B build)

(cd gopher-mcp.git/build/_deps && rm -rf *-src/.git *-src/tests *-src/test *-src/docs)

mkdir gopher-mcp/_deps

cp -a gopher-mcp.git/build/_deps/*-src gopher-mcp/_deps/

rm -rf gopher-mcp.git/build/

cat > gopher-mcp/FetchContent.cmake <<EOF
function(FetchContent_Declare)
endfunction()
function(FetchContent_MakeAvailable pkgName)
add_subdirectory(_deps/\${pkgName}-src)
endfunction()
EOF

sed -i 's/include(FetchContent)/include(FetchContent.cmake)/' gopher-mcp/CMakeLists.txt
sed -i 's/add_subdirectory.tests./#add_subdirectory(tests)/' gopher-mcp/_deps/nghttp2-src/CMakeLists.txt

git add gopher-mcp
