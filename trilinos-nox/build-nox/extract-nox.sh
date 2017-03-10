# Use this script to build a tarball out of the Trilinos repository that only
# contains packages which are necessary to build the Nox solver.
# Therefor you must have the Trilinos, the OMCompiler/3rdParty and the
# OMCompiler repository on your system. Use it to update Nox in OpenModelica
# occasionally.
#
# The following variables have to be set by the user first.

# Set Trilinos tag or branch to use.
RELEASE_TAG="master"

# If 'RELEASE_TAG' is a real release tag named 'trilinos-release-*', you can
# leave the 'RELEASE_VERSION' empty. It will be set automatically. If you choose
# another tag or even a branch, the variable must be specified manually.
RELEASE_VERSION="12.11"


# Set path to Git if not already in environment.
GIT="C:/Program Files/Git/bin/git.exe"

# Set path to Trilinos repository.
TRILINOS_LIBRARY_DIR="C:/Users/MSD1LO/Documents/Library/Trilinos"

# Set path to 3rdParty repository.
THIRD_PARTY_DIR="C:/Users/MSD1LO/Documents/Stuff/3rdParty"

# Set path to BLAS library.
BLAS_LIBRARY_DIRS="C:/Users/MSD1LO/Documents/Library/OpenModelica/OMDev/lib/3rdParty/Lapack/Lib"

# Set path to LAPACK library.
LAPACK_LIBRARY_DIRS="C:/Users/MSD1LO/Documents/Library/OpenModelica/OMDev/lib/3rdParty/Lapack/Lib"

################################################################################

# Get version number.
if [[ -z "$RELEASE_VERSION" ]]; then
	RELEASE_VERSION="${RELEASE_TAG#trilinos-release-}"
	RELEASE_VERSION="${RELEASE_VERSION//-/.}"
fi
echo "Using Trilinos release version $RELEASE_VERSION ..."

# Look for Git.
if [[ -z "$GIT" ]]; then
	which git > /dev/null 2>&1
	[[ $? -eq 0 ]] && GIT="$(which git)" || (echo "Git not found" && exit 1)
fi
echo "Using Git $GIT ..."

# Checkout release tag in Trilinos repository.
"$GIT" -C "$TRILINOS_LIBRARY_DIR" checkout "$RELEASE_TAG"

# Create temporary directory.
TEMP_DIR="$(date +%s)"
mkdir -p "$TEMP_DIR"
cd "$TEMP_DIR"

# Prepare for tarball generation.
cmake \
	-G "MSYS Makefiles" \
	-DGIT_EXEC:FILEPATH="$GIT" \
	-DBLAS_LIBRARY_DIRS="$BLAS_LIBRARY_DIRS" \
	-DLAPACK_LIBRARY_DIRS="$LAPACK_LIBRARY_DIRS" \
	-DTrilinos_ENABLE_ALL_PACKAGES=OFF \
	-DTrilinos_ENABLE_EXAMPLES:BOOL=OFF \
	-DBUILD_SHARED_LIBS:BOOL=ON \
	-DTrilinos_ENABLE_NOX:BOOL=ON \
	-DTrilinos_ENABLE_TESTS:BOOL=OFF \
	-DTrilinos_ASSERT_MISSING_PACKAGES:BOOL=OFF \
	-DXpetra_ENABLE_Tpetra:BOOL=OFF \
	-DTrilinos_ENABLE_ThyraCore:BOOL=OFF \
	-DNOX_ENABLE_ThyraEpetraExtAdapters:BOOL=OFF \
	-DTrilinos_ENABLE_Gtest:BOOL=OFF \
	-DTrilinos_ENABLE_Epetra:BOOL=OFF \
	-DTrilinos_ENABLE_Galeri:BOOL=OFF \
	-DTrilinos_ENABLE_AztecOO:BOOL=OFF \
	-DTrilinos_ENABLE_Ifpack:BOOL=OFF \
	-DTrilinos_ENABLE_KokkosCore:BOOL=OFF \
	-DTrilinos_EXTRA_LINK_FLAGS:STRING="-lws2_32" \
	-DTrilinos_CPACK_SOURCE_GENERATOR:STRING="TGZ" \
	-DCMAKE_INSTALL_PREFIX:STRING="build" \
	"$TRILINOS_LIBRARY_DIR"

# Build tarball.
make package_source -j7

# Extract files and move them to 'OMCompiler/3rdParty'.
tar -xf "trilinos-$RELEASE_VERSION-serial-$RELEASE_VERSION-Source.tar.gz" \
    -C "$THIRD_PARTY_DIR/trilinos-nox" --strip-components=1

# Clean up.
rm "trilinos-$RELEASE_VERSION-serial-$RELEASE_VERSION-Source.tar.gz"
cd ..
rm -r "$TEMP_DIR"

# Checkout previous branch in Trilinos repository.
"$GIT" -C "$TRILINOS_LIBRARY_DIR" checkout -

# Add include to 'CMakeLists.txt'.
sed -i'' 's%INCLUDE(${CMAKE_SOURCE_DIR}/ProjectName.cmake)%INCLUDE(${CMAKE_SOURCE_DIR}/ProjectName.cmake)\n\n# Include CMake file for building Nox.\nINCLUDE(${CMAKE_SOURCE_DIR}/build-nox/NoxLists.cmake)%' "$THIRD_PARTY_DIR/trilinos-nox/CMakeLists.txt"

# Fix 'get_hostname' error on Windows.
sed -i'' 's#gethostname (hostname, 255);#throw std::runtime_error("Function *gethostname* not supported on Windows.");#' "$THIRD_PARTY_DIR/trilinos-nox/packages/teuchos/parameterlist/src/Teuchos_XMLPerfTestArchive.cpp"
sed -i'' 's#gethostname(hostname, sizeof(hostname));#pr_error("Function *gethostname* not supported on Windows.");#' "$THIRD_PARTY_DIR/trilinos-nox/packages/ml/src/Utils/ml_utils.c"
sed -i'' 's#sprintf(buf, "Host: %s   PID: %d (mpi task %d)", hostname, getpid(),mypid);##' "$THIRD_PARTY_DIR/trilinos-nox/packages/ml/src/Utils/ml_utils.c"

# Fix 'expected unqualified-id' error.
sed -i'' 's#template BLAS<long int, std::complex<float> >;##' "$THIRD_PARTY_DIR/trilinos-nox/packages/teuchos/numerics/src/Teuchos_BLAS.cpp"
sed -i'' 's#template BLAS<long int, std::complex<double> >;##' "$THIRD_PARTY_DIR/trilinos-nox/packages/teuchos/numerics/src/Teuchos_BLAS.cpp"
sed -i'' 's#template BLAS<long int, float>;##' "$THIRD_PARTY_DIR/trilinos-nox/packages/teuchos/numerics/src/Teuchos_BLAS.cpp"
sed -i'' 's#template BLAS<long int, double>;##' "$THIRD_PARTY_DIR/trilinos-nox/packages/teuchos/numerics/src/Teuchos_BLAS.cpp"
