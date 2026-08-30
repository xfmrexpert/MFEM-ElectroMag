cmake_minimum_required(VERSION 3.14) # FetchContent requires 3.14+
project(MFEM-ElectroMag VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Build options
option(BUILD_TESTING "Build tests" ON)
option(ENABLE_NATIVE_ARCH "Enable -march=native for Release builds" OFF)
option(USE_OPENMP "Enable OpenMP for parallel assembly" ON)

# Compiler flags for different build types
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        add_compile_options(/O2 /DNDEBUG)
        if(ENABLE_NATIVE_ARCH)
            add_compile_options(/arch:AVX2)
        endif()
    else()
        if(ENABLE_NATIVE_ARCH)
            add_compile_options(-O3 -DNDEBUG -march=native)
        else()
            add_compile_options(-O3 -DNDEBUG)
        endif()
    endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(MSVC)
        add_compile_options(/Od /Zi)
        # MSVC supports AddressSanitizer (/fsanitize=address) but not UBSan.
        option(ENABLE_SANITIZERS "Enable address sanitizer" OFF)
        if(ENABLE_SANITIZERS)
            add_compile_options(/fsanitize=address)
        endif()
    else()
        add_compile_options(-O0 -g)
        # Enable sanitizers for Debug builds (optional, can be disabled)
        option(ENABLE_SANITIZERS "Enable address and undefined behavior sanitizers" ON)
        if(ENABLE_SANITIZERS)
            add_compile_options(-fsanitize=address -fsanitize=undefined)
            add_link_options(-fsanitize=address -fsanitize=undefined)
        endif()
    endif()
endif()

include(FetchContent)

# Make MFEM throw mfem::ErrorException on failed MFEM_ASSERT/MFEM_VERIFY
# instead of calling std::abort(). This lets us catch malformed-mesh errors
# raised from inside the mfem::Mesh constructor (e.g. CheckBdrElementOrientation
# indexing faces_info with -1 for orphan boundary edges) and report them
# uniformly. Must be set before FetchContent_MakeAvailable(mfem).
set(MFEM_USE_EXCEPTIONS ON CACHE BOOL "Throw mfem::ErrorException on MFEM errors" FORCE)

# Enable MFEM's internal debug checks (MFEM_ASSERT and the extra MFEM_VERIFY
# blocks guarded by #ifdef MFEM_DEBUG, e.g. the interior-face orientation
# check in Mesh::Finalize). Independent of CMAKE_BUILD_TYPE. Override with
# -DMFEM_DEBUG=OFF on the cmake command line if you need to match a release
# build's looser behavior. Must be set before FetchContent_MakeAvailable(mfem).
# set(MFEM_DEBUG ON CACHE BOOL "Enable MFEM internal debug checks" FORCE)

# Automatically download and build MFEM
FetchContent_Declare(
  mfem
  GIT_REPOSITORY https://github.com/mfem/mfem.git
  GIT_TAG        v4.9  # Pin to a stable version for reproducibility
  GIT_PROGRESS   TRUE  # Stream clone progress to the CMake output pane so
                       # configure does not appear to hang after
                       # "Detecting CXX compile features".
  EXCLUDE_FROM_ALL  # Don't build MFEM's miniapps/examples as part of "Build All"
)

# This creates the 'mfem' target available for linking.
# EXCLUDE_FROM_ALL above keeps MFEM out of the default build set so Visual
# Studio's Output window isn't dominated by MFEM compile messages; the 'mfem'
# library is still built on demand when mfem-electromag links it.
FetchContent_MakeAvailable(mfem)

# Eigen supplies the sparse Cholesky (SimplicialLDLT) behind the "direct"
# linear solver option. Eigen is header-only, so SOURCE_SUBDIR points at a
# directory that does not exist: that skips Eigen's own (heavy) CMake project
# entirely and just downloads the headers, which we add as an include dir.
# Only the MPL2-licensed subset is used; EIGEN_MPL2_ONLY makes any accidental
# use of an LGPL-licensed Eigen module a compile error, keeping this project
# safely MIT-licensed.
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        3.4.0  # Pin to a stable version for reproducibility
  GIT_PROGRESS   TRUE   # See MFEM block above for rationale.
  SOURCE_SUBDIR  do-not-configure-eigen
)
FetchContent_MakeAvailable(eigen)

add_library(eigen_headers INTERFACE)
target_include_directories(eigen_headers SYSTEM INTERFACE ${eigen_SOURCE_DIR})
target_compile_definitions(eigen_headers INTERFACE EIGEN_MPL2_ONLY)

# nlohmann/json supplies the configuration parser. Fetched rather than vendored
# so its version is pinned in one visible place, consistent with MFEM and Eigen.
# JSON_DIAGNOSTICS/JSON_DIAGNOSTIC_POSITIONS include the JSON pointer path and
# byte positions in exception messages so config errors (e.g. type_error.302)
# point at the offending field/line. They are set here as INTERFACE definitions
# on the target: the ABI tag depends on them, so carrying them with the target
# keeps every linked TU consistent automatically.
set(JSON_BuildTests OFF CACHE INTERNAL "")
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.12.0  # Matches the previously vendored src/json.hpp
  GIT_PROGRESS   TRUE     # See MFEM block above for rationale.
)
FetchContent_MakeAvailable(nlohmann_json)
target_compile_definitions(nlohmann_json
    INTERFACE JSON_DIAGNOSTICS=1 JSON_DIAGNOSTIC_POSITIONS=1)

add_executable(mfem-electromag src/main.cpp)
target_link_libraries(mfem-electromag PRIVATE mfem eigen_headers nlohmann_json::nlohmann_json)

# Project headers are included by bare name ("physics_solver.hpp"), so that
# neither src/ nor test/ depends on the relative path between them.
target_include_directories(mfem-electromag PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Build identity banner: regenerate build_info.hpp on every build so a stale
# executable is obvious at runtime. GenerateBuildInfo.cmake uses configure_file,
# which only rewrites the header when the commit/dirty state actually changes,
# so this does not trigger a rebuild on every invocation.
add_custom_target(build_info ALL
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_VERSION=${PROJECT_VERSION}
        -DGIT_DIR=${CMAKE_SOURCE_DIR}
        -DSRC=${CMAKE_SOURCE_DIR}/src/io/build_info.hpp.in
        -DDST=${CMAKE_BINARY_DIR}/generated/build_info.hpp
        -P ${CMAKE_SOURCE_DIR}/cmake/GenerateBuildInfo.cmake
    BYPRODUCTS ${CMAKE_BINARY_DIR}/generated/build_info.hpp
    COMMENT "Generating build info header"
    VERBATIM
)
add_dependencies(mfem-electromag build_info)
target_include_directories(mfem-electromag PRIVATE ${CMAKE_BINARY_DIR}/generated)

# OpenMP support
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(mfem-electromag PRIVATE OpenMP::OpenMP_CXX)
        message(STATUS "OpenMP found and enabled. Use OMP_NUM_THREADS to control parallelism.")
    else()
        message(WARNING "OpenMP requested but not found. Building without OpenMP support.")
    endif()
endif()

# Documentation with Doxygen
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
    message(STATUS "Doxygen found. Use 'make docs' to generate documentation.")
else()
    message(STATUS "Doxygen not found. Documentation generation will be skipped.")
endif()

# Testing
if(BUILD_TESTING)
    enable_testing()

    # Fetch Catch2
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.5.2  # Pin to stable version
        GIT_PROGRESS   TRUE    # See MFEM block above for rationale.
    )
    FetchContent_MakeAvailable(Catch2)

    # Create test executable
    add_executable(mfem_tests
        test/test_main.cpp
        test/test_axisymmetric_axis.cpp
        test/test_axisymmetric_quadrature.cpp
        test/test_config_validator.cpp
        test/test_input_parser.cpp
        test/test_mqs_port_operator.cpp
        test/test_solvers.cpp
        test/test_status_reporter.cpp
    )

    target_link_libraries(mfem_tests PRIVATE mfem eigen_headers nlohmann_json::nlohmann_json Catch2::Catch2WithMain)
    target_include_directories(mfem_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)

    # Register tests with CTest
    list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
    include(CTest)
    include(Catch)
    catch_discover_tests(mfem_tests)
endif()

