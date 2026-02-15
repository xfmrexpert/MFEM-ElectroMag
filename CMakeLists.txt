cmake_minimum_required(VERSION 3.14) # FetchContent requires 3.14+
project(MFEM-ElectroMag CXX)

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
    if(ENABLE_NATIVE_ARCH)
        add_compile_options(-O3 -DNDEBUG -march=native)
    else()
        add_compile_options(-O3 -DNDEBUG)
    endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0 -g)
    # Enable sanitizers for Debug builds (optional, can be disabled)
    option(ENABLE_SANITIZERS "Enable address and undefined behavior sanitizers" ON)
    if(ENABLE_SANITIZERS)
        add_compile_options(-fsanitize=address -fsanitize=undefined)
        add_link_options(-fsanitize=address -fsanitize=undefined)
    endif()
endif()

include(FetchContent)

# Automatically download and build MFEM
FetchContent_Declare(
  mfem
  GIT_REPOSITORY https://github.com/mfem/mfem.git
  GIT_TAG        v4.9  # Pin to a stable version for reproducibility
)

# This creates the 'mfem' target available for linking
FetchContent_MakeAvailable(mfem)

add_executable(mfem-electromag src/main.cpp)
target_link_libraries(mfem-electromag PRIVATE mfem)

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
    )
    FetchContent_MakeAvailable(Catch2)

    # Create test executable
    add_executable(mfem_tests
        test/test_main.cpp
        test/test_input_parser.cpp
        test/test_solvers.cpp
    )

    target_link_libraries(mfem_tests PRIVATE mfem Catch2::Catch2WithMain)

    # Register tests with CTest
    list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
    include(CTest)
    include(Catch)
    catch_discover_tests(mfem_tests)
endif()

