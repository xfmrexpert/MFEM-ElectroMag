cmake_minimum_required(VERSION 3.14) # FetchContent requires 3.14+
project(InductiveSolver CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# Automatically download and build MFEM
FetchContent_Declare(
  mfem
  GIT_REPOSITORY https://github.com/mfem/mfem.git
  GIT_TAG        v4.7  # Pin to a stable version for reproducibility
)

# This creates the 'mfem' target available for linking
FetchContent_MakeAvailable(mfem)

add_executable(solve_axi src/main.cpp)
target_link_libraries(solve_axi PRIVATE mfem)
