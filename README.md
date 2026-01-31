# MFEM-ElectroMag

A finite element solver for electromagnetic problems using MFEM (Modular Finite Element Methods). This solver supports electrostatic, magnetostatic, and magnetoquasistatic problems in both axisymmetric and planar geometries.

## Features

- **Electrostatics**: Solves for electric potential and field distributions
- **Magnetostatics**: Solves for magnetic vector potential and field distributions
- **Magnetoquasistatics**: Time-harmonic eddy current problems
- **Axisymmetric and Planar**: Supports both 2D coordinate systems
- **JSON Configuration**: Easy problem setup via JSON files
- **ParaView Output**: Direct visualization of results
- **Comprehensive Testing**: Unit and integration tests with Catch2
- **Configuration Validation**: Automatic checking of input files
- **Solver Factory Pattern**: Extensible architecture for new physics types
- **OpenMP Support**: Parallel assembly for improved performance
- **API Documentation**: Doxygen-generated documentation

## Dependencies

### Required

- **CMake** (≥ 3.14): Build system
- **C++ Compiler**: Supporting C++17 standard
  - GCC (≥ 7.0)
  - Clang (≥ 5.0)
  - MSVC (≥ 2017)
- **MFEM** (v4.7): Automatically downloaded and built by CMake

### Optional

- **SuiteSparse**: For direct sparse solvers (recommended for better performance)
- **HYPRE**: For advanced preconditioners and solvers
- **METIS**: For mesh partitioning
- **OpenMP**: For parallel assembly (usually included with compiler)
- **Doxygen**: For generating API documentation
- **Catch2**: For running tests (automatically downloaded)

## Installation Instructions

### Linux (Ubuntu/Debian)

```bash
# Install dependencies
sudo apt-get update
sudo apt-get install -y cmake g++ git

# Optional: Install SuiteSparse for better solver performance
sudo apt-get install -y libsuitesparse-dev

# Clone the repository
git clone https://github.com/xfmrexpert/MFEM-ElectroMag.git
cd MFEM-ElectroMag

# Build
mkdir build
cd build
cmake ..
make -j$(nproc)
```

### macOS

```bash
# Install dependencies using Homebrew
brew install cmake git

# Optional: Install SuiteSparse
brew install suite-sparse

# Clone the repository
git clone https://github.com/xfmrexpert/MFEM-ElectroMag.git
cd MFEM-ElectroMag

# Build
mkdir build
cd build
cmake ..
make -j$(sysctl -n hw.ncpu)
```

### Windows

Using Visual Studio:

```powershell
# Install CMake and Git
# Download from: https://cmake.org/download/ and https://git-scm.com/

# Clone the repository
git clone https://github.com/xfmrexpert/MFEM-ElectroMag.git
cd MFEM-ElectroMag

# Create build directory
mkdir build
cd build

# Generate Visual Studio project
cmake ..

# Build using CMake
cmake --build . --config Release

# Or open the generated .sln file in Visual Studio and build
```

Using MinGW/MSYS2:

```bash
# Install dependencies
pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-gcc git

# Clone and build
git clone https://github.com/xfmrexpert/MFEM-ElectroMag.git
cd MFEM-ElectroMag
mkdir build
cd build
cmake .. -G "MinGW Makefiles"
mingw32-make -j
```

## Building with MFEM Options

MFEM will be automatically downloaded and configured. To enable optional MFEM features:

```bash
# In the build directory
cmake .. -DMFEM_USE_SUITESPARSE=ON -DMFEM_USE_METIS=ON
make -j
```

## Usage

After building, the executable `mfem-electromag` will be in the `build` directory:

```bash
# Run with a configuration file
./mfem-electromag path/to/config.json

# Example: Run test cases
./mfem-electromag ../test/electrostatic_test.json
./mfem-electromag ../test/magnetostatic_test.json
./mfem-electromag ../test/mqs_test.json

# Run with OpenMP parallelism (if enabled)
OMP_NUM_THREADS=4 ./mfem-electromag config.json
```

### Running Tests

```bash
# In build directory
ctest --output-on-failure

# Or run test executable directly
./mfem_tests
```

### Generating Documentation

```bash
# In build directory (requires Doxygen)
make docs

# Documentation will be in docs/doxygen/html/index.html
```

## Configuration

Problems are configured using JSON files. See the `test/` directory for examples of:
- `electrostatic_test.json`: Electrostatic problem setup
- `magnetostatic_test.json`: Magnetostatic problem setup
- `mqs_test.json`: Magnetoquasistatic problem setup

See `examples/` directory for complete example problems with documentation:
- `simple_capacitor/`: Parallel plate capacitor with analytical validation
- `solenoid/`: Magnetostatic coil with field calculations
- `eddy_current/`: Time-harmonic eddy current analysis

For detailed mathematical formulation, see `docs/math_formulation.md`.

## Output

Results are saved in ParaView format (`.pvd`, `.vtu` files) and can be visualized using:
- [ParaView](https://www.paraview.org/): Open-source visualization application
- [VisIt](https://visit-dav.github.io/visit-website/): Scientific data visualization tool

## Project Structure

```
.
├── src/                           # Source files
│   ├── main.cpp                   # Main entry point
│   ├── physics_solver.hpp         # Base solver class
│   ├── electrostatic_solver.hpp   # Electrostatic solver
│   ├── magnetostatic_solver.hpp   # Magnetostatic solver
│   ├── magnetoquasistatic_solver.hpp  # MQS solver
│   ├── input_parser.hpp           # JSON configuration parser
│   ├── solver_factory.hpp         # Solver factory pattern
│   ├── config_validator.hpp       # Configuration validation
│   ├── constants.hpp              # Physical constants
│   ├── enums.hpp                  # Type enumerations
│   └── axisymmetric_*.hpp         # Custom axisymmetric integrators
├── test/                          # Test configurations and unit tests
│   ├── test_*.cpp                 # Unit tests (Catch2)
│   └── *.json                     # Test configurations
├── examples/                      # Complete example problems
│   ├── simple_capacitor/          # Capacitor example
│   ├── solenoid/                  # Solenoid example
│   └── eddy_current/              # Eddy current example
├── docs/                          # Documentation
│   ├── math_formulation.md        # Mathematical formulation
│   ├── solver_performance.md      # Performance optimization guide
│   ├── future_3d.md               # 3D extension roadmap
│   └── future_amr.md              # AMR implementation guide
├── .github/workflows/             # CI/CD configuration
├── CMakeLists.txt                 # CMake build configuration
├── Doxyfile                       # Doxygen configuration
├── LICENSE                        # MIT License
└── README.md                      # This file
```

## Design Decisions

### Header-Only Implementation

This project uses a header-only design for the following reasons:

**Advantages:**
- **Simplicity**: No separate compilation of library components
- **Fast development**: Changes don't require recompiling a library
- **Template flexibility**: All solver code can be templated without export issues
- **Easy integration**: Users just include headers, no linking needed
- **Inline optimization**: Compiler can optimize across all boundaries

**Trade-offs:**
- **Compilation time**: Each translation unit compiles all included headers
- **Code bloat**: Potential for larger binaries

**Mitigation:**
- Single `main.cpp` keeps compilation fast for this application
- MFEM (the heavy dependency) is pre-compiled separately
- For larger projects, consider explicit template instantiation

**When to reconsider:**
- If compilation time becomes problematic (> 1 minute)
- If building a shared library for other projects
- Monitor: We track compile times in CI

Current compile time: ~30 seconds (acceptable for rapid development)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Copyright (c) 2026 T. C. Raymond

## Performance

For large-scale simulations:
- Use `Release` build type: `cmake .. -DCMAKE_BUILD_TYPE=Release`
- Enable OpenMP: `cmake .. -DUSE_OPENMP=ON` (default)
- Use configurable solver parameters in JSON:
  ```json
  {
    "simulation": {
      "solver_tolerance": 1e-10,
      "solver_max_iter": 2000,
      "solver_print_level": 1
    }
  }
  ```
- For > 100k DOFs, consider using HYPRE preconditioners (see `docs/solver_performance.md`)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

**Development workflow:**
1. Run tests: `ctest` in build directory
2. Check code compiles without warnings: `-Wall -Wextra`
3. Ensure examples still run correctly
4. Update documentation for API changes

## Acknowledgments

This project uses the [MFEM](https://mfem.org) library for finite element computations and [nlohmann/json](https://github.com/nlohmann/json) for JSON parsing.
