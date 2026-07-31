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
cmake .. -DMFEM_USE_METIS=ON
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

### MQS Frequency Scenarios

Magnetoquasistatic frequency is defined on every scenario, not in the
`simulation` block. A scalar performs one solve:

```json
{"name": "60Hz", "frequency": 60.0, "excitations": []}
```

A sweep uses an inclusive linear or logarithmic range:

```json
{
  "name": "Sweep",
  "frequency": {"scale": "log", "start": 10.0, "stop": 1000.0, "points": 5},
  "excitations": []
}
```

`points` includes both endpoints; with `points: 1`, only `start` is solved.
Field analyses write one result per expanded point. For MQS
`coupling_matrix` analyses, each scenario defines a frequency point and the
solver writes a separate frequency-labeled resistance/inductance CSV pair.

For detailed mathematical formulation, see `docs/math_formulation.md`.

## Output

Results are saved in ParaView format (`.pvd`, `.vtu` files) and can be visualized using:
- [ParaView](https://www.paraview.org/): Open-source visualization application
- [VisIt](https://visit-dav.github.io/visit-website/): Scientific data visualization tool

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
