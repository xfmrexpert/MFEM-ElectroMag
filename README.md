# MFEM-ElectroMag

A finite element solver for electromagnetic problems using MFEM (Modular Finite Element Methods). This solver supports electrostatic, magnetostatic, and magnetoquasistatic problems in both axisymmetric and planar geometries.

## Features

- **Electrostatics**: Solves for electric potential and field distributions
- **Magnetostatics**: Solves for magnetic vector potential and field distributions  
- **Magnetoquasistatics**: Time-harmonic eddy current problems
- **Axisymmetric and Planar**: Supports both 2D coordinate systems
- **JSON Configuration**: Easy problem setup via JSON files
- **ParaView Output**: Direct visualization of results

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

After building, the executable `solve_axi` will be in the `build` directory:

```bash
# Run with a configuration file
./solve_axi path/to/config.json

# Example: Run test cases
./solve_axi ../test/electrostatic_test.json
./solve_axi ../test/magnetostatic_test.json
./solve_axi ../test/mqs_test.json
```

## Configuration

Problems are configured using JSON files. See the `test/` directory for examples of:
- `electrostatic_test.json`: Electrostatic problem setup
- `magnetostatic_test.json`: Magnetostatic problem setup
- `mqs_test.json`: Magnetoquasistatic problem setup

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
│   └── axisymmetric_*.hpp         # Custom axisymmetric integrators
├── test/                          # Test configurations
├── CMakeLists.txt                 # CMake build configuration
├── LICENSE                        # MIT License
└── README.md                      # This file
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Copyright (c) 2026 T. C. Raymond

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Acknowledgments

This project uses the [MFEM](https://mfem.org) library for finite element computations and [nlohmann/json](https://github.com/nlohmann/json) for JSON parsing.
