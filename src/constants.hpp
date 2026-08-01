// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

namespace Constants {
    // Physical constants
    constexpr double MU_0 = 4.0 * 3.141592653589793 * 1e-7;  // Permeability of free space [H/m]
    constexpr double EPSILON_0 = 8.854187817e-12;             // Permittivity of free space [F/m]
    constexpr double TWO_PI = 2.0 * 3.141592653589793;        // 2*pi

    // Numerical tolerances
    // NOTE: axis proximity is no longer governed by a fixed absolute tolerance.
    // See axisym::kRelativeGeometryTolerance in axisymmetric_mesh.hpp, which
    // scales with the mesh bounding box so the same model works in any units.

    // Default solver parameters
    constexpr double DEFAULT_SOLVER_TOLERANCE = 1e-12;
    constexpr int DEFAULT_SOLVER_MAX_ITER = 1000;
    constexpr int DEFAULT_SOLVER_PRINT_LEVEL = 1;
}
