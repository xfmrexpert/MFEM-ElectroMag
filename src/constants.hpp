// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

namespace Constants {
    // Physical constants
    constexpr double MU_0 = 4.0 * 3.141592653589793 * 1e-7;  // Permeability of free space [H/m]
    constexpr double EPSILON_0 = 8.854187817e-12;             // Permittivity of free space [F/m]
    constexpr double TWO_PI = 2.0 * 3.141592653589793;        // 2*pi

    // Numerical tolerances
    constexpr double AXIS_TOLERANCE = 0.0;                    // Disabled (using mesh shifting instead)
                                                              // Set to 0 to disable weighted residual near axis

    // Default solver parameters
    constexpr double DEFAULT_SOLVER_TOLERANCE = 1e-12;
    constexpr int DEFAULT_SOLVER_MAX_ITER = 1000;
    constexpr int DEFAULT_SOLVER_PRINT_LEVEL = 1;
}
