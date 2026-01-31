// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

/**
 * @brief Physics simulation types
 */
enum class PhysicsType {
    Electrostatic,
    Magnetostatic,
    Magnetoquasistatic
};

/**
 * @brief Boundary condition types
 */
enum class BoundaryType {
    Dirichlet,
    Neumann,
    Robin
};

/**
 * @brief Solver model types (geometry)
 */
enum class ModelType {
    Axisymmetric,
    Planar
};
