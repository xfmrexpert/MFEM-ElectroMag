// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

/**
 * @brief Physics formulation selected by "simulation.physics"
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
 * @brief Coordinate assumption / weak form ("simulation.geometry")
 */
enum class GeometryType {
    Axisymmetric,
    Planar
};

/// PhysicsType -> canonical JSON string (used by the factory and diagnostics).
inline const char* ToString(PhysicsType p) {
    switch (p) {
        case PhysicsType::Electrostatic:      return "electrostatics";
        case PhysicsType::Magnetostatic:      return "magnetostatics";
        case PhysicsType::Magnetoquasistatic: return "magnetoquasistatics";
    }
    return "unknown";
}
