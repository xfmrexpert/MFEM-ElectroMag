// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

/**
 * @brief Physics formulation selected by "simulation.physics"
 */
enum class PhysicsType {
    Electrostatics,
    Magnetostatics,
    Magnetoquasistatics
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
        case PhysicsType::Electrostatics:      return "electrostatics";
        case PhysicsType::Magnetostatics:      return "magnetostatics";
        case PhysicsType::Magnetoquasistatics: return "magnetoquasistatics";
    }
    return "unknown";
}
