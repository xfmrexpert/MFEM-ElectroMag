// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "axis_geometry.hpp"

/**
 * @brief Axisymmetric magnetic field relations in the vector potential.
 *
 * Physics, not geometry: these turn A_phi and its derivatives into field
 * quantities. They live next to the geometry header because the axis limit
 * below must use the same scale-relative tolerance as the classification.
 */
namespace axisym {

/**
 * @brief Axial component of B = curl(A_phi e_phi), with the regular axis limit.
 *
 * Away from the axis this is the plain B_z = dA/dr + A_phi/r. Regularity forces
 * A_phi -> 0 linearly as r -> 0, so l'Hopital gives the finite limit
 * B_z -> 2 dA/dr on the axis instead of dividing by zero.
 *
 * The axis test uses the same scale-relative @p tolerance as the geometry
 * classification in axis_geometry.hpp, so a coordinate that the radial scan
 * calls "on the axis" also takes the limit here. An exact `r == 0.0` test would
 * leave a mesh whose axis nodes carry round-off evaluating `A_phi / r` at a
 * near-zero radius, which is unbounded rather than merely inaccurate.
 *
 * Valid only for the constrained SOLUTION, never for individual basis functions
 * before essential elimination: shape functions need not vanish on the axis, so
 * no per-basis limit exists.
 */
inline double AxialFluxDensity(double A_phi, double dA_dr, double r,
							   double tolerance)
{
   return IsOnAxisGeometry(r, tolerance) ? (2.0 * dA_dr)
										 : (dA_dr + A_phi / r);
}

} // namespace axisym
