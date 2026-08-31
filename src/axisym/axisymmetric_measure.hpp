// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"
#include "../core/constants.hpp"

// -----------------------------------------------------------------------------
// The axisymmetric integration measure - the single definition in the project.
// -----------------------------------------------------------------------------
// Revolving a meridional (r, z) domain through the full azimuthal angle gives
// the volume element
//
//     dV = 2*pi*r dr dz,
//
// and the corresponding meridional boundary element 2*pi*r ds. EVERY
// axisymmetric integrator in this project applies the full measure via
// Axisymmetric::Measure(r) below. Nothing is factored out and restored later.
//
// The consequence, and the reason this convention was chosen: every assembled
// matrix, every right-hand side, and every quantity derived from them is in SI
// units and is directly comparable to a hand calculation. A stiffness entry is
// a real capacitance-like quantity, G_dc is a real conductance in siemens, and
// a flux linkage is real webers - with no per-call-site scale factor to
// remember. Adding a new derived quantity requires no knowledge of any
// normalization convention at all.
//
// NOT part of this measure, and deliberately left alone:
//   - AxisymmetricConductanceCoeff's 1/(2*pi*r), which is physical: it comes
//     from E_phi = V/(2*pi*r) for an azimuthal conductor, not from the volume
//     element.
namespace Axisymmetric
{
   // Geometric measure at radius r, excluding the quadrature weight and the
   // Jacobian determinant.
   //
   // Typed in mfem::real_t rather than double because the result is multiplied
   // straight into MFEM quadrature weights. MFEM is precision-configurable
   // (see config.hpp: real_t is float under MFEM_USE_SINGLE), so hardcoding
   // double here would silently widen every weight in a single-precision build.
   inline mfem::real_t Measure(mfem::real_t r)
   {
	  return static_cast<mfem::real_t>(Constants::TWO_PI) * r;
   }
}
