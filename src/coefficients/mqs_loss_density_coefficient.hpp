// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT
//
// Time-averaged Joule (eddy-current) loss density for the magnetoquasistatic
// formulation, using the peak (amplitude) phasor convention:
//
//     P = 0.5 * sigma * |E|^2,   E = E_drive - j*omega*A
//
// The 1/2 is the time average of a peak-amplitude phasor product, matching the
// convention used by the coupling-matrix extraction (R = Re(V/I)) and by the
// common commercial AC/DC and eddy-current tools.
//
// Splitting E into real and imaginary parts with A = A_re + j*A_im:
//
//     -j*omega*A = omega*A_im - j*omega*A_re
//
// so
//
//     Re E = E_drive_re + omega*A_im
//     Im E = E_drive_im - omega*A_re
//
// IMPORTANT (why E_drive cannot be dropped):
//   The simplification P = 0.5*sigma*omega^2*|A|^2 is valid ONLY where
//   E_drive = 0. That holds for conductive regions with no port unknown (a flux
//   shield or a steel brace), but NOT for massive terminals or open-current
//   regions, whose solved port voltage is exactly what enforces their current
//   constraint. Using the simplified form there silently under- or over-reports
//   the loss, so the drive term is carried explicitly.
//
// The per-attribute drive table is zero-initialised and populated only where a
// port unknown exists. Unported conductors therefore fall out of the general
// expression with no special case, which is why this one coefficient serves
// driven terminals, open-current regions, and passive conductors alike.
//
// Geometry enters in one place only. For an azimuthal massive conductor the
// drive field produced by a port voltage V is
//
//     E_phi = V / (2*pi*r)
//
// the same physical 1/r that appears in AxisymmetricConductanceCoeff. That
// factor is a property of the field, entirely distinct from the 2*pi*r volume
// measure supplied by the integrator when this density is integrated.

#pragma once

#include <vector>

#include "mfem.hpp"
#include "../core/constants.hpp"

class MqsLossDensityCoefficient : public mfem::Coefficient {
public:
	// @param sigma       Piecewise conductivity; also selects which regions
	//                    dissipate at all (sigma = 0 yields exactly zero).
	// @param a_re,a_im   Real and imaginary parts of the vector potential.
	// @param omega       Angular frequency [rad/s].
	// @param drive_re    Per-attribute real part of the drive amplitude, indexed
	// @param drive_im    by (attribute - 1). Both are zero for unported regions.
	// @param axisymmetric  When true the drive amplitude carries the physical
	//                      1/(2*pi*r) factor.
	MqsLossDensityCoefficient(mfem::Coefficient& sigma,
							  const mfem::GridFunction& a_re,
							  const mfem::GridFunction& a_im,
							  double omega,
							  std::vector<double> drive_re,
							  std::vector<double> drive_im,
							  bool axisymmetric)
		: sigma_(sigma), a_re_(a_re), a_im_(a_im), omega_(omega),
		  drive_re_(std::move(drive_re)), drive_im_(std::move(drive_im)),
		  axisymmetric_(axisymmetric) {
		MFEM_ASSERT(drive_re_.size() == drive_im_.size(),
					"Drive amplitude tables must have matching sizes");
	}

	double Eval(mfem::ElementTransformation& T,
				const mfem::IntegrationPoint& ip) override {
		T.SetIntPoint(&ip);

		const double sigma = sigma_.Eval(T, ip);
		// Non-conducting regions dissipate nothing. Returning early also keeps
		// the axis guard below from firing on regions that cannot dissipate.
		if (sigma <= 0.0) { return 0.0; }

		double e_re = 0.0;
		double e_im = 0.0;
		DriveField(T, ip, e_re, e_im);

		e_re += omega_ * a_im_.GetValue(T, ip);
		e_im -= omega_ * a_re_.GetValue(T, ip);

		return 0.5 * sigma * (e_re * e_re + e_im * e_im);
	}

private:
	// Drive field for the element's attribute, zero where no port unknown
	// exists. Attributes are 1-based in MFEM; an attribute beyond the table is
	// treated as unported rather than as an error, since the table is sized to
	// the mesh and a larger attribute cannot carry a port.
	void DriveField(mfem::ElementTransformation& T,
					const mfem::IntegrationPoint& ip,
					double& e_re, double& e_im) const {
		const int index = T.Attribute - 1;
		if (index < 0 || index >= static_cast<int>(drive_re_.size())) { return; }

		e_re = drive_re_[index];
		e_im = drive_im_[index];
		if (e_re == 0.0 && e_im == 0.0) { return; }

		if (!axisymmetric_) { return; }

		mfem::Vector pos;
		T.Transform(ip, pos);
		const double r = pos(0);
		// Only reachable for a ported region, and a massive port that reaches
		// the axis is rejected at assembly time because its DC conductance
		// integral diverges. Hitting this means that guard was bypassed.
		MFEM_VERIFY(r > 0.0,
			"Axisymmetric drive field V/(2*pi*r) is singular at r = "
			<< r << ". A driven azimuthal conductor cannot touch the "
			"symmetry axis.");

		const double scale = 1.0 / (Constants::TWO_PI * r);
		e_re *= scale;
		e_im *= scale;
	}

	mfem::Coefficient& sigma_;
	const mfem::GridFunction& a_re_;
	const mfem::GridFunction& a_im_;
	double omega_;
	std::vector<double> drive_re_;
	std::vector<double> drive_im_;
	bool axisymmetric_;
};
