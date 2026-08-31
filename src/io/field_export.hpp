// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT
//
// FieldExportSet: the single, format-agnostic description of WHAT a solver
// wants serialized for a scenario. It is built once in a post-solve step
// (PhysicsSolver::CollectExportFields) and handed to format-specific writers
// (field_writers.hpp) that decide HOW to sample/serialize it.
//
// A field is one of:
//   - Primary: a solution GridFunction, borrowed. Writers that can render it at
//     native (high) order do so (ParaView); others resample. The descriptor
//     names the field's ROLE, not its representation: a GridFunction may be
//     scalar or vector-valued (VDim > 1), and writers must branch on the
//     GridFunction rather than on the Kind to tell them apart.
//   - DerivedScalar:  a scalar Coefficient (material property, |B|, ...).
//   - DerivedVector:  a vector Coefficient (E = -grad V, B = curl A, ...).
//
// Derived coefficients are LAZY: nothing is projected here. Each writer projects
// them onto its own target mesh, because the formats deliberately sample on
// different meshes (ParaView: native; Gmsh: refined-to-linear). "Compute once"
// therefore means one field DEFINITION, not one numeric projection.
//
// The set OWNS the derived coefficients (and any intermediate coefficients they
// reference, parked via Own()) so they outlive the per-format projection. Raw
// pointers in the descriptors target heap-allocated coefficients, so the set is
// safely movable (returned by value from CollectExportFields).

#pragma once

#include <memory>
#include <string>
#include <vector>
#include "mfem.hpp"

struct FieldExport {
	enum class Kind { Primary, DerivedScalar, DerivedVector };

	Kind        kind;
	std::string name;
	mfem::GridFunction*      primary = nullptr; // Primary: borrowed solution GF
	mfem::Coefficient*       scalar  = nullptr; // DerivedScalar
	mfem::VectorCoefficient* vector  = nullptr; // DerivedVector
};

class FieldExportSet {
public:
	// A borrowed solution GridFunction. Caller guarantees it outlives the write.
	void AddPrimary(std::string name, mfem::GridFunction& gf) {
		fields_.push_back(
			{ FieldExport::Kind::Primary, std::move(name), &gf, nullptr, nullptr });
	}

	// A derived scalar field; the set takes ownership of the coefficient and
	// returns a reference so it can be wired into later fields.
	mfem::Coefficient& AddScalar(std::string name,
								 std::unique_ptr<mfem::Coefficient> coeff) {
		mfem::Coefficient& ref = *coeff;
		owned_scalars_.push_back(std::move(coeff));
		fields_.push_back(
			{ FieldExport::Kind::DerivedScalar, std::move(name), nullptr, &ref, nullptr });
		return ref;
	}

	// A derived scalar field backed by a coefficient the caller keeps owning
	// (e.g. a long-lived material-property member).
	void AddScalar(std::string name, mfem::Coefficient& borrowed) {
		fields_.push_back(
			{ FieldExport::Kind::DerivedScalar, std::move(name), nullptr, &borrowed, nullptr });
	}

	// A derived vector field; the set takes ownership and returns a reference so
	// it can feed a later field (e.g. a magnitude built from two vectors).
	mfem::VectorCoefficient& AddVector(std::string name,
									   std::unique_ptr<mfem::VectorCoefficient> coeff) {
		mfem::VectorCoefficient& ref = *coeff;
		owned_vectors_.push_back(std::move(coeff));
		fields_.push_back(
			{ FieldExport::Kind::DerivedVector, std::move(name), nullptr, nullptr, &ref });
		return ref;
	}

	// Park an intermediate coefficient that a later field references but that is
	// not itself exported (e.g. the gradient feeding a scaled-product vector).
	mfem::VectorCoefficient& Own(std::unique_ptr<mfem::VectorCoefficient> coeff) {
		mfem::VectorCoefficient& ref = *coeff;
		owned_vectors_.push_back(std::move(coeff));
		return ref;
	}
	mfem::Coefficient& Own(std::unique_ptr<mfem::Coefficient> coeff) {
		mfem::Coefficient& ref = *coeff;
		owned_scalars_.push_back(std::move(coeff));
		return ref;
	}

	const std::vector<FieldExport>& Fields() const { return fields_; }
	bool Empty() const { return fields_.empty(); }

private:
	std::vector<FieldExport> fields_;
	std::vector<std::unique_ptr<mfem::Coefficient>>       owned_scalars_;
	std::vector<std::unique_ptr<mfem::VectorCoefficient>> owned_vectors_;
};
