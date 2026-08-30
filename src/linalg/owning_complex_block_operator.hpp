// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include "mfem.hpp"

// -----------------------------------------------------------------------------
// A ComplexOperator that OWNS its real and imaginary parts.
//
// mfem::ComplexOperator references (or, with ownReal/ownImag, deletes raw) its
// two real operators, leaving the caller to manage their lifetime and the
// strict "real/imag parts must outlive the ComplexOperator" ordering by hand.
// This wrapper takes unique_ptr ownership of both parts and the ComplexOperator
// built over them, so a single object governs the whole complex apparatus and
// tears it down in the correct order automatically.
//
// Either part may be null (an absent real or imaginary block), matching
// ComplexOperator's own contract.
// -----------------------------------------------------------------------------
class OwningComplexBlockOperator
{
public:
	OwningComplexBlockOperator(
		std::unique_ptr<mfem::Operator> real_op,
		std::unique_ptr<mfem::Operator> imag_op,
		mfem::ComplexOperator::Convention convention =
			mfem::ComplexOperator::HERMITIAN)
		: real_part(std::move(real_op)),
		  imag_part(std::move(imag_op))
	{
		MFEM_VERIFY(real_part || imag_part,
					"ComplexOperator needs at least a real or imaginary part.");
		system = std::make_unique<mfem::ComplexOperator>(
			real_part.get(), imag_part.get(),
			/*ownReal=*/false, /*ownImag=*/false, convention);
	}

	// The assembled complex operator; call FormLinearSystem/Mult on it.
	mfem::ComplexOperator& GetOperator()             { return *system; }
	const mfem::ComplexOperator& GetOperator() const { return *system; }

	mfem::Operator* RealPart() { return real_part.get(); }
	mfem::Operator* ImagPart() { return imag_part.get(); }

private:
	// 'system' references the two parts, so it is declared last and destroyed
	// first (members destruct in reverse declaration order).
	std::unique_ptr<mfem::Operator> real_part;
	std::unique_ptr<mfem::Operator> imag_part;
	std::unique_ptr<mfem::ComplexOperator> system;
};
