// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <vector>
#include "mfem.hpp"
#include "bordered_block_operator.hpp"
#include "complex_block_layout.hpp"
#include "dense_columns_operator.hpp"
#include "owning_complex_block_operator.hpp"

// The assembled time-harmonic MQS operator with massive-conductor port coupling.
// M_omega_sigma already contains the omega*sigma scaling:
//
//   [ K + j M_omega_sigma    -C              ] [ A      ]   [ F_source    ]
//   [ -C^T                    -j/(omega G_dc) ] [ V_port ] = [ I/(j omega) ]
//
// The lower equation is the physical current balance
//
//     V_port/G_dc - j*omega*C^T*A = I,
//
// divided by j*omega to make the field/port coupling symmetric. Consequently
// the lower RHS is the scaled current phasor I/(j*omega), not the physical
// current phasor I itself. For a real in-phase current this value is -j*I/omega.
//
// K (curl-curl, 1/mu) and M_omega_sigma (omega*sigma mass) are the real/imag
// parts of the field SesquilinearForm; they are REFERENCED here, not owned (the
// caller keeps that form alive). The port load vectors (the columns of C) are
// owned here.
//
// This class owns only the MQS operator composition and its algebraic layout. It
// does not own scenario vectors, finite-element fields, or linear-solver policy.
// The wiring is delegated to reusable algebraic building blocks:
//   - DenseColumnsOperator       : C, the rank-N_ports border from port loads.
//   - BorderedBlockOperator      : the real [K, -C; -C^T, 0] and imaginary
//                                  [M_omega_sigma, 0; 0, -1/(omega G_dc)] blocks.
//   - OwningComplexBlockOperator : owns both blocks and the ComplexOperator.
// With zero massive ports the borders/corners are absent and the operator
// degenerates to K + j M_omega_sigma.
class MqsMassivePortOperator
{
public:
	MqsMassivePortOperator(
		int n_dofs,
		mfem::SparseMatrix& K,
		mfem::SparseMatrix& M_omega_sigma,
		std::vector<std::unique_ptr<mfem::Vector>> port_loads,
		const std::vector<mfem::real_t>& inverse_omega_conductance_diag)
		: layout(n_dofs, static_cast<int>(port_loads.size())),
		  port_loads(std::move(port_loads))
	{
		MFEM_VERIFY(
			static_cast<int>(inverse_omega_conductance_diag.size()) == layout.NPorts(),
			"Port self-admittance diagonal must have one entry per port.");

		std::unique_ptr<mfem::Operator> coupling;
		if (layout.NPorts() > 0)
		{
			mfem::Array<mfem::Vector*> column_ptrs(layout.NPorts());
			for (int p = 0; p < layout.NPorts(); ++p)
			{
				column_ptrs[p] = this->port_loads[p].get();
			}
			coupling = std::make_unique<DenseColumnsOperator>(layout.NDofs(), column_ptrs);
		}

		auto real_block = std::make_unique<BorderedBlockOperator>(
			layout.NDofs(), layout.NPorts(), K, std::move(coupling), -1.0,
			/*corner=*/nullptr);

		std::unique_ptr<mfem::Operator> port_corner;
		if (layout.NPorts() > 0)
		{
			auto inverse_omega_conductance =
				std::make_unique<mfem::DenseMatrix>(layout.NPorts(), layout.NPorts());
			*inverse_omega_conductance = 0.0;
			for (int p = 0; p < layout.NPorts(); ++p)
			{
				(*inverse_omega_conductance)(p, p) = inverse_omega_conductance_diag[p];
			}
			port_corner = std::move(inverse_omega_conductance);
		}

		auto imaginary_block = std::make_unique<BorderedBlockOperator>(
			layout.NDofs(), layout.NPorts(), M_omega_sigma,
			/*border=*/nullptr, 1.0, std::move(port_corner));

		complex_operator = std::make_unique<OwningComplexBlockOperator>(
			std::move(real_block), std::move(imaginary_block),
			mfem::ComplexOperator::HERMITIAN);
	}

	mfem::ComplexOperator& Operator() { return complex_operator->GetOperator(); }
	const mfem::ComplexOperator& Operator() const { return complex_operator->GetOperator(); }

	const ComplexPortLayout& Layout() const { return layout; }

	ComplexPortVectorView View(mfem::Vector& vector) const
	{
		return ComplexPortVectorView(vector, layout.NDofs(), layout.NPorts());
	}

	ConstComplexPortVectorView View(const mfem::Vector& vector) const
	{
		return ConstComplexPortVectorView(vector, layout.NDofs(), layout.NPorts());
	}

	mfem::Array<int> MakeEssentialTDofs(const mfem::Array<int>& field_tdofs) const
	{
		return ComplexEssentialTDofs(field_tdofs, layout.HalfSize());
	}

private:
	// complex_operator transitively references port_loads, so it is declared last
	// and destroyed first. layout is independent and remains the size authority.
	ComplexPortLayout layout;
	std::vector<std::unique_ptr<mfem::Vector>> port_loads;
	std::unique_ptr<OwningComplexBlockOperator> complex_operator;
};
