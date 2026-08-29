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
// M_sigma is assembled once and scaled by the active omega at operator-application
// time:
//
//   [ K + j M_omega_sigma    -C              ] [ A      ]   [ F_source    ]
//   [ -C^T                    -j G_dc/omega    ] [ V_port ] = [ I/(j omega) ]
//
// The lower equation is the physical current balance
//
//     G_dc*V_port - j*omega*C^T*A = I,
//
// divided by j*omega to make the field/port coupling symmetric. Consequently
// the lower RHS is the scaled current phasor I/(j*omega), not the physical
// current phasor I itself. For a real in-phase current this value is -j*I/omega.
//
// K (curl-curl, 1/mu) and M_sigma (sigma mass) are the real/imag
// parts of the field SesquilinearForm; they are REFERENCED here, not owned (the
// caller keeps that form alive). The port load vectors (the columns of C) are
// owned here.
//
// This class owns only the MQS operator composition and its algebraic layout. It
// does not own scenario vectors, finite-element fields, or linear-solver policy.
// The wiring is delegated to reusable algebraic building blocks:
//   - DenseColumnsOperator       : C, the rank-N_ports border from port loads.
//   - BorderedBlockOperator      : the real [K, -C; -C^T, 0] and imaginary
//                                  [omega M_sigma, 0; 0, -G_dc/omega] blocks.
//   - OwningComplexBlockOperator : owns both blocks and the ComplexOperator.
// With zero massive ports the borders/corners are absent and the operator
// degenerates to K + j M_omega_sigma.
class MqsMassivePortOperator
{
	class ScaledReferenceOperator : public mfem::Operator
	{
	public:
		ScaledReferenceOperator(mfem::Operator& referenced, mfem::real_t scale)
			: mfem::Operator(referenced.Height(), referenced.Width()),
			  referenced(referenced), scale(scale) {}

		void SetScale(mfem::real_t value) { scale = value; }

		void Mult(const mfem::Vector& x, mfem::Vector& y) const override
		{
			referenced.Mult(x, y);
			y *= scale;
		}

		void MultTranspose(const mfem::Vector& x, mfem::Vector& y) const override
		{
			referenced.MultTranspose(x, y);
			y *= scale;
		}

	private:
		mfem::Operator& referenced;
		mfem::real_t scale;
	};

	class ConductanceCornerOperator : public mfem::Operator
	{
	public:
		ConductanceCornerOperator(const std::vector<mfem::real_t>& conductances,
								  mfem::real_t omega)
			: mfem::Operator(static_cast<int>(conductances.size())),
			  conductances(conductances), omega(omega) {}

		void SetOmega(mfem::real_t value) { omega = value; }

		void Mult(const mfem::Vector& x, mfem::Vector& y) const override
		{
			for (int p = 0; p < Height(); ++p) {
				y(p) = -conductances[p] * x(p) / omega;
			}
		}

		void MultTranspose(const mfem::Vector& x, mfem::Vector& y) const override
		{
			Mult(x, y);
		}

	private:
		std::vector<mfem::real_t> conductances;
		mfem::real_t omega;
	};

public:
	MqsMassivePortOperator(
		int n_dofs,
		mfem::SparseMatrix& K,
		mfem::SparseMatrix& M_sigma,
		std::vector<std::unique_ptr<mfem::Vector>> port_loads,
		const std::vector<mfem::real_t>& conductances,
		mfem::real_t omega)
		: layout(n_dofs, static_cast<int>(port_loads.size())),
		  port_loads(std::move(port_loads)),
		  stiffness(K), sigma_mass(M_sigma), conductances(conductances),
		  active_omega(omega)
	{
		MFEM_VERIFY(omega > 0.0, "MQS angular frequency must be positive.");
		MFEM_VERIFY(
			static_cast<int>(conductances.size()) == layout.NPorts(),
			"Port conductance data must have one entry per port.");

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
			auto frequency_corner =
				std::make_unique<ConductanceCornerOperator>(conductances, omega);
			conductance_corner = frequency_corner.get();
			port_corner = std::move(frequency_corner);
		}

		auto frequency_mass = std::make_unique<ScaledReferenceOperator>(M_sigma, omega);
		auto imaginary_block = std::make_unique<BorderedBlockOperator>(
			layout.NDofs(), layout.NPorts(), *frequency_mass,
			/*border=*/nullptr, 1.0, std::move(port_corner));
		owned_scaled_mass = std::move(frequency_mass);

		complex_operator = std::make_unique<OwningComplexBlockOperator>(
			std::move(real_block), std::move(imaginary_block),
			mfem::ComplexOperator::HERMITIAN);
	}

	void SetOmega(mfem::real_t value)
	{
		MFEM_VERIFY(value > 0.0, "MQS angular frequency must be positive.");
		active_omega = value;
		owned_scaled_mass->SetScale(value);
		if (conductance_corner) conductance_corner->SetOmega(value);
	}

	// Assemble the packed real matrix that ComplexOperator applies matrix-free,
	// in the HERMITIAN layout its Mult() implements:
	//
	//     [ R  -I ] [ x_r ]   [ y_r ]
	//     [ I   R ] [ x_i ] = [ y_i ]
	//
	// where R = [K, -C; -C^T, 0] and I = [omega*M_sigma, 0; 0, -G_dc/omega] are
	// this operator's real and imaginary halves. Building it explicitly is the
	// price of admission for a direct solve: the block operators are apply-only,
	// so a factorization needs the coefficients gathered into one CSR matrix.
	//
	// The result is frequency-dependent through both omega-scaled blocks, so it
	// must be rebuilt (and refactored) whenever SetOmega() changes the frequency.
	std::unique_ptr<mfem::SparseMatrix> AssemblePackedMatrix() const
	{
		const int h = layout.HalfSize();
		const int n = layout.NDofs();
		auto packed = std::make_unique<mfem::SparseMatrix>(2 * h, 2 * h);

		// Real diagonal blocks: R at (0,0) and (h,h).
		AddSparseBlock(*packed, stiffness, 0, 0, 1.0);
		AddSparseBlock(*packed, stiffness, h, h, 1.0);

		// Imaginary blocks: -I at (0,h) and +I at (h,0). The sigma mass carries
		// the omega scaling that ScaledReferenceOperator applies at Mult time.
		AddSparseBlock(*packed, sigma_mass, 0, h, -active_omega);
		AddSparseBlock(*packed, sigma_mass, h, 0, active_omega);

		// Port border C occupies the field/port off-diagonals of R with sign -1,
		// mirrored into both real diagonal blocks.
		for (int p = 0; p < layout.NPorts(); ++p)
		{
			const mfem::Vector& column = *port_loads[p];
			for (int d = 0; d < n; ++d)
			{
				const mfem::real_t value = column(d);
				if (value == 0.0) { continue; }
				packed->Add(d, n + p, -value);
				packed->Add(n + p, d, -value);
				packed->Add(h + d, h + n + p, -value);
				packed->Add(h + n + p, h + d, -value);
			}

			// Conductance corner -G_dc/omega sits in the imaginary block only,
			// so it lands in the two off-diagonal halves with opposite signs.
			const mfem::real_t corner = -conductances[p] / active_omega;
			packed->Add(n + p, h + n + p, -corner);
			packed->Add(h + n + p, n + p, corner);
		}

		packed->Finalize();
		return packed;
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
	// Scatter scale*block into `target` at the given row/column origin. The
	// block is assumed finalized (CSR), which the SesquilinearForm guarantees.
	static void AddSparseBlock(mfem::SparseMatrix& target,
							   const mfem::SparseMatrix& block,
							   int row_offset, int col_offset,
							   mfem::real_t scale)
	{
		const int* I = block.GetI();
		const int* J = block.GetJ();
		const mfem::real_t* data = block.GetData();
		for (int r = 0; r < block.Height(); ++r)
		{
			for (int k = I[r]; k < I[r + 1]; ++k)
			{
				target.Add(row_offset + r, col_offset + J[k], scale * data[k]);
			}
		}
	}

	// complex_operator transitively references port_loads, so it is declared last
	// and destroyed first. layout is independent and remains the size authority.
	ComplexPortLayout layout;
	std::vector<std::unique_ptr<mfem::Vector>> port_loads;
	// Referenced (not owned) field matrices, retained so the packed matrix can be
	// reassembled at each new frequency.
	mfem::SparseMatrix& stiffness;
	mfem::SparseMatrix& sigma_mass;
	std::vector<mfem::real_t> conductances;
	mfem::real_t active_omega;
	ConductanceCornerOperator* conductance_corner = nullptr;
	std::unique_ptr<ScaledReferenceOperator> owned_scaled_mass;
	std::unique_ptr<OwningComplexBlockOperator> complex_operator;
};
