// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <vector>
#include "mfem.hpp"
#include "dense_columns_operator.hpp"
#include "bordered_block_operator.hpp"
#include "owning_complex_block_operator.hpp"

// The assembled time-harmonic MQS operator with massive-conductor port coupling:
//
//   [ K + jwM    -C   ] [ A      ]   [ F_source ]
//   [ -C^T       jG_s ] [ V_port ] = [ I_port   ]
//
// K (curl-curl, 1/mu) and M (omega*sigma mass) are the real/imag parts of the
// field SesquilinearForm; they are REFERENCED here, not owned (the caller keeps
// that form alive). The port load vectors (the columns of C) are owned here.
//
// The wiring is delegated to reusable building blocks:
//   - DenseColumnsOperator        : C, the rank-N_ports border from port loads.
//   - BorderedBlockOperator       : the real [K, -C; -C^T, 0] and imaginary
//                                   [M, 0; 0, G_s] saddle blocks; each owns its
//                                   border/corner and the inner BlockOperator.
//   - OwningComplexBlockOperator  : owns both saddle blocks and the
//                                   ComplexOperator built over them.
// With zero massive ports both borders/corner are absent and this degenerates
// to (K + jwM); the same plumbing serves 0..N ports with no separate path.
//
// The system is solved as one real vector of size FullSize() laid out
// [Re_Mesh, Re_Port, Im_Mesh, Im_Port]; HalfSize() exposes the Re/Im split so
// callers can pack the RHS and extract the solution.
class PortCoupledComplexSystem
{
public:
	// K, M: real/imag field matrices (referenced). port_loads: one column of C
	// per massive port (ownership taken). g_scaled_diag: matching diagonal port
	// self-admittance entries, sized to port_loads.size().
	PortCoupledComplexSystem(int n_dofs,
							 mfem::SparseMatrix& K,
							 mfem::SparseMatrix& M,
							 std::vector<std::unique_ptr<mfem::Vector>> port_loads,
							 const std::vector<mfem::real_t>& g_scaled_diag)
		: n_dofs(n_dofs),
		  n_ports(static_cast<int>(port_loads.size())),
		  port_loads(std::move(port_loads))
	{
		MFEM_VERIFY(static_cast<int>(g_scaled_diag.size()) == n_ports,
					"Port self-admittance diagonal must have one entry per port.");

		// Real saddle block: [ K  -C ; -C^T  0 ]. C is the dense-columns border
		// built from the owned port loads; the -1 sign is folded into the block
		// coefficient by BorderedBlockOperator. Corner is absent (zero block).
		std::unique_ptr<mfem::Operator> C_re;
		if (n_ports > 0)
		{
			mfem::Array<mfem::Vector*> form_ptrs(n_ports);
			for (int i = 0; i < n_ports; ++i) { form_ptrs[i] = this->port_loads[i].get(); }
			C_re = std::make_unique<DenseColumnsOperator>(n_dofs, form_ptrs);
		}
		auto real_block = std::make_unique<BorderedBlockOperator>(
			n_dofs, n_ports, K, std::move(C_re), -1.0, /*corner=*/nullptr);

		// Imag saddle block: [ M  0 ; 0  G_s ]. No border; the corner is the
		// diagonal port self-admittance matrix.
		std::unique_ptr<mfem::Operator> G_im;
		if (n_ports > 0)
		{
			auto G_scaled = std::make_unique<mfem::DenseMatrix>(n_ports, n_ports);
			*G_scaled = 0.0;
			for (int i = 0; i < n_ports; ++i) { (*G_scaled)(i, i) = g_scaled_diag[i]; }
			G_im = std::move(G_scaled);
		}
		auto imag_block = std::make_unique<BorderedBlockOperator>(
			n_dofs, n_ports, M, /*border=*/nullptr, 1.0, std::move(G_im));

		// Own the two saddle blocks and the ComplexOperator over them.
		system = std::make_unique<OwningComplexBlockOperator>(
			std::move(real_block), std::move(imag_block),
			mfem::ComplexOperator::HERMITIAN);
	}

	// The assembled complex block operator; call FormLinearSystem on it.
	mfem::ComplexOperator& GetOperator() { return system->GetOperator(); }

	int NDofs()    const { return n_dofs; }
	int NPorts()   const { return n_ports; }
	int HalfSize() const { return n_dofs + n_ports; }      // one (Re or Im) half
	int FullSize() const { return 2 * (n_dofs + n_ports); }

private:
	// 'system' (transitively, the real saddle block's DenseColumnsOperator)
	// references the owned port_loads, so port_loads is declared first and
	// outlives it (members destruct in reverse declaration order).
	int n_dofs;
	int n_ports;
	std::vector<std::unique_ptr<mfem::Vector>> port_loads;   // pointed-to by C
	std::unique_ptr<OwningComplexBlockOperator> system;
};
