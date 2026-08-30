// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include "mfem.hpp"

// -----------------------------------------------------------------------------
// A 2x2 "bordered" (saddle-point) operator over a large primary block A bordered
// by a few global unknowns:
//
//     [ A        s*C ]
//     [ s*C^T    D   ]
//
// A is the primary/field operator (REFERENCED, kept alive by the caller). C is
// the rectangular border (n_primary x n_border) and D is the small corner block
// (n_border x n_border); BOTH are OWNED here. The scalar s (border_sign) is
// folded into the off-diagonal block coefficients so the border itself need not
// be pre-scaled (the real saddle block uses s = -1).
//
// Either the border or the corner may be null:
//   - null border  -> the (0,1) and (1,0) blocks are zero,
//   - null corner  -> the (1,1) block is zero,
//   - n_border == 0 -> the operator degenerates to just A, with no extra code path.
//
// This bundles all the auxiliary pieces (border, its transpose, the corner) with
// the inner BlockOperator under a SINGLE owner, removing the "declaration order
// is load-bearing" reasoning that hand-built block systems otherwise require.
// It derives from mfem::Operator and delegates the action to the inner block, so
// it can be used anywhere an Operator is expected (e.g. as the real or imaginary
// part of a ComplexOperator).
// -----------------------------------------------------------------------------
class BorderedBlockOperator : public mfem::Operator
{
public:
	BorderedBlockOperator(int n_primary, int n_border,
						  mfem::Operator& field,
						  std::unique_ptr<mfem::Operator> border,
						  mfem::real_t border_sign,
						  std::unique_ptr<mfem::Operator> corner)
		: mfem::Operator(n_primary + n_border),
		  border(std::move(border)),
		  corner(std::move(corner))
	{
		MFEM_VERIFY(field.Height() == n_primary && field.Width() == n_primary,
					"Field block must be n_primary x n_primary.");
		if (this->border)
		{
			MFEM_VERIFY(this->border->Height() == n_primary &&
						this->border->Width()  == n_border,
						"Border block must be n_primary x n_border.");
			border_transpose =
				std::make_unique<mfem::TransposeOperator>(this->border.get());
		}
		if (this->corner)
		{
			MFEM_VERIFY(this->corner->Height() == n_border &&
						this->corner->Width()  == n_border,
						"Corner block must be n_border x n_border.");
		}

		block_offsets.SetSize(3);
		block_offsets[0] = 0;
		block_offsets[1] = n_primary;
		block_offsets[2] = n_primary + n_border;

		block = std::make_unique<mfem::BlockOperator>(block_offsets);
		block->SetBlock(0, 0, &field);
		if (this->border)
		{
			block->SetBlock(0, 1, this->border.get(),       border_sign);
			block->SetBlock(1, 0, border_transpose.get(),   border_sign);
		}
		if (this->corner)
		{
			block->SetBlock(1, 1, this->corner.get());
		}
	}

	void Mult(const mfem::Vector& x, mfem::Vector& y) const override
	{
		block->Mult(x, y);
	}

	void MultTranspose(const mfem::Vector& x, mfem::Vector& y) const override
	{
		block->MultTranspose(x, y);
	}

	mfem::BlockOperator& Block() { return *block; }

private:
	// Declaration order is reverse of destruction order. 'block' references the
	// external field plus the members below it, so it is declared last and torn
	// down first; 'border_transpose' references 'border', so it follows it.
	mfem::Array<int> block_offsets;
	std::unique_ptr<mfem::Operator> border;             // owned (may be null)
	std::unique_ptr<mfem::TransposeOperator> border_transpose; // refs border
	std::unique_ptr<mfem::Operator> corner;             // owned (may be null)
	std::unique_ptr<mfem::BlockOperator> block;         // refs field + above
};
