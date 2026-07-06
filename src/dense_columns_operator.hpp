// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"

// -----------------------------------------------------------------------------
// A tall-and-thin (height x num_columns) operator stored as its columns:
//
//     C = [ c_0 | c_1 | ... | c_{n-1} ],   each c_i a Vector of length height.
//
// Mult maps a small coefficient vector x (size num_columns) to the weighted
// column sum  y = sum_i x_i c_i, and MultTranspose returns the inner products
// y_i = c_i . x. This is the generic rank-num_columns coupling used to border a
// large field operator with a few global unknowns (port voltages, circuit
// constraints, Lagrange multipliers, ...). The column vectors are referenced,
// not owned; the caller keeps them alive.
// -----------------------------------------------------------------------------
class DenseColumnsOperator : public mfem::Operator
{
public:
	// height: length of each column (the large/field dimension).
	// columns: one Vector pointer per column; the array is copied (the pointers
	// and the vectors they reference must outlive this operator).
	DenseColumnsOperator(int height, const mfem::Array<mfem::Vector*>& columns)
		: mfem::Operator(height, columns.Size()), columns(columns)
	{
		for (int i = 0; i < columns.Size(); ++i)
		{
			MFEM_ASSERT(columns[i] != nullptr, "Null column vector.");
			MFEM_ASSERT(columns[i]->Size() == height,
						"Column vector size must match operator height.");
		}
	}

	// y = C x = sum_i x_i * column_i   (size height)
	void Mult(const mfem::Vector& x, mfem::Vector& y) const override
	{
		MFEM_ASSERT(x.Size() == width,  "Input size must equal the column count.");
		MFEM_ASSERT(y.Size() == height, "Output size must equal the column length.");

		y = 0.0;
		for (int i = 0; i < width; ++i)
		{
			y.Add(x(i), *columns[i]);
		}
	}

	// y = C^T x,  y_i = column_i . x   (size num_columns)
	void MultTranspose(const mfem::Vector& x, mfem::Vector& y) const override
	{
		MFEM_ASSERT(x.Size() == height, "Input size must equal the column length.");
		MFEM_ASSERT(y.Size() == width,  "Output size must equal the column count.");

		for (int i = 0; i < width; ++i)
		{
			y(i) = (*columns[i]) * x;
		}
	}

private:
	mfem::Array<mfem::Vector*> columns;   // referenced, not owned
};
