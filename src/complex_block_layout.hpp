// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"

// -----------------------------------------------------------------------------
// Helpers for the packed real vector used to solve a complex, port-coupled
// block system. ComplexOperator stores the real part first, then the imaginary
// part, so a system with N_DOFs field unknowns and N_Ports port unknowns is
// solved as ONE real vector laid out:
//
//     [ Re_Mesh (N_DOFs) | Re_Port (N_Ports) | Im_Mesh (N_DOFs) | Im_Port (N_Ports) ]
//       \------------- half size = N_DOFs + N_Ports -------------/
//
// These classes name the four slots so callers never hand-roll the "d + h" /
// "h + N_DOFs + p" arithmetic that is so easy to get wrong.
// -----------------------------------------------------------------------------

// Index arithmetic shared by the mutable and const views. Carries no data
// beyond the two sizes, so it is cheap to copy and safe to embed.
class ComplexPortLayout
{
public:
	ComplexPortLayout(int n_dofs, int n_ports)
		: n_dofs(n_dofs), n_ports(n_ports) {}

	int NDofs()    const { return n_dofs; }
	int NPorts()   const { return n_ports; }
	int HalfSize() const { return n_dofs + n_ports; }      // one (Re or Im) half
	int FullSize() const { return 2 * (n_dofs + n_ports); }

	// Packed indices for each named slot.
	int ReMeshIndex(int d) const { return d; }
	int RePortIndex(int p) const { return n_dofs + p; }
	int ImMeshIndex(int d) const { return HalfSize() + d; }
	int ImPortIndex(int p) const { return HalfSize() + n_dofs + p; }

protected:
	int n_dofs;
	int n_ports;
};

// Mutable, non-owning view over a packed solver vector. Use it to assemble the
// RHS and to lift essential values into the initial guess.
class ComplexPortVectorView : public ComplexPortLayout
{
public:
	ComplexPortVectorView(mfem::Vector& v, int n_dofs, int n_ports)
		: ComplexPortLayout(n_dofs, n_ports), v(v)
	{
		MFEM_ASSERT(v.Size() == FullSize(),
					"Vector size does not match the packed complex+port layout.");
	}

	mfem::real_t& ReMesh(int d) { return v(ReMeshIndex(d)); }
	mfem::real_t& RePort(int p) { return v(RePortIndex(p)); }
	mfem::real_t& ImMesh(int d) { return v(ImMeshIndex(d)); }
	mfem::real_t& ImPort(int p) { return v(ImPortIndex(p)); }

	mfem::real_t ReMesh(int d) const { return v(ReMeshIndex(d)); }
	mfem::real_t RePort(int p) const { return v(RePortIndex(p)); }
	mfem::real_t ImMesh(int d) const { return v(ImMeshIndex(d)); }
	mfem::real_t ImPort(int p) const { return v(ImPortIndex(p)); }

private:
	mfem::Vector& v;
};

// Read-only, non-owning view over a packed solver vector. Use it to extract the
// solved field DOFs back out of the monolithic solution.
class ConstComplexPortVectorView : public ComplexPortLayout
{
public:
	ConstComplexPortVectorView(const mfem::Vector& v, int n_dofs, int n_ports)
		: ComplexPortLayout(n_dofs, n_ports), v(v)
	{
		MFEM_ASSERT(v.Size() == FullSize(),
					"Vector size does not match the packed complex+port layout.");
	}

	mfem::real_t ReMesh(int d) const { return v(ReMeshIndex(d)); }
	mfem::real_t RePort(int p) const { return v(RePortIndex(p)); }
	mfem::real_t ImMesh(int d) const { return v(ImMeshIndex(d)); }
	mfem::real_t ImPort(int p) const { return v(ImPortIndex(p)); }

private:
	const mfem::Vector& v;
};

// Map scalar-space essential true DOFs onto the packed complex layout's
// essential-DOF list. Each real-space essential DOF d constrains BOTH its real
// copy (index d) and its imaginary copy (index d + half_size). Field DOFs are
// the only essential ones (ports are never constrained), so only the mesh block
// is duplicated. half_size is N_DOFs + N_Ports, i.e. ComplexPortLayout::HalfSize().
inline mfem::Array<int> ComplexEssentialTDofs(const mfem::Array<int>& ess_real,
											  int half_size)
{
	mfem::Array<int> ess_complex;
	ess_complex.Reserve(2 * ess_real.Size());
	for (int k = 0; k < ess_real.Size(); ++k)
	{
		const int d = ess_real[k];
		ess_complex.Append(d);              // Re copy
		ess_complex.Append(d + half_size);  // Im copy
	}
	return ess_complex;
}
