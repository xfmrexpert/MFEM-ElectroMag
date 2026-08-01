// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"

/**
 * @brief Integrates an axisymmetric natural boundary load.
 *
 * Assembles integral(q * v * r ds) on the meridional boundary. The global
 * 2*pi factor is omitted consistently with the other axisymmetric integrators.
 */
class AxisymmetricBoundaryLFIntegrator : public mfem::LinearFormIntegrator
{
private:
	mfem::Coefficient* load;

public:
	explicit AxisymmetricBoundaryLFIntegrator(
		mfem::Coefficient& q, const mfem::IntegrationRule* ir = nullptr)
		: mfem::LinearFormIntegrator(ir), load(&q)
	{
		MFEM_ASSERT(load != nullptr, "Coefficient cannot be null");
	}

	void AssembleRHSElementVect(const mfem::FiniteElement& el,
								mfem::ElementTransformation& trans,
								mfem::Vector& elvect) override
	{
		const int nd = el.GetDof();
		elvect.SetSize(nd);
		elvect = 0.0;

		mfem::Vector shape(nd);
		mfem::Vector pos(trans.GetSpaceDim());
		const mfem::IntegrationRule* ir = GetIntegrationRule(el, trans);

		for (int i = 0; i < ir->GetNPoints(); ++i)
		{
			const mfem::IntegrationPoint& ip = ir->IntPoint(i);
			trans.SetIntPoint(&ip);
			trans.Transform(ip, pos);

			const double r = pos(0);
			const double value = load->Eval(trans, ip);
			const double weight = ip.weight * trans.Weight() * r * value;

			el.CalcShape(ip, shape);
			elvect.Add(weight, shape);
		}
	}

protected:
	const mfem::IntegrationRule* GetDefaultIntegrationRule(
		const mfem::FiniteElement& trial_fe,
		const mfem::FiniteElement&,
		const mfem::ElementTransformation& trans) const override
	{
		const int order = trial_fe.GetOrder() + trans.Order() + trans.OrderW();
		if (trial_fe.Space() == mfem::FunctionSpace::rQk)
		{
			return &mfem::RefinedIntRules.Get(trial_fe.GetGeomType(), order);
		}
		return &mfem::IntRules.Get(trial_fe.GetGeomType(), order);
	}
};
