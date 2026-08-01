// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "mfem.hpp"
#include "problem_config.hpp"

struct MarkedBoundaryCondition {
	mfem::Array<int> Marker;
	BoundaryCondition Condition;
};
