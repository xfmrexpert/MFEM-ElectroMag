// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <string>
#include "enums.hpp"

enum Quantity {Voltage, Current};

struct Region {
    std::vector<int> AttributeIds;
    int Material;
};

struct Port {
    int Region;
};

struct Material {
    double Conductivity;
    double RelPermittivity;
    double RelPermeability;
};

struct Source {
    double CurrentDensity;
    std::vector<int> Markers;
};

struct BoundaryCondition {
    std::string type;      // "Dirichlet", "Neumann", "Robin"
    std::vector<int> marker;
    double value;
    double robin_coeff;    // For Robin BCs: alpha * u + beta * du/dn = value

    BoundaryCondition(const std::string& t, const std::vector<int>& m, double v, double rc = 0.0)
        : type(t), marker(m), value(v), robin_coeff(rc) {}
};

struct ProblemConfig {
    int Order = 1;
    double Frequency = 60.0;
    ::ModelType ModelType = ::ModelType::Planar;
    double SolverTolerance;
    int SolverMaxIter;
    int SolverPrintLevel;
    std::string MeshPath;
    std::vector<Region> Regions;
    std::vector<Material> Materials;
    std::vector<Port> Ports;
    std::vector<BoundaryCondition> BoundaryConditions;
    std::vector<Source> Sources;
};