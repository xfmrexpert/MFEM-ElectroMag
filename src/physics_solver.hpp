// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include "mfem.hpp"
#include "json.hpp"

using json = nlohmann::json;

class PhysicsSolver {
protected:
    mfem::Mesh &mesh;
    json &config;

    // Helper for generating consistent output names
    std::string GetOutputName(const std::string& base) {
        return base;
    }

public:
    PhysicsSolver(mfem::Mesh &m, json &c) : mesh(m), config(c) {}
    
    // Virtual destructor is essential for unique_ptr polymorphism
    virtual ~PhysicsSolver() = default;

    virtual void Setup() = 0;
    virtual void Run() = 0;
    virtual void Save() = 0;
};