// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include "mfem.hpp"
#include "json.hpp"
#include "problem_config.hpp"

using json = nlohmann::json;

/**
 * @brief Base class for physics solvers using MFEM
 *
 * @warning The mesh and config references must outlive this solver instance.
 *          Do not destroy the mesh or config objects before the solver is done.
 */
class PhysicsSolver {
protected:
    mfem::Mesh &mesh;
    const json &config_json;
    ProblemConfig config;

public:
    PhysicsSolver(mfem::Mesh &m, const json &c) : mesh(m), config_json(c) {}
    
    // Virtual destructor is essential for unique_ptr polymorphism
    virtual ~PhysicsSolver() = default;

    virtual void Setup() = 0;
    virtual void Run() = 0;
    virtual void SaveScenario(const std::string& scenario_name) = 0;
    virtual void SaveStudy() = 0;

    mfem::Array<int> MarkerFromAttrs(const std::vector<int>& attrs) const {
        mfem::Array<int> m(mesh.bdr_attributes.Max());
        m = 0;
        for (int a : attrs) {
            if (a > 0 && a <= m.Size()) m[a - 1] = 1;
        }
        return m;
    }

    // attr id (1-based) -> material index (0-based); -1 if no region claims it.
    std::vector<int> BuildAttrToMaterial() const {
        std::vector<int> m(mesh.attributes.Max() + 1, -1);
        for (const auto& region : config.Regions)
            for (int a : region.AttributeIds)
                if (a > 0 && a <= mesh.attributes.Max())
                    m[a] = region.Material;
        return m;
    }
};