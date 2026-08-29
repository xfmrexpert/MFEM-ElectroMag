// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include "../src/config_validator.hpp"
#include "../src/boundary_validation.hpp"
#include <limits>

namespace {

bool HasError(const ConfigValidator& validator, const std::string& field);

json ValidConfig() {
	return json{
		{"simulation", {
			{"physics_type", "electrostatics"},
			{"geometry_type", "planar"},
			{"analysis_type", "field"},
			{"mesh", "test.msh"},
			{"order", 1}
		}},
		{"entity_groups", json::array({
			{{"name", "Domain"}, {"dim", 2}, {"attribute_ids", {1}}},
			{{"name", "Boundary"}, {"dim", 1}, {"attribute_ids", {1}}}
		})},
		{"regions", json::array({
			{{"entity_group", "Domain"}, {"material", "Dielectric"}}
		})},
		{"materials", json::array({
			{{"name", "Dielectric"}, {"properties", {{"epsilon_r", 1.0}}}}
		})},
		{"terminals", json::array({
			{{"name", "Drive"}, {"excitation_type", "voltage"}, {"entity_group", "Boundary"}}
		})},
		{"boundaries", json::array({
			{{"type", "Dirichlet"}, {"entity_group", "Boundary"}, {"value", 0.0}}
		})},
		{"scenarios", json::array({
			{{"name", "Driven"}, {"excitations", json::array({
				{{"terminal", "Drive"}, {"value", 1.0}}
			})}}
		})}
	};
}

TEST_CASE("BoundaryConditionValidator distinguishes essential and weak closures",
		  "[boundary_validation]") {
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		1, 1, mfem::Element::QUADRILATERAL, true, 1.0, 1.0);
	for (int i = 0; i < mesh.GetNBE(); ++i) {
		mesh.GetBdrElement(i)->SetAttribute(i + 1);
	}

	mfem::H1_FECollection fec(1, mesh.Dimension());
	mfem::FiniteElementSpace fespace(&mesh, &fec);
	BoundaryConditionValidator validator(mesh, fespace);
	const std::unordered_map<std::string, mfem::Array<int>> no_terminals;

	auto marker = [&](int attribute) {
		mfem::Array<int> result(mesh.GetNBE());
		result = 0;
		result[attribute - 1] = 1;
		return result;
	};
	int adjacent_attribute = -1;
	mfem::Array<int> first_vertices;
	mesh.GetBdrElementVertices(0, first_vertices);
	for (int boundary = 1; boundary < mesh.GetNBE() && adjacent_attribute < 0;
		 ++boundary) {
		mfem::Array<int> candidate_vertices;
		mesh.GetBdrElementVertices(boundary, candidate_vertices);
		for (int i = 0; i < first_vertices.Size(); ++i) {
			for (int j = 0; j < candidate_vertices.Size(); ++j) {
				if (first_vertices[i] == candidate_vertices[j]) {
					adjacent_attribute = boundary + 1;
				}
			}
		}
	}
	REQUIRE(adjacent_attribute > 0);

	const BoundaryCondition dirichlet_zero(
		BoundaryConditionType::Dirichlet, "First", 0.0);
	const BoundaryCondition dirichlet_one(
		BoundaryConditionType::Dirichlet, "Second", 1.0);
	const BoundaryCondition neumann(
		BoundaryConditionType::Neumann, "Second", 2.0);

	SECTION("allows a Dirichlet and Neumann boundary junction") {
		const std::vector<MarkedBoundaryCondition> closures = {
			{marker(1), dirichlet_zero}, {marker(adjacent_attribute), neumann}
		};
		REQUIRE_NOTHROW(validator.ValidateBoundaryConditions(
			closures, no_terminals));
	}

	SECTION("rejects different Dirichlet values at a shared corner") {
		const std::vector<MarkedBoundaryCondition> closures = {
			{marker(1), dirichlet_zero},
			{marker(adjacent_attribute), dirichlet_one}
		};
		REQUIRE_THROWS(validator.ValidateBoundaryConditions(
			closures, no_terminals));
	}

	SECTION("rejects duplicate assignments on one boundary attribute") {
		const std::vector<MarkedBoundaryCondition> closures = {
			{marker(1), dirichlet_zero}, {marker(1), neumann}
		};
		REQUIRE_THROWS(validator.ValidateBoundaryConditions(
			closures, no_terminals));
	}
}

TEST_CASE("ConfigValidator validates boundary values and Robin metadata",
		  "[config_validator][boundaries]") {
	SECTION("accepts Neumann without a Robin coefficient") {
		json config = ValidConfig();
		config["boundaries"][0]["type"] = "Neumann";
		config["boundaries"][0]["value"] = -2.5;

		ConfigValidator validator;
		REQUIRE(validator.Validate(config));
	}

	SECTION("rejects non-finite boundary values") {
		json config = ValidConfig();
		config["boundaries"][0]["value"] =
			std::numeric_limits<double>::infinity();

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "boundaries[0].value"));
	}

	SECTION("requires a finite coefficient for Robin") {
		json missing = ValidConfig();
		missing["boundaries"][0]["type"] = "Robin";

		ConfigValidator missing_validator;
		REQUIRE_FALSE(missing_validator.Validate(missing));
		REQUIRE(HasError(missing_validator,
						 "boundaries[0].robin_coefficient"));

		json non_finite = missing;
		non_finite["boundaries"][0]["robin_coefficient"] =
			std::numeric_limits<double>::infinity();

		ConfigValidator finite_validator;
		REQUIRE_FALSE(finite_validator.Validate(non_finite));
		REQUIRE(HasError(finite_validator,
						 "boundaries[0].robin_coefficient"));
	}

	SECTION("rejects a Robin coefficient on non-Robin boundaries") {
		json config = ValidConfig();
		config["boundaries"][0]["robin_coefficient"] = 1.0;

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "boundaries[0].robin_coefficient"));
	}
}

json ValidMqsConfig() {
	json config = ValidConfig();
	config["simulation"]["physics_type"] = "magnetoquasistatics";
	config["materials"][0]["properties"] = {{"mu_r", 1.0}, {"sigma", 0.0}};
	config["scenarios"][0]["frequency"] = 60.0;
	return config;
}

bool HasError(const ConfigValidator& validator, const std::string& field) {
	for (const auto& error : validator.GetErrors()) {
		if (error.field == field) return true;
	}
	return false;
}

} // namespace

TEST_CASE("ConfigValidator validates open-current MQS regions",
		  "[config_validator][mqs][regions]") {
	SECTION("accepts a conductive open-current region") {
		json config = ValidMqsConfig();
		config["materials"][0]["properties"]["sigma"] = 1.0e6;
		config["regions"][0]["current_constraint"] = "open";
		config["terminals"] = json::array();
		config["scenarios"][0]["excitations"] = json::array();

		ConfigValidator validator;
		REQUIRE(validator.Validate(config));
	}

	SECTION("rejects unsupported values and physics") {
		json config = ValidConfig();
		config["regions"][0]["current_constraint"] = "fixed";

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "regions[0].current_constraint"));
	}

	SECTION("rejects zero conductivity") {
		json config = ValidMqsConfig();
		config["regions"][0]["current_constraint"] = "open";
		config["terminals"] = json::array();
		config["scenarios"][0]["excitations"] = json::array();

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "regions[0].current_constraint"));
	}

	SECTION("rejects overlap with a massive terminal") {
		json config = ValidMqsConfig();
		config["materials"][0]["properties"]["sigma"] = 1.0e6;
		config["regions"][0]["current_constraint"] = "open";
		config["terminals"] = json::array({
			{{"name", "Port"}, {"excitation_type", "current"},
			 {"conductor_type", "massive"}, {"entity_group", "Domain"}}
		});
		config["scenarios"][0]["excitations"] = json::array();

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "regions[0].current_constraint"));
	}
}

TEST_CASE("ConfigValidator requires every meshed domain attribute to have a region",
		  "[config_validator][regions]") {
	// Two domain attributes in the mesh; the base config's single region covers
	// only attribute 1. An uncovered attribute gets no material, contributing zero
	// stiffness and making the assembled system singular, so it must be an error.
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
		2, 1, mfem::Element::QUADRILATERAL, true, 2.0, 1.0);
	mesh.SetAttribute(0, 1);
	mesh.SetAttribute(1, 2);
	mesh.SetAttributes();

	SECTION("rejects a domain attribute no region claims") {
		json config = ValidConfig();

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config, &mesh));
		REQUIRE(HasError(validator, "regions"));
	}

	SECTION("accepts full coverage") {
		json config = ValidConfig();
		config["entity_groups"].push_back(
			{{"name", "Domain2"}, {"dim", 2}, {"attribute_ids", {2}}});
		config["regions"].push_back(
			{{"entity_group", "Domain2"}, {"material", "Dielectric"}});

		ConfigValidator validator;
		REQUIRE(validator.Validate(config, &mesh));
	}

	SECTION("is skipped without a mesh") {
		json config = ValidConfig();

		ConfigValidator validator;
		REQUIRE(validator.Validate(config));
	}
}

TEST_CASE("ConfigValidator validates material names",
		  "[config_validator][materials]") {
	SECTION("requires a name") {
		json config = ValidConfig();
		config["materials"][0].erase("name");

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "materials[0].name"));
	}

	SECTION("rejects duplicate names") {
		json config = ValidConfig();
		config["materials"].push_back(
			{{"name", "Dielectric"}, {"properties", {{"epsilon_r", 3.0}}}});

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "materials[1].name"));
	}

	SECTION("rejects a region referencing an undeclared material") {
		json config = ValidConfig();
		config["regions"][0]["material"] = "Missing";

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "regions[0].material"));
	}

	SECTION("rejects a non-string material reference") {
		json config = ValidConfig();
		config["regions"][0]["material"] = 1;

		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "regions[0].material"));
	}
}

TEST_CASE("ConfigValidator accepts the canonical schema", "[config_validator]") {
	ConfigValidator validator;
	REQUIRE(validator.Validate(ValidConfig()));
	REQUIRE(validator.GetErrors().empty());
}

TEST_CASE("ConfigValidator rejects wrong JSON types without throwing", "[config_validator]") {
	json config = ValidConfig();
	config["simulation"]["order"] = "first";
	config["entity_groups"][0]["attribute_ids"][0] = 1.5;
	config["scenarios"][0]["excitations"][0]["floating"] = "false";

	ConfigValidator validator;
	REQUIRE_NOTHROW(validator.Validate(config));
	REQUIRE_FALSE(validator.Validate(config));
	REQUIRE(HasError(validator, "simulation.order"));
	REQUIRE(HasError(validator, "entity_groups[0].attribute_ids[0]"));
	REQUIRE(HasError(validator, "scenarios[0].excitations[0].floating"));
}

TEST_CASE("ConfigValidator enforces canonical simulation field names", "[config_validator]") {
	json config = ValidConfig();
	config["simulation"].erase("physics_type");
	config["simulation"]["physics"] = "electrostatics";
	config["simulation"]["geometry"] = "planar";
	config["simulation"]["results_file"] = "result.msh";

	ConfigValidator validator;
	REQUIRE_FALSE(validator.Validate(config));
	REQUIRE(HasError(validator, "simulation.physics"));
	REQUIRE(HasError(validator, "simulation.geometry"));
	REQUIRE(HasError(validator, "simulation.results_file"));
}

TEST_CASE("ConfigValidator enforces canonical terminal field names", "[config_validator]") {
	json config = ValidConfig();
	config["terminals"][0].erase("excitation_type");
	config["terminals"][0]["excitation"] = "voltage";

	ConfigValidator validator;
	REQUIRE_FALSE(validator.Validate(config));
	REQUIRE(HasError(validator, "terminals[0].excitation"));
	// The legacy spelling gets the actionable rename message on its own; the
	// generic "missing required field" error would only add noise here.
	REQUIRE_FALSE(HasError(validator, "terminals[0].excitation_type"));
}

TEST_CASE("ConfigValidator requires terminal excitation_type", "[config_validator]") {
	json config = ValidConfig();
	config["terminals"][0].erase("excitation_type");

	ConfigValidator validator;
	REQUIRE_FALSE(validator.Validate(config));
	REQUIRE(HasError(validator, "terminals[0].excitation_type"));
}

TEST_CASE("ConfigValidator rejects duplicate names", "[config_validator]") {
	SECTION("entity groups") {
		json config = ValidConfig();
		config["entity_groups"].push_back(
			{{"name", "Domain"}, {"dim", 2}, {"attribute_ids", {2}}});
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "entity_groups[2].name"));
	}

	SECTION("terminals") {
		json config = ValidConfig();
		config["terminals"].push_back(
			{{"name", "Drive"}, {"excitation_type", "voltage"}, {"entity_group", "Boundary"}});
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "terminals[1].name"));
	}

	SECTION("scenarios") {
		json config = ValidConfig();
		config["scenarios"].push_back(
			{{"name", "Driven"}, {"excitations", json::array()}});
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "scenarios[1].name"));
	}
}

TEST_CASE("ConfigValidator checks entity group references and dimensions", "[config_validator]") {
	SECTION("unknown group") {
		json config = ValidConfig();
		config["regions"][0]["entity_group"] = "Missing";
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "regions[0].entity_group"));
	}

	SECTION("wrong group dimension") {
		json config = ValidConfig();
		config["terminals"][0]["entity_group"] = "Domain";
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "terminals[0].entity_group"));
	}

	SECTION("unknown terminal") {
		json config = ValidConfig();
		config["scenarios"][0]["excitations"][0]["terminal"] = "Missing";
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "scenarios[0].excitations[0].terminal"));
	}
}

TEST_CASE("ConfigValidator rejects invalid ranges and enum values", "[config_validator]") {
	json config = ValidConfig();
	config["simulation"]["order"] = 0;
	config["simulation"]["export_refine"] = 0;
	config["simulation"]["amr"] = {
		{"max_iterations", -1},
		{"error_fraction", 1.5},
		{"error_tolerance", -1.0}
	};
	config["entity_groups"][0]["dim"] = 3;
	config["materials"][0]["properties"]["sigma"] = -1.0;
	config["terminals"][0]["conductor_type"] = "unknown";

	ConfigValidator validator;
	REQUIRE_FALSE(validator.Validate(config));
	REQUIRE(HasError(validator, "simulation.order"));
	REQUIRE(HasError(validator, "simulation.export_refine"));
	REQUIRE(HasError(validator, "simulation.amr.max_iterations"));
	REQUIRE(HasError(validator, "simulation.amr.error_fraction"));
	REQUIRE(HasError(validator, "simulation.amr.error_tolerance"));
	REQUIRE(HasError(validator, "entity_groups[0].dim"));
	REQUIRE(HasError(validator, "materials[0].properties.sigma"));
	REQUIRE(HasError(validator, "terminals[0].conductor_type"));
}

TEST_CASE("ConfigValidator enforces MQS scenario frequencies", "[config_validator][mqs]") {
	SECTION("accepts scalar and sweep frequencies") {
		json config = ValidMqsConfig();
		config["scenarios"].push_back({
			{"name", "Sweep"},
			{"frequency", {{"scale", "log"}, {"start", 10.0},
						   {"stop", 1000.0}, {"points", 3}}},
			{"excitations", json::array()}
		});
		ConfigValidator validator;
		REQUIRE(validator.Validate(config));
	}

	SECTION("requires frequency on every scenario") {
		json config = ValidMqsConfig();
		config["scenarios"][0].erase("frequency");
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "scenarios[0].frequency"));
	}

	SECTION("rejects obsolete simulation frequency") {
		json config = ValidMqsConfig();
		config["simulation"]["frequency"] = 60.0;
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "simulation.frequency"));
	}

	SECTION("rejects invalid sweep values") {
		json config = ValidMqsConfig();
		config["scenarios"][0]["frequency"] = {
			{"scale", "octave"}, {"start", -1.0}, {"stop", -2.0}, {"points", 0}
		};
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "scenarios[0].frequency.scale"));
		REQUIRE(HasError(validator, "scenarios[0].frequency.start"));
		REQUIRE(HasError(validator, "scenarios[0].frequency.stop"));
		REQUIRE(HasError(validator, "scenarios[0].frequency.points"));
	}

	SECTION("rejects descending ranges") {
		json config = ValidMqsConfig();
		config["scenarios"][0]["frequency"] = {
			{"scale", "linear"}, {"start", 100.0}, {"stop", 10.0}, {"points", 2}
		};
		ConfigValidator validator;
		REQUIRE_FALSE(validator.Validate(config));
		REQUIRE(HasError(validator, "scenarios[0].frequency.stop"));
	}

	SECTION("does not require frequency for other physics") {
		ConfigValidator validator;
		REQUIRE(validator.Validate(ValidConfig()));
	}
}
