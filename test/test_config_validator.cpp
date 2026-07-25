// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include "../src/config_validator.hpp"

namespace {

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
			{{"entity_group", "Domain"}, {"material", 1}}
		})},
		{"materials", json::array({
			{{"properties", {{"epsilon_r", 1.0}}}}
		})},
		{"terminals", json::array({
			{{"name", "Drive"}, {"excitation", "voltage"}, {"entity_group", "Boundary"}}
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

bool HasError(const ConfigValidator& validator, const std::string& field) {
	for (const auto& error : validator.GetErrors()) {
		if (error.field == field) return true;
	}
	return false;
}

} // namespace

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
			{{"name", "Drive"}, {"excitation", "voltage"}, {"entity_group", "Boundary"}});
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
	config["simulation"]["frequency"] = -60.0;
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
	REQUIRE(HasError(validator, "simulation.frequency"));
	REQUIRE(HasError(validator, "simulation.export_refine"));
	REQUIRE(HasError(validator, "simulation.amr.max_iterations"));
	REQUIRE(HasError(validator, "simulation.amr.error_fraction"));
	REQUIRE(HasError(validator, "simulation.amr.error_tolerance"));
	REQUIRE(HasError(validator, "entity_groups[0].dim"));
	REQUIRE(HasError(validator, "materials[0].properties.sigma"));
	REQUIRE(HasError(validator, "terminals[0].conductor_type"));
}
