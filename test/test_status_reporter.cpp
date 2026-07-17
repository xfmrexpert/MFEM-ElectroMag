// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include "../src/status_reporter.hpp"

#include <sstream>
#include <stdexcept>
#include <string>

TEST_CASE("StatusReporter validates verbosity values", "[status_reporter]") {
	REQUIRE(StatusReporter::VerbosityFromInt(0) == StatusReporter::Verbosity::Status);
	REQUIRE(StatusReporter::VerbosityFromInt(1) == StatusReporter::Verbosity::Solver);
	REQUIRE(StatusReporter::VerbosityFromInt(2) == StatusReporter::Verbosity::Diagnostics);
	REQUIRE_THROWS_AS(StatusReporter::VerbosityFromInt(-1), std::invalid_argument);
	REQUIRE_THROWS_AS(StatusReporter::VerbosityFromInt(3), std::invalid_argument);
}

TEST_CASE("StatusReporter gates solver and diagnostic output", "[status_reporter]") {
	std::ostringstream output;
	StatusReporter reporter(output);

	reporter.Diagnostic("hidden diagnostic");
	REQUIRE_FALSE(reporter.SolverOutputEnabled());
	REQUIRE(reporter.SolverPrintLevel(2) == 0);
	REQUIRE(output.str().empty());

	reporter.SetVerbosity(StatusReporter::Verbosity::Solver);
	reporter.Diagnostic("still hidden");
	REQUIRE(reporter.SolverOutputEnabled());
	REQUIRE(reporter.SolverPrintLevel(2) == 2);
	REQUIRE(output.str().empty());

	reporter.SetVerbosity(StatusReporter::Verbosity::Diagnostics);
	reporter.Diagnostic("visible diagnostic");
	REQUIRE(output.str() == "visible diagnostic\n");
}

TEST_CASE("StatusReporter times completed and failed operations", "[status_reporter]") {
	std::ostringstream output;
	StatusReporter reporter(output);

	{
		auto operation = reporter.Start("successful step");
	}

	try {
		auto operation = reporter.Start("failing step");
		throw std::runtime_error("failure");
	}
	catch (const std::runtime_error&) {
	}

	const std::string text = output.str();
	REQUIRE(text.find("Starting successful step...\n") != std::string::npos);
	REQUIRE(text.find("Completed successful step in ") != std::string::npos);
	REQUIRE(text.find("Starting failing step...\n") != std::string::npos);
	REQUIRE(text.find("Failed failing step in ") != std::string::npos);
}
