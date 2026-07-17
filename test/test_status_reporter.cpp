// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include "../src/status_reporter.hpp"

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

TEST_CASE("StatusReporter emits independently parseable JSON Lines", "[status_reporter]") {
	std::ostringstream output;
	std::ostringstream errors;
	StatusReporter reporter(output, errors);
	reporter.SetFormat(StatusReporter::Format::JsonLines);
	reporter.SetVerbosity(StatusReporter::Verbosity::Solver);

	reporter.Status("ready");
	reporter.Diagnostic("hidden diagnostic");
	reporter.Warning("warning text");
	reporter.Error("error text");
	reporter.SolverOutput() << "iteration 1\nresidual 1e-6\n" << std::flush;
	{
		auto operation = reporter.Start("test operation");
	}

	std::vector<nlohmann::json> events;
	std::istringstream lines(output.str());
	for (std::string line; std::getline(lines, line);) {
		REQUIRE_FALSE(line.empty());
		events.push_back(nlohmann::json::parse(line));
	}

	REQUIRE(events.size() == 7);
	REQUIRE(events[0] == nlohmann::json{{"event", "message"}, {"level", "status"}, {"message", "ready"}});
	REQUIRE(events[1]["level"] == "warning");
	REQUIRE(events[2]["level"] == "error");
	REQUIRE(events[3]["level"] == "solver");
	REQUIRE(events[3]["message"] == "iteration 1");
	REQUIRE(events[4]["message"] == "residual 1e-6");
	REQUIRE(events[5]["event"] == "operation");
	REQUIRE(events[5]["state"] == "started");
	REQUIRE(events[6]["state"] == "completed");
	REQUIRE(events[6]["elapsed_seconds"].is_number());
	REQUIRE(errors.str().empty());
}

TEST_CASE("StatusReporter emits JSON diagnostics only at highest verbosity", "[status_reporter]") {
	std::ostringstream output;
	StatusReporter reporter(output);
	reporter.SetFormat(StatusReporter::Format::JsonLines);

	reporter.Diagnostic("hidden");
	REQUIRE(output.str().empty());

	reporter.SetVerbosity(StatusReporter::Verbosity::Diagnostics);
	reporter.Diagnostic("visible");

	const nlohmann::json event = nlohmann::json::parse(output.str());
	REQUIRE(event["event"] == "message");
	REQUIRE(event["level"] == "diagnostic");
	REQUIRE(event["message"] == "visible");
}
