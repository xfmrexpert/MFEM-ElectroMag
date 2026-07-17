// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

class StatusReporter {
public:
	enum class Verbosity {
		Status = 0,
		Solver = 1,
		Diagnostics = 2
	};

	class Operation {
	public:
		Operation(StatusReporter& reporter, std::string name)
			: reporter_(reporter),
			  name_(std::move(name)),
			  start_(Clock::now()),
			  exception_count_(std::uncaught_exceptions()) {
			reporter_.Write("Starting", name_);
		}

		Operation(const Operation&) = delete;
		Operation& operator=(const Operation&) = delete;
		Operation(Operation&&) = delete;
		Operation& operator=(Operation&&) = delete;

		~Operation() noexcept {
			const double seconds =
				std::chrono::duration<double>(Clock::now() - start_).count();
			reporter_.Write(std::uncaught_exceptions() > exception_count_
								? "Failed"
								: "Completed",
							name_, seconds);
		}

	private:
		using Clock = std::chrono::steady_clock;

		StatusReporter& reporter_;
		std::string name_;
		Clock::time_point start_;
		int exception_count_;
	};

	explicit StatusReporter(std::ostream& output = std::cout)
		: output_(output) {}

	static StatusReporter& Global() {
		static StatusReporter reporter;
		return reporter;
	}

	static Verbosity VerbosityFromInt(int value) {
		switch (value) {
			case 0: return Verbosity::Status;
			case 1: return Verbosity::Solver;
			case 2: return Verbosity::Diagnostics;
			default:
				throw std::invalid_argument("verbosity must be 0, 1, or 2");
		}
	}

	void SetVerbosity(Verbosity verbosity) { verbosity_ = verbosity; }
	Verbosity GetVerbosity() const { return verbosity_; }

	bool SolverOutputEnabled() const {
		return verbosity_ >= Verbosity::Solver;
	}

	bool DiagnosticsEnabled() const {
		return verbosity_ >= Verbosity::Diagnostics;
	}

	int SolverPrintLevel(int configured_level) const {
		return SolverOutputEnabled() ? configured_level : 0;
	}

	Operation Start(std::string name) {
		return Operation(*this, std::move(name));
	}

	void Status(const std::string& message) {
		output_ << message << '\n';
	}

	void Diagnostic(const std::string& message) {
		if (DiagnosticsEnabled()) {
			output_ << message << '\n';
		}
	}

private:
	void Write(const char* state, const std::string& name) noexcept {
		try {
			output_ << state << ' ' << name << "...\n";
		}
		catch (...) {
		}
	}

	void Write(const char* state, const std::string& name, double seconds) noexcept {
		try {
			std::ostringstream elapsed;
			elapsed << std::fixed << std::setprecision(3) << seconds;
			output_ << state << ' ' << name << " in " << elapsed.str() << " s\n";
		}
		catch (...) {
		}
	}

	std::ostream& output_;
	Verbosity verbosity_ = Verbosity::Status;
};
