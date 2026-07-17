// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <cctype>
#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <utility>
#include "json.hpp"

class StatusReporter {
public:
	enum class Format {
		Text,
		JsonLines
	};

	enum class Verbosity {
		Status = 0,
		Solver = 1,
		Diagnostics = 2
	};

private:
	class MessageBuffer : public std::streambuf {
	public:
		explicit MessageBuffer(StatusReporter& reporter) : reporter_(reporter) {}

	protected:
		int overflow(int ch) override {
			if (ch == traits_type::eof()) {
				return traits_type::not_eof(ch);
			}
			Append(static_cast<char>(ch));
			return ch;
		}

		std::streamsize xsputn(const char* data, std::streamsize size) override {
			for (std::streamsize i = 0; i < size; ++i) {
				Append(data[i]);
			}
			return size;
		}

		int sync() override {
			FlushLine();
			return 0;
		}

	private:
		void Append(char ch) {
			if (ch == '\n') {
				FlushLine();
			}
			else if (ch != '\r') {
				line_ += ch;
			}
		}

		void FlushLine() {
			if (!line_.empty()) {
				reporter_.SolverMessage(line_);
				line_.clear();
			}
		}

		StatusReporter& reporter_;
		std::string line_;
	};

public:

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

	explicit StatusReporter(std::ostream& output = std::cout,
		std::ostream& error_output = std::cerr)
		: output_(output),
		  error_output_(error_output),
		  solver_buffer_(*this),
		  solver_output_(&solver_buffer_) {}

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
	void SetFormat(Format format) { format_ = format; }
	Format GetFormat() const { return format_; }
	bool IsMachineReadable() const { return format_ == Format::JsonLines; }

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
		WriteMessage("status", message);
	}

	void Diagnostic(const std::string& message) {
		if (DiagnosticsEnabled()) {
			WriteMessage("diagnostic", message);
		}
	}

	void Warning(const std::string& message) {
		WriteMessage("warning", message, true);
	}

	void Error(const std::string& message) {
		WriteMessage("error", message, true);
	}

	void SolverMessage(const std::string& message) {
		WriteMessage("solver", message);
	}

	std::ostream& SolverOutput() { return solver_output_; }

private:
	void WriteJson(const nlohmann::json& event) {
		output_ << event.dump() << '\n' << std::flush;
	}

	void WriteMessage(const char* level, const std::string& message,
		bool use_error_output = false) {
		if (IsMachineReadable()) {
			WriteJson({
				{"event", "message"},
				{"level", level},
				{"message", message}
			});
		}
		else {
			std::ostream& stream = use_error_output ? error_output_ : output_;
			stream << message << '\n' << std::flush;
		}
	}

	void Write(const char* state, const std::string& name) noexcept {
		try {
			if (IsMachineReadable()) {
				WriteJson({
					{"event", "operation"},
					{"state", state == std::string("Starting") ? "started" : state},
					{"name", name}
				});
			}
			else {
				output_ << state << ' ' << name << "...\n" << std::flush;
			}
		}
		catch (...) {
		}
	}

	void Write(const char* state, const std::string& name, double seconds) noexcept {
		try {
			if (IsMachineReadable()) {
				std::string machine_state = state;
				machine_state[0] = static_cast<char>(std::tolower(
					static_cast<unsigned char>(machine_state[0])));
				WriteJson({
					{"event", "operation"},
					{"state", machine_state},
					{"name", name},
					{"elapsed_seconds", seconds}
				});
			}
			else {
				std::ostringstream elapsed;
				elapsed << std::fixed << std::setprecision(3) << seconds;
				output_ << state << ' ' << name << " in " << elapsed.str() << " s\n"
					<< std::flush;
			}
		}
		catch (...) {
		}
	}

	std::ostream& output_;
	std::ostream& error_output_;
	Verbosity verbosity_ = Verbosity::Status;
	Format format_ = Format::Text;
	MessageBuffer solver_buffer_;
	std::ostream solver_output_;
};
