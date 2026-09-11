// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>
#include <highfive/H5File.hpp>
#include "mfem.hpp"

namespace matrix_io {

class CouplingMatrixWriter {
public:
	CouplingMatrixWriter(const std::filesystem::path& path,
						 const std::vector<std::string>& terminals,
						 const std::string& physics,
						 const std::string& geometry)
		: file_(path.string(), HighFive::File::Truncate),
		  group_(file_.createGroup("/coupling")),
		  terminal_count_(terminals.size()) {
		file_.createAttribute("schema_version", 1);
		group_.createAttribute("physics_type", physics);
		group_.createAttribute("geometry_type", geometry);
		group_.createDataSet("terminal_names", terminals);
	}

	void WriteMatrix(const std::string& name, const mfem::DenseMatrix& matrix,
					 const std::string& units) {
		std::vector<double> values;
		AppendMatrix(matrix, values);
		WriteValues(name, {terminal_count_, terminal_count_}, values, units);
	}

	void WriteFrequencies(const std::vector<double>& frequencies) {
		MFEM_VERIFY(!frequencies.empty(), "Coupling sweep must contain a frequency.");
		double previous_frequency = 0.0;
		for (double frequency : frequencies) {
			MFEM_VERIFY(std::isfinite(frequency) && frequency > previous_frequency,
				"Coupling frequencies must be positive, unique, and ascending.");
			previous_frequency = frequency;
		}
		auto frequency_dataset = group_.createDataSet("frequency_hz", frequencies);
		frequency_dataset.createAttribute("units", std::string("Hz"));
	}

	void WriteMatrixSeries(const std::string& name,
						   const std::vector<const mfem::DenseMatrix*>& matrices,
						   const std::string& units) {
		MFEM_VERIFY(group_.exist("frequency_hz") &&
			group_.getDataSet("frequency_hz").getElementCount() == matrices.size(),
			"Matrix series must match the frequency axis.");
		std::vector<double> values;
		for (const auto* matrix : matrices) {
			MFEM_VERIFY(matrix, "Matrix series must not contain null matrices.");
			AppendMatrix(*matrix, values);
		}
		WriteValues(name, {matrices.size(), terminal_count_, terminal_count_}, values, units);
	}

private:
	HighFive::File file_;
	HighFive::Group group_;
	std::size_t terminal_count_;

	void AppendMatrix(const mfem::DenseMatrix& matrix,
					  std::vector<double>& values) const {
		MFEM_VERIFY(matrix.Height() == static_cast<int>(terminal_count_) &&
			matrix.Width() == static_cast<int>(terminal_count_),
			"Coupling matrix dimensions must match the terminal names.");
		// MFEM stores columns contiguously; HDF5/C# use the last axis fastest.
		for (int row = 0; row < matrix.Height(); ++row) {
			for (int column = 0; column < matrix.Width(); ++column) {
				values.push_back(static_cast<double>(matrix(row, column)));
			}
		}
	}

	void WriteValues(const std::string& name,
					 const std::vector<std::size_t>& dimensions,
					 const std::vector<double>& values,
					 const std::string& units) {
		auto quantity_group = group_.createGroup(name);
		auto dataset = quantity_group.createDataSet<double>("values", HighFive::DataSpace(dimensions));
		dataset.write_raw(values.data());
		dataset.createAttribute("units", units);
	}
};

} // namespace matrix_io
