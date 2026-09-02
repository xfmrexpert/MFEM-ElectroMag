#include "mfem.hpp"
#include <iostream>

int main(int argc, char** argv) {
	mfem::Mesh mesh(argv[1], 1, 1);
	std::cout << "dim              = " << mesh.Dimension() << "\n";
	std::cout << "NE               = " << mesh.GetNE() << "\n";
	std::cout << "NBE              = " << mesh.GetNBE() << "\n";
	std::cout << "attributes.Max   = "
			  << (mesh.attributes.Size() ? mesh.attributes.Max() : 0) << "\n";
	std::cout << "bdr_attr.Size    = " << mesh.bdr_attributes.Size() << "\n";
	std::cout << "bdr_attr.Max     = "
			  << (mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0) << "\n";

	auto names = mesh.attribute_sets.GetAttributeSetNames();
	std::cout << "domain set names = " << names.size() << "\n";
	for (const auto& n : names) std::cout << "   '" << n << "'\n";
	auto bnames = mesh.bdr_attribute_sets.GetAttributeSetNames();
	std::cout << "bdr set names    = " << bnames.size() << "\n";
	for (const auto& n : bnames) std::cout << "   '" << n << "'\n";
	return 0;
}
