#include "reader.hpp"
#include <fstream>

std::vector<double> reader(std::filesystem::path file_path, const size_t size)
{
	std::vector<double> data(size);
	std::ifstream fin(file_path);
	fin.read((char*)data.data(), sizeof(double) * size);
}