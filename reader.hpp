#pragma once
#include <vector>
#include <filesystem>

std::vector<double> reader(std::filesystem::path file_path, const size_t size);