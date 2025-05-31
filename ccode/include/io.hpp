//Assume the sparse matrix came in triplet form in a csv file
#ifndef __IO_h__
#define __IO_h__
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <fstream>
#include <sstream>

using namespace Eigen;

template <typename T>
T parse_value(const std::string& str) {
	std::istringstream iss(str);
	T val;
	iss >> val;
	if (iss.fail()) {
		throw std::runtime_error("failed to parse value: " + str);
	}
	return val;
}

template <typename T>
SparseMatrix<T> sparse_csv(const std::string& filename) {
	std::vector<Triplet<T>> triplets;
	std::ifstream file(filename);
	if(!file.is_open()) {
		std::cout << "Failed to open " << filename << std::endl;
		exit(0);
	}
	std::string line;
	int max_row = 0, max_col = 0;
	std::string srow, scol, sval;
	int row, col;
	T val;
	while(std::getline(file, line)) {
		std::getline(line, srow, ',');
		std::getline(line, scol, ',');
		std::getline(line, sval, ',');

		int row = parse_value<int>(srow);
		int col = parse_value<int>(scol);
		T val = parse_value<T>(sval);

		triplets.emplace_back(row, col, val);
		if (row > max_row) max_row = row;
		if (col > max_col) max_col = col;
	}
	SparseMatrix<T> mat(max_row+1, max_col+1);
	mat.setFromTriplets(triplets.begin(), triplets.end());
	return mat;
}
/*
   int main() {
   SparseMatrix<double> sparseMat = loadCSVasSparse("data.csv");

   std::cout << "Sparse matrix size: " << sparseMat.rows()
   << "x" << sparseMat.cols() << std::endl;
   std::cout << "Non-zero entries: " << sparseMat.nonZeros() << std::endl;

   return 0;
   }
 */

#endif
