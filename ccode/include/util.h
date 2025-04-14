#ifndef __UTIL_h__
#define __UTIL_h__
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <memory>
#include <cmath>
#include <vector>

using namespace Eigen;
template<typename T> using SpMat = SparseMatrix<T>;
template<typename T> using Trip = Triplet<T>;
template<typename T> using Vec = Matrix<T, Dynamic, 1>;

//print SpMat
template<typename T>
void print_matrix(SpMat<T>& A, int n) {
	int num_col = A.outerSize();
	for (int i = 0; i < num_col; i++) {
		for (typename SpMat<T>::InnerIterator iter(A, i); iter; iter++) {
			std::cout<<"("<<iter.row()<<", "<<iter.col()<<")"
				<<iter.value()<<std::endl;
		}
	}
}

//deep copy submatrix from a SpMat type
template<typename T>
SpMat<T> submatrix_cpy(SpMat<T>& A, int row, int col, int row_blc, int
		col_blc) {
	std::vector<Trip<T>> trip_list;
	trip_list.reserve(std::round(row_blc/row)*A.nonZeros());

	for (int j = col; j < col+col_blc; j++) {
		for (typename SpMat<T>::InnerIterator iter(A,j); iter.row() < row+row_blc;
				iter++) {
			trip_list.push_back(Trip<T>(iter.row(), iter.col(),
						iter.value()));
			SpMat<T> copy(row_blc, col_blc);
			copy.setFromTriplets(trip_list.begin(), trip_list.end());
			return copy;
		}
	}
}


#endif
