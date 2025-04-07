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

}

//deep copy submatrix from a SpMat type
template<typename T>
SpMat<T> submatrix_cpy(SpMat<T>& A, int row, int col, int row_blc, int
col_blc) {

}


#endif
