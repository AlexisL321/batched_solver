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


//deep copy submatrix from a SpMat type



#endif
