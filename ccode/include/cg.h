#ifndef __CG_h__
#define __CG_h__
#include <vector>
#include <cmath>
#include <iostream>
#include <Eigen/Sparse>

using namespace Eigen;
template<typename T> using SpMat = SparseMatrix<T>;
template<typename T> using Trip = Triplet<T>;
template<typename T> using Vec = Matrix<T, Dynamic, 1>;

template<typename T>
class CG {
	private:
	int max_iter = 1000;
	public:
		double tol;
		int n; //length of x
		vector<T> x_init;
	public:
		CG(int n) : n(n), x_init(n, T()) {} 
		~CG()= default;
		void Solve

}

#endif
