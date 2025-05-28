#ifndef __CG_h__
#define __CG_h__
#include <vector>
#include <cmath>
#include <iostream>
#include <Eigen/Sparse>
#include "util.h"

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
		int num_iter;
		std::shared_ptr<Vec<T>> x;
		std::unique_ptr<Vec<T>> x_init;
		const SpMat<T>* A;
		const Vec<T>* b;
		//Vec<T>* x_init;
		//Vec<T>* x; //stores final result x
		double res = 10000; //stores residual
		std::unique_ptr<std::vector<double>> res_list;//stores list of residuals
		Vec<T> p;
		Vec<T> r;
		T alpha = 0;
		T beta = 0;
	public:
		CG(int n, SpMat<T>& A, Vec<T>& b, Vec<T>& x_init, int max_iter = 1000) {
			this->n = n;
			this->x_init = std::make_unique<Vec<T>>(n);
			this->x = std::make_shared<Vec<T>>(n);
			*(this->x_init) = *x_init;
			*(this->x) = *x_init;
			this->A = &A;
			this->b = &b;
			this->max_iter = max_iter;

			res_list = std::make_unique<std::vector<double>>();
			tol = 1e-6 * this->x_init.norm();
			p = Vec<T>::Zeros(n);
			r = Vec<T>::Zeros(n);
		}
		~CG()= default;

		void Solve(){
			r = *b - *A * (*x_init);
			p = r;
			num_iter = 0;
			while (num_iter < max_iter) {
				Ap = (*A)*p;
				rrT = r.norm();
				alpha = rrT/(p.transpose()*Ap);
				x = x + alpha * p;
				r_ = r - alpha * Ap;
				res = r_.norm();
				res_list.push_back(res);
				if (res < tol) break;
				beta = res/rrT;
				p = r_ + beta*p;
				num_iter ++;
			}
		}
};

#endif
