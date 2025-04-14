#ifndef __JACOBI_h__
#define __JACOBI_h__
#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>
#include <Eigen/Cholesky>
#include "util.h"
//#include <Accelerate/Accelerate.h>


using namespace Eigen;
template<typename T> using SpMat = SparseMatrix<T>;
template<typename T> using Trip = Triplet<T>;
template<typename T> using Vec = Matrix<T, Dynamic, 1>;

template<typename T>
class BlcJacobi {
	private:
		int max_iter = 1000;
	public:
		int n; //matrix size
		int blc_size;
		int num_blc;
		int num_iter;
		double tol;
		//TODO
		std::shared_ptr<Vec<T>> x;
		std::unique_ptr<Vec<T>> x_init;
		const SpMat<T>* A;
		const Vec<T>* b;
		//Vec<T>* x_init;
		//Vec<T>* x; //stores final result x
		double res = 10000; //stores residual
							//Matrix<double, Dynamic, 1> res_list; 
		std::unique_ptr<std::vector<double>> res_list;//stores list of residuals
		std::string ordering;

	public:
		BlcJacobi(int n, int blc_size, int num_blc, Vec<T>& x_init,
				SpMat<T>& A, Vec<T>& b, int max_iter=1000, std::string
				ordering="AMD") {
			if (blc_size * num_blc != n) {
				std::cout<<"block size*number of blocks is not matrix size"<<std::endl;
				exit(1);
			}
			num_iter = 0;
			this->n = n;
			this->blc_size = blc_size;
			this->num_blc = num_blc;
			this->max_iter = max_iter;
			res_list = std::make_unique<std::vector<double>>();

			//TODO
			this->x_init = std::make_unique<Vec<T>>(n);
			this->x = std::make_shared<Vec<T>>(n);
			//this->A = std::make_shared<SpMat<T>>(n, n);
			//this->b = std::make_shared<Vec<T>>(n);
			this->A = A;
			this->b = b;
			//this->x = new Vec<T>(n);
			//this->x_init = new Vec<T>(n);
			*(this->x_init) = *x_init;
			*(this->x) = *x_init;
			tol = 1e-6 * this->x_init.norm();
			this->ordering = ordering;
		}
		~BlcJacobi(){
			//delete x_init;
			//delete x;
		}
		void Jacobi(){
			while (res > tol && num_iter < max_iter) {
				for (int i = 0; i < num_blc; i++) {
					Vec<T> update(blc_size);
					update = b->block<blc_size, 1>(i*blc_size, 0);
					for (int j = 0; j < num_blc; j++) {
						if (i != j) {
							update -= A->block<blc_size, blc_size>(i*blc_size, 
									j*blc_size) * x->block<blc_size, 1>(j*blc_size, 0);

						}
					}
					Solve_s(i, update);

					res = (*x-*x_init).norm();
					*x_init = *x;
					res_list->push_back(res);
				}
			}
		}
		void Solve_s(int i, Vec<T>& RHS){//sparse cholesky
			SimplicialLLT<SpMat<T>> solver;
			SpMat<T> block_mat = A->block(i*blc_size, i*blc_size, blc_size,
					blc_size).eval();
			Vec<T> x_blc(blc_size);

			if (ordering == "N") {//natural ordering
				NaturalOrdering<int> natural;
				solver.analyzePattern(block_mat, natural);
			}
			else if (ordering == "AMD") {//AMD ordering
				AMDOrdering<int> amd;
				solver.analyzePattern(block_mat, amd);
			}
			if (solver.info() != Success) {//Success is an enum in Core module
				std::cout<<"solver.analyzePattern() failed."<<std::endl;
				exit(0);
			}
			solver.factorize(block_mat);
			if (solver.info() != Success) {//Success is an enum in Core module
				std::cout<<"solver.factorize(block_mat) failed."<<std::endl;
				exit(0);
			}
			x_blc = sovler.sovle(RHS);
			//solver.compute(block_mat);
			if (solver.info() != Success) {
				std::cout<<"solver.solve(RHS) failed."<<std::endl;
				exit(0);
			}
			//update this->x with x_blc
			x->segment(i*blc_size, blc_size) = x_blc;
		}

		//dense cholesky
		void Solve_d(int i, Vec<T>& RHS){
			Matrix<T, blc_size, blc_size> blc_mat = A->block(i*blc_size,
					i*blc_size, blc_size, blc_size).eval().toDense();
			LLT<Matrix<T, blc_size, blc_size>> solver(blc_size);
			solver.compute(blc_mat);
			if (solver.info() != Success) {//Success is an enum in Core module
				std::cout<<"dense solver.compute(block_mat) failed."<<std::endl;
				exit(0);
			}
			Vec<T> x_blc(blc_size);
			x_blc = sovler.solve(RHS);
			if (solver.info() != Success) {//Success is an enum in Core module
				std::cout<<"dense solver.solve(block_mat) failed."<<std::endl;
				exit(0);
			}
			x->segment(i*blc_size, blc_size) = x_blc;
		}

#endif
