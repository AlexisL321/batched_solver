#ifndef __MG_h__
#define __MG_h__
#include <vector>
#include <cmath>
#include <memory>
#include "util.h"
#include "io.h"
#include <Eigen/Sparse>

using namespace Eigen;
template<typename T> using SpMat = SparseMatrix<T>;
template<typename T> using Vec = Matrix<T, Dynamic, 1>;

template<typename T>
class MG{
	public:
		int lvls;
		int v_cycles;
		int n;
		double tol;
		Vec<T>* b;
		std::shared_ptr<std::vector<SpMat<T>> discretizations;
		std::shared_ptr<std::vector<Vec<T>> RHS;
		std::shared_ptr<std::vector<Vec<T>> X;
		Vec<T>* x_init;
		int blc_size;
		int num_blc;
	public:
		MG(int lvls, int v_cycles, int n, SpMat<T>& A, Vec<T>& b, Vec<T>& x_init){
			this->lvls = lvls;
			this->v_cycles = v_cycles;
			this->n = n;
			discretizations = std::make_shared<SpMat<T>>(lvls);
			RHS = std::make_shared<Vec<T>>(lvls);
			X = std::make_shared<Vec<T>>(lvls);
			init_matrices(A, n, lvls); //initialize discretization matrices for different levels

			this->x_init = &x_init;
			this->b = &b;
			RHS->at(0) = *b;

		}
	private:
		//initialize discretization matrices at all levels by reading in csv files
		void init_matrices(SpMat<T>& A, int n, int lvls){
			//call a function in io.h to read in csv files
		}

	public:
		void solve(int blc_size){
			this->blc_size = blc_size;
			this->num_blc = std::ceil(n/blc_size);
			for (int i = 0; i < v_cycles; i++) {
				for (int j = 0; j < lvls; j++) {
					smooth(j);
					//calculate residual r

					//restrict r and store it as b in RHS at j+1

				}

			}

		}

		void smooth(int j) {
			int divide = std::pow(4, j+1);
			SpMat<T>* A = this->discretizations->at(j);
			Vec<T>* b = this->RHS->at(j);
			Vec<T>* x = this->X->at(j);

			//Jacobi smoother using functions in Jacobi class

			//store the relaxed x at jth in X vector

		}

		void prolongation(){

		}

		void restriction(){

		}
};

#endif
