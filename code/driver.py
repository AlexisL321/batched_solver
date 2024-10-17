import pointwise
import load_data
import block
import scipy as sp
import scipy.sparse as sps
import matplotlib.pyplot as plt
import cProfile
import numpy as np


def main():
	# load matrices
	path1 = '../data/2D_Square_Domain_Structured_131072.mat'
	path2 = '../data/2D_Square_Domain_Structured_524288.mat'
	path3 = '../data/2D_Square_Domain_UnStructured_26508.mat'
	path4 = '../data/2D_Square_Domain_UnStructured_5408.mat'
	#mat1, mat2, mat3, mat4 = load_data.load_all(path1, path2, path3, path4)
	#num_row = mat1.shape[0]
	#num_col = mat1.shape[1]

	# randomly generate a sparse array b
	#b = sps.random(num_row, 1, density=0.01, format='csr')

	# initial guess of x
	#x_init = sps.random(num_col, 1, density=0.01, format='csr')

	# test with a smaller matrix
	n_test = 1000
	mat1 = sps.random(n_test, n_test, density=0.3, format='csc')
	for i in range(n_test):
		mat1[i,i] += n_test*0.2
	#print(mat1)
	x_init = sps.random(n_test, 1, density=0.3, format='csc')
	x_final = np.random.rand(n_test, 1)
	for i in range (n_test):
		if x_final[i] != 0:
			x_final[i] += 0.1
	b = mat1.dot(x_final)

	# solve for x
	#list_res_j, num_iter_j = pointwise.jacobi(x_init, mat1, b)
	#list_res_g, num_iter_g = pointwise.GaussSeidel(x_init, mat1, b)

	# try block algorithms
	bsize = 10
	num_blc = 100
	solver = 'sparse_chol'
	list_res_j, num_iter_j = block.blc_jacobi(x_init, mat1, b, bsize,
		num_blc, solver)
	list_res_g, num_iter_g = block.blc_GaussSeidel(x_init, mat1, b, bsize,
	num_blc, solver)

	# plot the res and iteration
	plt.plot(list_res_j, label='jacobi', color='red')
	plt.plot(list_res_g, label='GS', color='blue')

	plt.ylabel('residual')
	plt.xlabel('iteration')
	#plt.xlim(0, 1000)
	plt.legend()
	plt.show()


if __name__=="__main__":
	main()
