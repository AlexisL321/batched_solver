# This file is for pointwise Jacobi method
import numpy as np
import scipy as sp
import scipy.sparse as sps
import scipy.sparse.linalg as lg

def jacobi (x_init, A, b):
	print("called jacobi")
	n = A.shape[1]
	res = 1000
	tol = 1e-6 * np.linalg.norm(b)
	print(tol)
	x = x_init.copy()
	num_iter = 0
	list_res = []
	while res > tol:
		x_init = x.copy()
		for i in range(n):
			update = b[i]
			for j in range(n):
				if j != i and A[i,j] != 0 and x_init[j] != 0:
					update -= A[i, j]*x_init[j]
				#print("update:", update)
			x[i] = (1/A[i,i]) * update

		num_iter += 1
		if num_iter > 1000:
			break
		res = lg.norm(x-x_init)
		print(num_iter, ": ",res)
		#print(x-x_init)
		list_res.append(res)
	
	return list_res, num_iter

def GaussSeidel (x_init, A, b):
	print("called gauss seidel")
	n = A.shape[1]
	res = 1000
	tol = 1e-6 * np.linalg.norm(b)
	print(tol)
	x = x_init.copy()
	num_iter = 0
	list_res = []
	while res > tol:
		x_init = x.copy()
		for i in range(n):
			update = b[i]
			for j in range(n):
				if j != i and A[i,j] != 0 and x_init[j] != 0:
					update -= A[i, j]*x[j]
			x[i] = (1/A[i,i]) * update

		num_iter += 1
		if num_iter > 1000:
			break
		res = lg.norm(x-x_init)
		print(num_iter, ": ", res)
		list_res.append(res)
	
	return list_res, num_iter

