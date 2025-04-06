import scipy as sp
import scipy.sparse as sps
import cProfile
import numpy as np
import scipy.linalg.lapack as lapack
import sksparse.cholmod as cholmod

def blc_jacobi(x_init, A, b, bsize, num_blc, solver, max_iter=1000):
	print("called block jacobi")
	n = A.shape[1]
	if bsize * num_blc < n:
		try:
			raise NameError("wrong block size and block number")
		except NameError:
			print("wrong block size and block numer")
			raise
	res = 1000
	tol = 0
	if sps.issparse(b):
		tol = 1e-6 * sps.linalg.norm(b)
	else:
		tol = 1e-6 * np.linalg.norm(b)
	#print(tol)
	x = x_init.copy()
	num_iter = 0 
	list_res = []
	while res > tol:
		x_init = x.copy()
		for i in range(num_blc):
			update = b[i*bsize:(i+1)*bsize].copy()
			for j in range(num_blc):
				if j != i:
					update -= A[i*bsize:(i+1)*bsize,
						j*bsize:(j+1)*bsize].dot(x_init[j*bsize:\
						(j+1)*bsize])
					#print("update size: ", update.shape)
			update =\
				solve(A[i*bsize:(i+1)*bsize,i*bsize:(i+1)*bsize],
				update, bsize, solver)
			x[i*bsize:(i+1)*bsize] = update.copy()
			#print(x[j*bsize:(j+1)*bsize-1])

		num_iter += 1
		if num_iter > max_iter:
			break
		res = sps.linalg.norm(x-x_init)
		print(num_iter, ": ",res)
		list_res.append(res)
	return x, list_res, num_iter

def blc_GaussSeidel(x_init, A, b, bsize, num_blc, solver, max_iter=1000):
	print("called block Gauss-Seidel")
	n = A.shape[1]
	if bsize * num_blc < n:
		try:
			raise NameError("wrong block size and block number")
		except NameError:
			print("wrong block size and block numer")
			raise
	res = 1000
	tol = 0 
	if sps.issparse(b):
		tol = 1e-6*sps.linalg.norm(b)
	else:
		tol = 1e-6 * np.linalg.norm(b)
	#print(tol)
	x = x_init.copy()
	num_iter = 0 
	list_res = []
	while res > tol:
		x_init = x.copy()
		for i in range(num_blc):
			update = b[i*bsize:(i+1)*bsize].copy()
			for j in range(num_blc):
				if j != i:
					update -= A[i*bsize:(i+1)*bsize,
						j*bsize:(j+1)*bsize].dot(x[j*bsize:\
						(j+1)*bsize])
					#print("x slice size: ", x[j*bsize:(j+1)*bsize].shape)
			update =\
				solve(A[i*bsize:(i+1)*bsize,i*bsize:(i+1)*bsize],
				update, bsize, solver)
			x[i*bsize:(i+1)*bsize] = update.copy()

		num_iter += 1
		if num_iter > max_iter:
			break
		res = sps.linalg.norm(x-x_init)
		print(num_iter, ": ",res)
		list_res.append(res)
	return x, list_res, num_iter

# Sovle for Ax=b using a designated solver
def solve(A, b, size, solver):	
	if solver == 'invert':
		if size > 3:
			raise NameError("can not do direct invert, matrix size bigger",
			" than 3")
		A_inv = inverse(A.copy(), size)
		return sps.csc_matrix(np.dot(A_inv, b))
	if solver == 'dense_chol':
		if sps.issparse(A):
			A = A.toarray()
		return sps.csc_matrix(dense_chol_solve(A.copy(), b.copy(), size))
	if solver == 'sparse_chol':
		if not sps.issparse(A):
			# use csc to be consistent with sksparse.cholmod
			A = sps.csc_matrix(A)
		return sps.csc_matrix(sparse_chol_solve(A.copy(), b.copy(), size))


# hard code inverse for small matrix size
def inverse(A, size):
	if size == 3:
		if A.shape[0] != 3:
			raise NameError("matrix size is not 3.")
		return_mat = np.zeros((3,3))
		a1 = A[0,0]
		a2 = A[0,1]
		a3 = A[0,2]
		b1 = A[1,0]
		b2 = A[1,1]
		b3 = A[1,2]
		c1 = A[2,0]
		c2 = A[2,1]
		c3 = A[2,2]
		det = a1*b2*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1 - a2*b1*c3 -\
			b3*c2*a1
		return_mat[0,0] = (b2*c3 - b3*c2)/det
		return_mat[1,0] = (b1*c3 - b3*c1)*(-1)/det
		return_mat[2,0] = (b1*c2 - b2*c1)/det
		return_mat[0,1] = (a2*c3 - a3*c2)*(-1)/det
		return_mat[1,1] = (a1*c3 - a3*c1)/det
		return_mat[2,1] = (a1*c2 - a2*c1)*(-1)/det
		return_mat[0,2] = (a2*b3 - a3*b2)/det
		return_mat[1,2] = (a1*b3 - a3*b1)*(-1)/det
		return_mat[2,2] = (a1*b2 - a2*b1)/det
		return return_mat
	if size == 2:
		if A.shape[0] != 2:
			raise NameError("matrix size is not 2.")
		return_mat = np.zeros((2,2))
		a = A[0,0]
		b = A[0,1]
		c = A[1,0]
		d = A[1,1]
		det = a*d - b*c
		return_mat[0,0] = d/det
		return_mat[0,1] = -b/det
		return_mat[1,0] = -c/det
		return_mat[1,1] = a/det
		return return_mat


def dense_chol_solve(A, b, size):
	if sps.issparse(A):
		A = A.toarray()
	L = sp.linalg.cholesky(A, lower=True)
	x = sp.linalg.cho_solve((L, True), b)
	#print(b.shape)
	#print("in chol_solve: ", x.shape)
	return x

def sparse_chol_solve(A, b, size):
	A = sps.csc_matrix(A)
	factor = cholmod.cholesky(A)
	x = factor(b)
	return x

def blc_solver(x_init, A, b, bsize, num_blc, solver, blc_solver, max_iter=1000):
    if solver == 'jacobi' or solver == 'Jacobi' or solver == 'j':
        return blc_jacobi(x_init, A, b, bsize, num_blc, blc_solver, max_iter=max_iter)
    if solver == 'GaussSeidel' or solver == 'gaussseidel' or solver == \
    'gs':
        return blc_GaussSeidel(x_init, A, b, bsize, num_blc, blc_solver, max_iter=max_iter)
