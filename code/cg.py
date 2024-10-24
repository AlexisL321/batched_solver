# Conjugate Gradient Method
import numpy as np
import scipy as sp
import scipy.sparse as sps

import pointwise
import block


# Conjugate Gradient method
def CG(x_init, A, b):
	x0 = x_init.copy()
	r0 = sps.csc_matrix(b - A.dot(x0))
	p0 = r0
	k = 0
	rk = r0
	pk = p0
	xk = x0
	#print(sp.linalg.norm(b).shape)
	tol = 0
	if sps.issparse(b):
		tol = 1e-6*sps.linalg.norm(b)
	else:
		tol = 1e-6*np.linalg.norm(b)
	print(tol)
	res_list = []
	res = 0

	while True:
		alpha_k = rk.T.dot(rk)/pk.T.dot(A.dot(pk))
		xk1 = xk + pk.multiply(alpha_k)
		rk1 = rk - A.dot(pk).multiply(alpha_k)
		res = sps.linalg.norm(rk1)
		res_list.append(res)
		if res < tol:
			break
		beta_k = rk1.T.dot(rk1)/rk.T.dot(rk)
		pk1 = rk1 + pk.multiply(beta_k)

		print(k, ": ",res)
		k = k+1
		xk = xk1
		rk = rk1
		pk = pk1
		
	return xk, res_list, k

def CG_smoother(x_init, A, b, smoother, num_iter):
	x0 = x_init.copy()
	r0 = sps.csc_matrix(b - A.dot(x0))
	p0 = r0
	k = 0
	rk = r0
	rk1_p = rk
	pk = p0
	xk = x0
	#print(sp.linalg.norm(b).shape)
	tol = 0
	if sps.issparse(b):
		tol = 1e-6*sps.linalg.norm(b)
	else:
		tol = 1e-6*np.linalg.norm(b)
	print(tol)
	res_list = []
	res = 0

	while True:
		alpha_k = rk1_p.T.dot(rk)/pk.T.dot(A.dot(pk))
		xk1 = xk + pk.multiply(alpha_k)
		rk1 = rk - A.dot(pk).multiply(alpha_k)
		res = sps.linalg.norm(rk1)
		res_list.append(res)
		if res < tol:
			break

		# preconditioning using pointwise Jacobi
		#b_p = np.zeros(b.shape)
		rk1_p,_,_ = pointwise.pointwise_solver(rk1, A, rk1, smoother,
			max_iter=num_iter)

		beta_k = rk1_p.T.dot(rk1)/rk1_p.T.dot(rk)
		pk1 = rk1_p + pk.multiply(beta_k)

		print(k, ": ",res)
		k = k+1
		xk = xk1
		rk = rk1_p
		pk = pk1
	
	return xk, res_list, k

def CG_smoother_blc(x_init, A, b, bsize, num_blc, solver, num_iter):
	x0 = x_init.copy()
	r0 = sps.csc_matrix(b - A.dot(x0))
	p0 = r0
	k = 0
	rk = r0
	rk1_p = rk
	pk = p0
	xk = x0
	#print(sp.linalg.norm(b).shape)
	tol = 0
	if sps.issparse(b):
		tol = 1e-6*sps.linalg.norm(b)
	else:
		tol = 1e-6*np.linalg.norm(b)
	print(tol)
	res_list = []
	res = 0

	while True:
		alpha_k = rk1_p.T.dot(rk)/pk.T.dot(A.dot(pk))
		xk1 = xk + pk.multiply(alpha_k)
		rk1 = rk - A.dot(pk).multiply(alpha_k)
		res = sps.linalg.norm(rk1)
		res_list.append(res)
		if res < tol:
			break

		# preconditioning using pointwise Jacobi
		#b_p = np.zeros(b.shape)
		rk1_p,_,_ = block.blc_jacobi(rk1, A, rk1, bsize, num_blc, solver,
			max_iter=num_iter)

		beta_k = rk1_p.T.dot(rk1)/rk1_p.T.dot(rk)
		pk1 = rk1_p + pk.multiply(beta_k)

		print(k, ": ",res)
		k = k+1
		xk = xk1
		# why rk =rk1 doesn't converge?
		rk = rk1_p
		pk = pk1
		
	return xk, res_list, k

