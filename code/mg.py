import numpy as np
import scipy as sp
import scipy.sparse as sps

import block
import pointwise

class Grid:
	def __init__(self, rhs, x):
		self.f = rhs.copy()
		self.v = x.copy()

def mg(x_init, A, b, smoother_iter, num_level, num_cycle):
	x = x_init
	tol = 0
    if sps.issparse(b):
        tol = 1e-6*sps.linalg.norm(b)
    else:
        tol = 1e-6*np.linalg.norm(b)
	print("tol:", tol)
	res_list = []
	res = 0

	# multigrid
	grids = [Grid(b, x_init) for _ in range(2)]
	A_array = [A]
	for cycle in range(num_cycle):
		# go down in the V-cycle
		for level in range(num_level):
			i1 = level % 2 # the current level (where relaxation is done)
			i2 = (level+1) % 2 # where coarser level grid is stored
			# relaxation on current grid level
			grids[i1].v = relax(A_array[level], grids[i1].v, grids[i1].f, smoother_iter)
			# restriction (current->coarser)
			if cycle == 0:
				A_restrict = restrict(A_array[level])
				A_array.append(A_restrict)
			grids[i2].f = restrict(grids[i1].f - \
				A_array[level+1].dot(grids[i1].v))

		# At the coarsest grid->direct solve
		i1 = num_level % 2 # where the coarsest grid is stored
		grids[i1].v = direct_solve(grids[i1].v, grids[i1].f)

		# go up in the V-cycle
		for level in range(num_level):
			i1 = (num_level-1-level) % 2 # the current level
			i2 = (num_level-level) % 2 # the coarser level
			# prolongation (coarser->current)
			grids[i1].v = grids[i1].v + prolong(grids[i2].v)
			# relax on current grid leve
			grids[i1].v = relax(A_array[num_level-1-level], grids[i1].v, grids[i1].f, smoother_iter)

		x = grids[0].v
		error = b - A.dot(x)

		if sps.issparse(error):
			res = sps.linalg.norm(error)
		else:
			res = np.linalg.norm(error)
		res_list.append(res)

		if res < tol:
			break
	
	return x, res_list, cycle


def relax(A, x, b, num_iter, block=False):
	if !block: # if using pointwise
		x_relaxed,_,_ = pointwise.jacobi(x.copy(), A, b, max_iter=num_iter)
		return x_relaxed
	else:
		bsize = 10
		num_blc = 10
		solver = 'sparse_chol'
		x_relaxed,_,_ = block.blc_jacobi(x.copy(), A, b, bsize, num_blc,
			solver, max_iter = num_iter)
		return x_relaxed

def restrict(x):
	
def prolong(x):
	
