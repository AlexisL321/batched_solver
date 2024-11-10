import numpy as np
import scipy as sp
import scipy.sparse as sps
import scipy.sparse.linalg as lg
import load_data
import cProfile
import functools
import pstats

def profile(func):
	@functools.wraps(func)
	def wrapper(*args, **kwargs):
		profiler = cProfile.Profile()
		profiler.enable()
		result = func(*args, **kwargs)
		profiler.disable()
		profiler.dump_stats(f"../profile/{func.__name__}.prof")
		p = pstats.Stats(f"../profile/{func.__name__}.prof")
		with open(f"../profile/{func.__name__}.txt", "w") as f:
			p.stream = f
			p.sort_stats("time").print_stats()
		return result
	return wrapper

@profile
def direct_spsolve(A, b, reorder=None):
	return lg.spsolve(A, b, permc_spec=reorder)

def main():
	reordering = input("reordering:")
	if reordering == "n" or reordering == "no":
		reordering = None
	
	path1 = '../data/2D_Square_Domain_Structured_524288.mat'
	mat1 = load_data.load_data(path1)
	num_row = mat1.shape[0]
	b = sps.random(num_row, 1, density=0.01, format='csc')
	x = direct_spsolve(mat1, b, reorder = reordering)
	verify = mat1.dot(x).reshape(-1,1) #(16224,)->(16224,1)
	#print(b.shape, verify.shape)
	error = np.linalg.norm(b-verify)
	return error

if __name__ == "__main__":
	result = main()
	print(result)
