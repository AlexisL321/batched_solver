# This file is for function to load data in mat format
import scipy as sp
import scipy.io as io
import matplotlib.pyplot as plt
import scipy.sparse as sparse
import numpy as np
import scipy.sparse.linalg as lg

def load_data(path):
	data = io.loadmat(path)['Aglobal']
	return data

def load_all(path1, path2, path3, path4):
	data1 = io.loadmat(path1)['Aglobal']
	data2 = io.loadmat(path2)['Aglobal']
	data3 = io.loadmat(path3)['Aglobal']
	data4 = io.loadmat(path4)['Aglobal']
	#print(data1.keys())
	return data1, data2, data3, data4

def visualize(matrix):
	sp_mat = sparse.csr_matrix(matrix)
	#plt.figure(figsize=(5,5))
	plt.spy(sp_mat, markersize=1)
	plt.title("sparse structure")
	plt.show()

# Does not work, too slow
def condition_num(matrix):
	evalue, evector = lg.eigsh(matrix, which='LM')
	evalue2, evector2 = lg.eigsh(matrix, which='SM')
	evalue = abs(evalue)
	evalue2 = abs(evalue2)
	return evalue.max()/evalue2.min()

def postorder(matrix):
	# try post-ordering
	return 0

def chmk_order(matrix):
	# try Cuthill-McKee ordering
	return 0

def main():
	path1 = '../data/2D_Square_Domain_Structured_131072.mat'
	path2 = '../data/2D_Square_Domain_Structured_524288.mat'
	path3 = '../data/2D_Square_Domain_UnStructured_26508.mat'
	path4 = '../data/2D_Square_Domain_UnStructured_5408.mat'
	mat1, mat2, mat3, mat4 = load_all(path1, path2, path3, path4)
	#visualize(mat1)
	#visualize(mat2)
	#visualize(mat3)
	#visualize(mat4)

if __name__ == "__main__":
	main()
