import numpy as np
cimport numpy as np
import cython

#cython functions to evaluate cluster averages
#  significantly better scaling than computing
#  cluster averages in python alone

DTYPE=np.float64

ctypedef np.float_t DTYPE_t

#@cython.cdivision(True)

#vectorized distance calculation w/ periodic boundaries
def periodic_dist(np.ndarray x0,np.ndarray x1):
	assert x0.dtype == DTYPE and x0.dtype == DTYPE
	cdef np.ndarray dimensions= np.ones((3))
	cdef np.ndarray delta = np.zeros(x0.shape[0])
	cdef np.ndarray dist = np.zeros(x1.shape[0])
	delta = np.abs(x0 - x1)
	delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
	dist = np.sqrt((delta ** 2).sum(axis=-1))
	return dist

def raster_2d(np.ndarray sites,np.ndarray verts, str d_occ, np.ndarray sc_mat,np.ndarray scpositions, dict var_map, np.ndarray numbers):
	assert verts.dtype == DTYPE and sc_mat.dtype == DTYPE and scpositions.dtype==DTYPE and sites.dtype == DTYPE
	cdef float dist_tol =0.001
	cdef float xyz, value
	cdef int i,j,ind,count,matchcount, ci,cp,cxyz,index,vt_count
	cdef np.ndarray p, vts
	cdef np.ndarray vert_occ = np.zeros(verts.shape[0])
	cdef np.ndarray vert_inds = np.zeros(verts.shape[0])
	cdef np.ndarray tmp = np.zeros((3))
	cdef np.ndarray clst_pos = np.zeros((int(1/round(np.prod(np.diag(sc_mat)),8)),verts.shape[0],3))
	cdef np.ndarray vert = np.zeros((3))
	matchcount = 0
	ci=0
	#raster cluster over the 2d cell
	for i in range(int(1/sc_mat[0][0])):
		for j in range(int(1/sc_mat[1][1])):
			tmp = (i * sc_mat[0]) + scpositions
			tmp = (j * sc_mat[1]) + tmp
			tmp = np.mod(tmp,  1.)
			cp =0
			for p in tmp:
				cxyz =0
				for xyz in p:
					clst_pos[ci][cp][cxyz] = p[cxyz]
					cxyz +=1
				cp +=1
			ci+=1

	# goes through and gets occupation of vertices
	for vts in clst_pos:
		vert_inds = np.zeros(verts.shape[0])#,dtype=np.int)
		vert_occ = np.zeros(verts.shape[0],dtype=np.int)
		vt_count=0
		for vert in vts:

			index=0
			for index in range(int(np.shape(sites)[0])):
				
				value = np.linalg.norm( vert - sites[index])
				if round(value,3) < dist_tol:
					vert_inds[vt_count] = int(index)
			vt_count +=1
		count = 0
		for ind in vert_inds:
			vert_occ[count] = var_map[numbers[ind]]
			count += 1
		# TODO: change instead to count all unique cluster labelings
		#print (''.join(str(int(k)) for k in vert_occ),d_occ)
		if ''.join(str(int(k)) for k in vert_occ) == d_occ:
			matchcount +=1
	
	return matchcount/(1/np.prod(np.diag(sc_mat)))

def raster_3d(np.ndarray sites,np.ndarray verts, str d_occ, np.ndarray sc_mat,np.ndarray scpositions, dict var_map, np.ndarray numbers):
	assert verts.dtype == DTYPE and sc_mat.dtype == DTYPE and scpositions.dtype==DTYPE and sites.dtype == DTYPE
	cdef float dist_tol =0.001
	cdef float xyz, value
	cdef int i,j,k,ind,count,matchcount, ci,cp,cxyz,index,vt_count
	cdef np.ndarray p, vts
	cdef np.ndarray vert_occ = np.zeros(verts.shape[0])
	cdef np.ndarray vert_inds = np.zeros(verts.shape[0])
	cdef np.ndarray tmp = np.zeros((3))
	cdef np.ndarray clst_pos = np.zeros((int(1/round(np.prod(np.diag(sc_mat)),8)),verts.shape[0],3))
	cdef np.ndarray vert = np.zeros((3))
	matchcount = 0
	ci=0
	#raster cluster over the 3d cell
	for i in range(int(1/sc_mat[0][0])):
		for j in range(int(1/sc_mat[1][1])):
			for k in range(int(1/sc_mat[2][2])):
				tmp = (i * sc_mat[0]) + scpositions
				tmp = (j * sc_mat[1]) + tmp
				tmp = (k * sc_mat[2]) + tmp
				tmp = np.mod(tmp,  1.)
				cp =0
				for p in tmp:
					cxyz =0
					for xyz in p:
						clst_pos[ci][cp][cxyz] = p[cxyz]
						cxyz +=1
					cp +=1
				ci+=1

	# goes through and gets occupation of vertices
	for vts in clst_pos:
		vert_inds = np.zeros(verts.shape[0])#,dtype=np.int)
		vert_occ = np.zeros(verts.shape[0],dtype=np.int)
		vt_count=0
		for vert in vts:
			index=0
			for index in range(int(np.shape(sites)[0])):
				
				value = np.linalg.norm( vert - sites[index])
				if round(value,3) < dist_tol:
					vert_inds[vt_count] = int(index)
			vt_count +=1
		count = 0
		for ind in vert_inds:
			vert_occ[count] = var_map[numbers[ind]]
			count += 1
		# TODO: change instead to count all unique cluster labelings
		#print (''.join(str(int(k)) for k in vert_occ),d_occ)
		if ''.join(str(int(k)) for k in vert_occ) == d_occ:
			matchcount +=1
	
	return matchcount/(1/np.prod(np.diag(sc_mat)))
