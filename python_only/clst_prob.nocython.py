from ase.io import read,write
from ase import Atoms,Atom, geometry
import time
from spglib import *
from lib_basis import *
#from clst_avg import *
import numpy as np
import collections
import itertools


#vectorized distance calculation w/ periodic boundaries
def periodic_dist(x0, x1, dimensions=np.array([1.,1.,1.])):
	delta = np.abs(x0 - x1)
	delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
	return np.sqrt((delta ** 2).sum(axis=-1))

def all_clst_labels(variables,clst_size):
	#returns a list of all possible occupation labels of a cluster
	# given the set of possible spin variables and # of vertices
	net = []
	clst_comps = list(itertools.combinations_with_replacement(variables,clst_size))
	for comp in clst_comps:
		permutes =itertools.permutations(comp)
		for perm in permutes:
			p = ','.join(str(k) for k in perm)
			if p not in net:
				net.append(p)
	return net

def get_m_labels(fcn_indices,clst_size):
	#returns m vector (sanchez 1984) for clusters with symmetric
	# occupations. if they are NOT symmetric, use all_clst_labels instead
	net = []
	clst_comps = list(itertools.combinations_with_replacement(fcn_indices,clst_size))
	for comp in clst_comps:
		p = ','.join(str(k) for k in comp)
		net.append(p)
	return net

def get_layer(atoms,atom):
	atomlayers= geometry.get_layers(atoms,(0,0,1),tolerance=0.01)
	return atomlayers[0][atom.index]

def get_top_struct(atoms):
	atoms_lst = []
	layers = geometry.get_layers(atoms,(0,0,1),tolerance=0.1)
	lays = list(set(layers[0]))
	for atom in atoms:
		if get_layer(atoms,atom) in lays[-5:]:
			atoms_lst.append(atom)
	new_atoms = Atoms(atoms_lst)
	new_atoms.set_cell(np.add(atoms.get_cell(), [0,0,0.]) )
	new_atoms.center(axis=2)
	new_atoms.set_pbc([True,True,True])
	return new_atoms

class cluster:
	def __init__(self, vertices,spin_variables, **kwargs):
		# cluster vertices in fractional coords
		self.vertices = vertices
		self.nverts = len(vertices)
		self.variables= spin_variables
		# m vectors from sanchez 1984 - default to 1
		self.mvec =[','.join(str(int(i)) for i in np.ones(len(vertices)))]
		# possible occupations of clusters
		self.occ = all_clst_labels(spin_variables,len(vertices))
		#initialize dictionary of cluster occupations
		self.occs = { oc : 0  for oc in self.occ}
		return None
	def copy(self):
		return self

	def add_mvec(self,vec):
		vecs = self.mvec
		if vec not in vecs:
			vecs.append(vec)
		self.mvec = vecs


class prim_cell:
	def __init__(self, prim_vectors, sites, numbers, zsubs, sublattices):
		# prim_vectors: np.array, cell vectors
		# sites: list of lists or np.ndarray, lattice sites in crystal coordinates for atomic basis
		# nsubs: integer, number of sublattices for the bulk material(s) (e.g. 1 for pure FCC Pt and 2 for NaCl in rocksalt)
		# zsubs: list of lists, layer indices for layers comprising sublattices
		#   (To be used for constructing composite clusters near an interface
		#	-assumes interface is along z-direction !

		spg_cell = (prim_vectors, sites, numbers)
		
		self.prim_vectors = prim_vectors
		self.sites = sites
		self.sym= get_symmetry(spg_cell,symprec=1.e-3)['rotations']

		# TODO implement sublattice detection
		self.nsubs = 1
		self.zsubs = zsubs
		self.sublattices = sublattices
		self.space = get_spacegroup(spg_cell,symprec=1e-5)

		return None
	def sym_2d(self):
		#returns 2D symmetry opperations (no inversion/reflection/rotation 
		#   across c-axis) if your slab norm is oriented along a or b
		#   it is best to reorient along c-axis
		opps = [ opp for opp in self.sym if opp[2][2] == 1]
		return opps

	def set_sym(self,inp):
		self.sym = inp

class cluster_cell:
	def __init__(self, vectors, sites,numbers, prim):
		# vectors: np.array, cell vectors
		# sites: list of lists or np.ndarray, lattice sites in crystal coordinates for atomic basis
		# prim: prim_cell class, primitive cell used to generate supercell

		#attempt to refine cell with SPGLIB
		spg_tmp = (vectors, sites, numbers)  #[list(set(numbers))[0]]*len(numbers) )
		vectors, sites, numbers_tmp = refine_cell(spg_tmp,symprec=1.e-3)
		
		
		self.vectors = vectors
		self.sites = sites
		self.prim = prim
		self.numbers = numbers_tmp
		self.clusters = []

		#mapping of atomic numbers to spin variables
		# will need to be changed to incorporate sublattices
		var_map = {}
		unique_nums = list(set(numbers))
		vrs = get_variables(len(unique_nums))

		for ind, num in enumerate(unique_nums):
			var_map[np.float64(num)] = np.float64(vrs[ind])
		self.var_map = var_map
		# transformation matrix for generating supercell
		p = prim.prim_vectors
		t = np.matmul(np.linalg.inv(p.T), vectors.T)
		# NOTE may need to introduce tolerance factor!
		t = t.astype(int)
		if t[2][2] == 0:
			t[2][2] = 1
		self.tmatrix = t

		return None

	def add_cluster(self, clust):
		if len(clust.variables) ==2:
			self.clusters.append(clust)
		else:
			possible_mvecs = get_m_labels( range(len(clust.variables))[1:], len(clust.vertices) )
			for vec in possible_mvecs:
				clust.add_mvec(vec)
			print (clust)
			self.clusters.append(clust)
	#NOTE for now, crystal structures that require translational operations to check for equivalent clusters are not implemented

	def cluster_avg(self,is_slab=False,dist_tol = 0.001):
		#This assumes an (m x n x l) supercell!
		if is_slab:
			sym_ops = self.prim.sym_2d()
		elif not is_slab:
			sym_ops = self.prim.sym

		clusters = self.clusters.copy()
		sc_mat = np.linalg.inv(self.tmatrix)
	
		count = 0
		for cluster in clusters:
			clst_dat = {
			'phi'  : {m:None for m in cluster.mvec },
			'occs' : {j:None for j in  cluster.occ} } 
			all_phis= {j:[] for j in cluster.occ}
			all_occs= {j:[] for j in cluster.occ}
			calced_occs= { j:0 for j in cluster.occ}
			for opp in sym_ops:
				positions = cluster.vertices.copy()

				# loop over vertices
				tmpsc = []
				for p in positions:
					pt = np.transpose(np.array([p]))
					transformed = np.matmul(opp,pt)
					tmpsc.append(np.transpose(transformed)[0])

				##convert cluster vertices to supercell basis
				scpositions = np.matmul( tmpsc , sc_mat)
				#raster the cluster over the crystal
				if is_slab:
					clst_pos = []
					#raster cluster over the 2d cell
					for i in range(int(1/sc_mat[0][0])):
						for j in range(int(1/sc_mat[1][1])):
							tmp = (i * sc_mat[0]) + scpositions
							tmp = (j * sc_mat[1]) + tmp
							tmp = np.mod(tmp,  1.)
							clst_pos.append(tmp)

					pdists = [ [ periodic_dist(vert,self.sites) for vert in verts] for verts in clst_pos]

					vert_inds = [ [[ i for i ,val in enumerate(v) if val < dist_tol][0] for v in p] for p in pdists]

					vert_occ =[  [self.var_map[self.numbers[ind]] for ind in vind] for vind in vert_inds]
					for v_oc in vert_occ:
						lab= ','.join(str(int(k)) for k in v_oc)
						calced_occs[lab] +=1

				else:

					clst_pos = []
					#raster cluster over the 3d cell
					for i in range(int(1/sc_mat[0][0])):
						for j in range(int(1/sc_mat[1][1])):
							for k in range(int(1/sc_mat[2][2])):
								tmp = (i * sc_mat[0]) + scpositions
								tmp = (j * sc_mat[1]) + tmp
								tmp = (k * sc_mat[2]) + tmp
								tmp = np.mod(tmp,  1.)
								clst_pos.append(tmp)
					pdists = [ [ periodic_dist(vert,self.sites) for vert in verts] for verts in clst_pos]

					vert_inds = [ [[ i for i ,val in enumerate(v) if val < dist_tol][0] for v in p] for p in pdists]

					vert_occ =[  [self.var_map[self.numbers[ind]] for ind in vind] for vind in vert_inds]
					for v_oc in vert_occ:
						lab= ','.join(str(int(k)) for k in v_oc)
						calced_occs[lab] +=1
			#TODO extend cluster averaging beyond integer mxnxl supercells
			for oc in cluster.occ:
				clst_dat['occs'][oc] = np.average(calced_occs[oc]/( len(sym_ops) * np.prod(np.diag(self.tmatrix))))
			for vec in cluster.mvec:
				weighted_phis = []
				for oc in cluster.occ:
					ints = [int(i) for i in oc.split(',')]
					phi = phi_t( [int(i) for i in vec.split(',')], np.array(ints), len(cluster.variables) )
					avg_prb = clst_dat['occs'][oc]
					weighted_phis.append(phi * avg_prb)

				clst_dat['phi'][vec] = np.average(weighted_phis)

			count +=1
		return  clst_dat


def from_ase(atoms,zsubs):
	def get_layer(atoms,atom):
		atomlayers= geometry.get_layers(atoms,(0,0,1),tolerance=0.01)
		return atomlayers[0][atom.index]
	sublattices = {}
	layer_by_atom = geometry.get_layers(atoms,(0,0,1),tolerance=0.01)[0]
	layers = list(set(layer_by_atom))
	if collections.Counter(zsubs[0]) == collections.Counter(layers):
		vectors = np.array([v for v in atoms.get_cell()])
		sites = np.array(atoms.get_scaled_positions())
		numbers = atoms.get_atomic_numbers()
		sublattices[0] = [at.index for at in atoms]

	else:
		vectors = np.array([v for v in atoms.get_cell()])
		sites = np.array(atoms.get_scaled_positions())
		numbers = atoms.get_atomic_numbers()
		for ind, lat in enumerate(zsubs):
			sublat = []
			for at in atoms:
				if get_layer(atoms,at) in sublat:
					sublat.append(at.index)
			sublattices[ind] = sublat
			
	return ( vectors,sites, numbers,sublattices)
