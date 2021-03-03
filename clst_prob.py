from ase.io import read,write
from ase import Atoms,Atom, geometry
from spglib import *
from lib_basis import *
from clst_avg import *
import numpy as np
import collections
import itertools

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
		self.mvec =','.join(str(int(i)) for i in np.ones(len(vertices)))
		# possible occupations of clusters
		self.occ = all_clst_labels(spin_variables,len(vertices))
		#initialize dictionary of cluster occupations
		self.occs = { oc : 0  for oc in self.occ}
		return None
	def copy(self):
		return self

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
		#spg_tmp = (vectors, sites, numbers)  #[list(set(numbers))[0]]*len(numbers) )
		#vectors, sites, numbers_tmp = refine_cell(spg_tmp,symprec=1.e-3)
		
		
		self.vectors = vectors
		self.sites = sites
		self.prim = prim
		self.numbers = numbers
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

	def cluster_avg(self,is_slab=False,dist_tol = 0.01):
		#This assumes an (m x n x l) supercell!
		#TODO implement sqrt(m) x sqrt(n) x l (general cells)
		if is_slab:
			sym_ops = self.prim.sym_2d()
		elif not is_slab:
			sym_ops = self.prim.sym

		clusters = self.clusters.copy()
		sc_mat = np.linalg.inv(self.tmatrix)
	
		clst_dat = {i : {
		'phi'  : None,
		'occs' : {} } for i in range(len(clusters)) }
		count = 0
		for cluster in clusters:
			all_phis= {j:[] for j in cluster.occ}
			all_occs= {j:[] for j in cluster.occ}
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
				calced_occs = None
				#raster the cluster over the crystal
				if is_slab:
					calced_occs = raster_2d(self.sites,np.array(cluster.vertices.copy()), cluster.occs.copy(), sc_mat, scpositions, self.var_map, self.numbers)
				else:
					calced_occs = raster_3d(self.sites,np.array(cluster.vertices.copy()), cluster.occs.copy(), sc_mat, scpositions, self.var_map, self.numbers)
				#append cluster probabilities to the all_occs dict	
				for oc in cluster.occ:
					#cluster probability at given occupation
					#TODO replace denominator with nsites for general supercells
					prb = calced_occs[oc] /(  np.prod(np.diag(self.tmatrix)))
					all_occs[oc].append( prb )
					#TODO calculate for all cluster basis functions 
					# in ternary + systems
					ints = [int(i) for i in oc.split(',')]
					mvec = cluster.mvec
					phi = phi_t( [int(i) for i in mvec.split(',')], np.array(ints), len(cluster.variables) )
					all_phis[oc].append(phi)
			phis = [np.average(all_phis[k]) * np.average(all_occs[k]) for k in cluster.occ]
			clst_dat[count]['phi'] = np.average(phis)
			for oc in cluster.occ:
				clst_dat[count]['occs'][oc] = np.average(all_occs[oc])
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
