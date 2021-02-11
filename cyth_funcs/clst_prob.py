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
			p = ''.join(str(k) for k in perm)
			if p not in net:
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
	def __init__(self, vertices,d_occ, **kwargs):
		# spin_variables : np.array of possible spin variables for sites in a cluster
		# cluster vertices in fractional coords (np.array)
		self.vertices = vertices
		try:
			# desired cluster occupation if applicable ('str')
			self.d_occ = d_occ 
		except KeyError:
			pass

		if kwargs['spin_variables']!= None:
			self.occ = all_clst_labels(kwargs['spin_variables'],len(vertices))
	
		return None

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

	def add_cluster(self, cluster):
		self.clusters.append(cluster)

	#NOTE for both averaging fcns, I need to add a check to see if symmetry permuted clusters are already accounted for by another operation(s) to get the true multiplcity.
	# E.g. for now, crystal structures that require translational operations to check for equivalent clusters are not implemented

	def cluster_avg_2d(self,dist_tol = 0.01):
		#NOTE this assumes an (m x n x 1) supercell!
		clusters = self.clusters
		sc_mat = np.linalg.inv(self.tmatrix)
	
		all_pos = []
		all_match = []
		all_phi = []
		all_sub = []
		for cluster in clusters:
			count = 0
			for opp in self.prim.sym_2d():
				clst_pos = []
				clst_phi = []
				clst_sub = []
				positions = cluster.vertices.copy()

				# loop over vertices
				tmpsc = []
				for p in positions:
					pt = np.transpose(np.array([p]))
					transformed = np.matmul(opp,pt)
					tmpsc.append(np.transpose(transformed)[0])#.T)

				##convert cluster vertices to supercell basis
				scpositions = np.matmul( tmpsc , sc_mat)
				#raster the cluster over the surface
				mtch = raster_2d(self.sites,np.array(cluster.vertices.copy()), cluster.d_occ, sc_mat, scpositions, self.var_map, self.numbers)
				all_match.append(mtch)
				#correlation function disabled
				all_phi.append(0.)

		return  (all_match , all_phi)


	def cluster_avg_3d(self,dof,dist_tol = 0.01):
		#NOTE this assumes an (m x n x l) supercell!
		clusters = self.clusters
		sc_mat = np.linalg.inv(self.tmatrix)
	
		all_pos = []
		all_match = []
		all_phi = []
		for cluster in clusters:
			count = 0
			# now use all operations
			for opp in self.prim.sym:
				clst_pos = []
				clst_phi = []
				matchcount = 0
				positions = cluster.vertices.copy()

				# loop over vertices
				tmpsc = []
				for p in positions:
					pt = np.transpose(np.array([p]))
					transformed = np.matmul(opp,pt)
					tmpsc.append(np.transpose(transformed)[0])

				
			#   ##convert cluster vertices to supercell basis
				scpositions = np.matmul( tmpsc , sc_mat)

				#raster cluster over the 3d cell
				mtch = raster_3d(self.sites,np.array(cluster.vertices.copy()), cluster.d_occ, sc_mat, scpositions, self.var_map, self.numbers)
				#TODO implement correlation functions
				clst_phi.append(0)
				all_match.append(mtch)
				all_phi.append(np.average(clst_phi))

		return  (all_match , all_phi )
			

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
		for ind, lat in enumerate(zsubs):
			sublat = []
			for at in atoms:
				if get_layer(atoms,at) in sublat:
					sublat.append(at.index)
			sublattices[ind] = sublat
	# returns lattice vectors, crystal positions, atomic numbers, and sublattice 
	#    indices for an ASE atoms object	
	return ( vectors,sites, numbers,sublattices)
