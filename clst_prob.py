from ase.io import read,write
from ase import Atoms,Atom, geometry
from spglib import *
from lib_basis import *
import numpy as np
import collections
import sys


def get_layer(atoms,atom):
	atomlayers= geometry.get_layers(atoms,(0,0,1),tolerance=0.01)
	return atomlayers[0][atom.index]

# Gets composition of a lattice: up to four components
def get_comp(atoms,layers):

	countats = [ atom.number for atom in atoms if get_layer(atoms,atom) in layers ]
	unique = list(set(countats))
	compA = countats.count(unique[0])/len(countats)
	try:
		compB = countats.count(unique[1])/len(countats)
	except IndexError:
		compB = 0.
	try: 
		compC = countats.count(unique[2])/len(countats)
	except IndexError:
		compC = 0.
	try:
		compD = countats.count(unique[3])/len(countats)
	except IndexError:
		compD = 0.

	atomic_conc = [compA, compB, compC, compD]
	return compA
	
	
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
	def __init__(self, vertices, **kwargs):
		# cluster vertices in fractional coords
		self.vertices = vertices

		# desired cluster occupation if applicable
		# ordered vector of occupation variables 
		try:
			self.d_occ = ''.join(str(k) for k in kwargs['d_occ'])
		except KeyError:
			pass

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
		#prim_vectors, sites, numbers = refine_cell(spg_cell,symprec=1.e-3)
		self.prim_vectors = prim_vectors
		self.sites = sites
		# assumes 1 primary sublattice (change to 2 for rocksalt structure, etc.)
		self.sym= get_symmetry(spg_cell,symprec=1.e-3)['rotations']
		self.nsubs = 1
		# number of sublattices
		self.zsubs = zsubs
		self.sublattices = sublattices
		#self.sym = get_spacegroup(spg_cell,symprec=1e-5)

		return None
	# reduces symmetry operations to 2d
	def sym_2d(self):
		opps = [ opp for opp in self.sym if opp[2][2] == 1]
		#spacegroup = int(self.sym.split()[-1].strip('(').strip(')'))
		#opps = sym_2d_lib(spacegroup)
		#print (opps)
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
			var_map[num] = vrs[ind]
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
	# E.g. for now, we cannot do crystal structures that require translational operations to check for equivalent clusters

	def cluster_avg_2d(self):
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
				matchcount = 0
				positions = cluster.vertices.copy()

				# loop over vertices
				tmpsc = []
				for p in positions:
					pt = np.transpose(np.array([p]))
					transformed = np.matmul(opp,pt)
					tmpsc.append(np.transpose(transformed)[0])#.T)

				#Uncomment to print out symmetry operations in a cell
				
				'''#get atoms object for plotting
				ats = read(sys.argv[1])
				ats_sc = ats.get_scaled_positions().tolist()
				count +=1
				numbers = ats.get_atomic_numbers().tolist()
				for ptmp in np.matmul(np.add([1.,1.,0], tmpsc),sc_mat):
					numbers.append(17)
					print (ptmp)
					ats_sc.append(np.mod(ptmp,1.0))

				new_ats = Atoms(numbers)
				new_ats.set_atomic_numbers(numbers)
				new_ats.set_cell(ats.get_cell())

				new_ats.set_scaled_positions(ats_sc)
				#atoms.wrap()
				write('opp_%d.cif' % count, new_ats)'''
				
			#   ##convert cluster vertices to supercell basis
				scpositions = np.matmul( tmpsc , sc_mat)

				#raster cluster over the 2d cell
				for i in range(int(1/sc_mat[0][0])):
					for j in range(int(1/sc_mat[1][1])):
						tmp = np.add(np.multiply(i,sc_mat[0]) , scpositions)
						tmp = np.add(np.multiply(j,sc_mat[1]) , tmp)
						tmp = np.mod(tmp,1.)
						clst_pos.append(tmp)
				# goes through and gets occupation of vertices
				# This part is slow. Needs to be optimized
				for verts in clst_pos:
					vert_occ = []
					for vert in verts:
						for ind, site in enumerate(self.sites):
							dist = np.linalg.norm(np.subtract(site,vert))
							#NOTE this tolerance may lead to errors in huge supercells
							if dist <= 1.e-2:
								vert_occ.append(self.var_map[self.numbers[ind]])


					if len(vert_occ) == len(cluster.vertices): # "Mismatched cluster"
						# calculates cluster correlation function
						phi = phi_1(np.array(vert_occ), 2)
						clst_phi.append(phi)

						if cluster.d_occ != None:
							if ''.join(str(k) for k in vert_occ) == cluster.d_occ:
								matchcount +=1
						else:
							print ('Warning: Define desired occupation to calculate SRO parameter')


				all_match.append(matchcount/len(clst_phi))
				all_phi.append(np.average(clst_phi))

		return  (all_match , all_phi)

	def cluster_avg_3d(self,dof):
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
					tmpsc.append(np.transpose(transformed)[0])#.T)

				#Uncomment to print out symmetry operations in a cell
				
				'''#get atoms object for plotting
				ats = read(sys.argv[1])
				ats_sc = ats.get_scaled_positions().tolist()
				count +=1
				numbers = ats.get_atomic_numbers().tolist()
				for ptmp in np.matmul(np.add([1.,1.,0], tmpsc),sc_mat):
					numbers.append(17)
					print (ptmp)
					ats_sc.append(np.mod(ptmp,1.0))

				new_ats = Atoms(numbers)
				new_ats.set_atomic_numbers(numbers)
				new_ats.set_cell(ats.get_cell())

				new_ats.set_scaled_positions(ats_sc)
				#atoms.wrap()
				write('opp_%d.cif' % count, new_ats)'''
				
			#   ##convert cluster vertices to supercell basis
				scpositions = np.matmul( tmpsc , sc_mat)
			#   ##wrap back into cell
				##scpositions = np.mod(scpositions,1.0)

				#raster cluster over the 3d cell
				for i in range(int(1/sc_mat[0][0])):
					for j in range(int(1/sc_mat[1][1])):
						for k in range(int(1/sc_mat[2][2])):
							tmp = np.add(np.multiply(i,sc_mat[0]) , scpositions)
							tmp = np.add(np.multiply(j,sc_mat[1]) , tmp)
							tmp = np.add(np.multiply(k,sc_mat[2]) , tmp)
							tmp = np.mod(tmp,1.)
							clst_pos.append(tmp)
				# goes through and gets occupation of vertices
				# This part is slow. Needs to be optimized
				for verts in clst_pos:
					vert_occ = []
					for vert in verts:
						for ind, site in enumerate(self.sites):
							dist = np.linalg.norm(np.subtract(site,vert))
							#NOTE this tolerance may lead to errors in huge supercells
							if dist <= 1.e-2:
								vert_occ.append(self.var_map[self.numbers[ind]])

					# TODO implement wykoff positions to avoid this statement
					if len(vert_occ) == len(cluster.vertices): # "Mismatched cluster"0
				
						#get basis function
						phi = phi_1(np.array(vert_occ), 2)
						clst_phi.append(phi)

						if cluster.d_occ != None:
							if ''.join(str(k) for k in vert_occ) == cluster.d_occ:
								matchcount +=1
						else:
							print ('Warning: Define desired occupation to calculate SRO parameter')

				all_match.append(matchcount/len(clst_phi))
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
			
	return ( vectors,sites, numbers,sublattices)
