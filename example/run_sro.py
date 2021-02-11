import sys
sys.path.append('/media/jgoff/03560cbf-785d-4eb7-81a9-2e4e509058fd/software/cluster_order_master')
from clst_prob import *
from multiprocessing import pool, cpu_count
from mchammer import *
from subprocess import call
import json
import os
import gc
import glob
from ase.io import read,write

atoms = read('atoms.cif')

#list of lists: each list containing atomic layer 
#  indices in each sublattice
zsubs= [list(range(5))]

#exctracts lattice sites,vectors and atom types from ASE atoms objects
prim_vectors, sites, numbers, sublattices = from_ase( atoms , zsubs )

# define primitive cell object
x = prim_cell(prim_vectors,sites,numbers, zsubs = zsubs, sublattices=sublattices)

#2D surface symmetry operations
sym2d = x.sym_2d()

#Vertices of the cluster in the primitive cell basis
clst_verts =[
[0.  , 0.  , 0.66290],
[0.33333333 , 0.66666667 , 0.58145],
[0.66666667 , 0.33333333 , 0.50000],
[0.00000 , 0.00000 , 0.41855]
]

#occupation of cluster in spin variables
# spin variables are assigned lowest spin variable 
# to lowest atomic number
occs= [-1,-1,-1,-1]

#cluster probabilities for a single microstate
import time
def single_prob(args):
	itr = args['itr']
	atoms = args['atoms']
	vectors, sites2, numbers2, sublattices2 = from_ase(atoms, zsubs)
	#define the cell the cluster will be rastered over "cluster_cell" object
	y = cluster_cell(vectors,sites2, numbers2, x)
	# define and append the cluster to the cluster_cell object
	clst =  cluster(clst_verts,d_occ=''.join(str(i) for i in occs))
	y.add_cluster(clst)
	# get cluster probability
	probabilities, clst_functions = y.cluster_avg_2d()
	return {'%05d' % itr:probabilities}

# evaluate cluster probabilities over 100 parallel MC chains
chain_length = 44000
for k in range(100):
	try:
		#this is the icet trajectory file
		#any other trajectory file can be replaced here
		dc = DataContainer.read( 'myrun_c_0.080_%d.dc' % k )
		name = 'example_%d'% k 

		if not os.path.isfile('%s.json' % name):
			icet_traj = dc.get('trajectory')

			args = []

			for m in range(chain_length)[1:]:
				args.append( { 'atoms':icet_traj[m],'itr': m } )
			ncpus = int(cpu_count())
			pl = pool.Pool(processes=ncpus)
			results = pl.map(single_prob, args)

			#Probabilities are stored in a JSON file
			res_dict = {}
			for i in range(chain_length)[1:]:
				res_dict[i] = None
			for result in results:
				for key,value in result.items():
					res_dict[int(key)] = value
			with open('%s.json' % name, 'w') as writejson:
				json.dump(res_dict, writejson, sort_keys=True, indent=2)

			#Clean up memory
			del(icet_traj)
			del(dc)
			del(results)
			del(res_dict)
			del(args)
			del(pl)
			gc.collect()

		else:
			pass
	except (FileNotFoundError,IndexError):
		print ('WARNING: Missing or corrupted trajectory', name,k)


