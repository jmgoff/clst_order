# clst_order

## Summary

Post-processing tools for calculating short-range order (SRO) parameters in
substitutional lattice systems. These tools can be used to calculate the
Warren-Cowley pair paramters as well as the cluster order parameters defined in
the COP paper. The occurrence of arbitrary multi-point SRO motifs are systematically
quantified by direct calculation normalized cluster probabilities. The software
is written in python/cython for front-end accessiblity, and is compatible with
any crystal structure file compatible [with
ASE](https://wiki.fysik.dtu.dk/ase/ase/io/io.html). This includes trajectory
files from LAMMPS, RASPA, VASP, Quantum Espresso, etc. Trajectories from Monte
Carlo simulations in cluster expansion software packages such as CASM and ICET
can also be handled. If you are quantifying SRO in Monte Carlo simulations based
on cluster expansions, see the ICET branch for a fast on-the-fly cluster order
parameter calculator.

## Installation

Add the software directory to your pythonpath then run the clean.sh script

<pre><code>./clean.sh </code></pre>

### cython (recommended)


Link the appropriate gcc compiler in your bin. For example in an anaconda bin
directory (if using anaconda 3):

<pre><code>ln -s  x86_64-conda_cos6-linux-gnu-gcc gcc  </code></pre>

Add numpy libraries to the build. (can be found with np.__file__ after importing
numpy as np in a python script)

<pre><code>export
CFLAGS=-I/path/to/python/packages/site-packages/numpy/core/include</code></pre>

<pre><code>python setup.py build_ext --inplace</code></pre>

### Python only setup

<pre><code>./python_only.sh </code></pre>


## Usage

Python software is currently set up for single sublattice 2D/quasi-2D surfaces
as well as 3D single sublattice substitutional crystals within (m x n x l)
supercells. Functionality is being added for multiple sublattices in 2D and 3D
crystals, and for evaluation of correlation functions. A minimalistic example is
given below for a 5-layer AB surface alloy system with ASE Atoms object inputs.
It calculates the probability of an A-A-A-A cluster, which can then be divided
by marginal probabilities to obtain cluster order parameters. Currently, one
cluster probability with one occupation is measured per calculation. This will
soon be changed to the probabilities of all occupations in a single cluster.

<pre><code>from clst_prob import *
from lib_basis import *
prim_vectors, p_sites, p_numbers, p_sublattices = from_ase( prim_atoms , zsubs )

p = prim_cell(prim_vectors, p_sites, p_numbers, zsubs = [[0,1,2,3,4]],
sublattices=p_sublattices)

vectors, sites, numbers, sublattices = from_ase(supercell_atoms, zsubs)

clst_verts =[ [0.  , 0.  , 0.66290], 
[0.33333333 , 0.66666667 , 0.58145],
[0.66666667 , 0.33333333 , 0.50000], 
[0.00000 , 0.00000 , 0.41855] ]

#define all cluster vectors at these vertices with 2 degrees of freedom
clst =  cluster(clst_verts,get_variables(2))

supercell = cluster_cell(vectors,sites,numbers,p)
supercell.add_cluster(clst)

print( y.cluster_avg(is_slab=True))</code></pre>

Detailed descriptions of variables are provided in the clst_prob.py
script and example execution is provided in the examples folder. This
includes an example for parallelized execution of the software for several large
trajectory files.
