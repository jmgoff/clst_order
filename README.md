# clst_order
Multi-point order parameters used to quantify chemical ordering in crystal systems.

This library written in python 3.7 can be used to quantify chemical ordering in substitutional
crystal systems on single sublattices. Future developments will extend this to crystal systems
with multiple non-exchanging sublattices, and systems with 3+ components. In the current
version, many-point ordering parameters can be calculated for binary alloy systems with both
3D and 2D symmetry, or the code can be used to calculate the Warren-Cowley pair parameters.
Due to it's interface with the Atomistic Simulation Environment package, it can be used to
evaluate chemical ordering from a variety of simulation types and software packages as a
post-processing tool.

Requirements:

Python 3.7.1 or later
spglib v. 1.16.0 (should be compatible with many versions though - this was the version used for the tests)
ASE v. 3.16.0 or later

Both required modules can be installed via pip:

pip install ase
pip install spglib

see examples for usage
