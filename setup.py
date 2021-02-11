from setuptools import Extension,setup
from Cython.Build import cythonize
import numpy
import glob

from subprocess import call

call('export CFLAGS=-I/home/jgoff/anaconda3/lib/python3.7/site-packages/numpy/core/include/',shell = True)
call('cp ./cyth_funcs/* .', shell = True)

include_path = [numpy.get_include()]

for f in glob.glob('*.pyx'):
	prefix = f.split('.')[0]
	extensions =(
	Extension(prefix, [f],
		#typical build wit NPY_1_6 works. have not tested 1_7 API version
		#define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
		include_dirs=include_path,
	),
	)


	setup(
		ext_modules=cythonize(extensions,language_level=3),
	)
