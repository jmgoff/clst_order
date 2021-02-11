import numpy as np
cimport numpy as np
import cython

# python functions used to define cluster basis using
# different site basis functions obeying scalar prod.
# defined in Sanchez et al. 1984

DTYPE=np.float64

ctypedef np.float_t DTYPE_t

@cython.cdivision(True)

def get_variables(int dof):
	cdef int i,remove
	cdef np.ndarray variables
	if dof %2 ==0:
		i = dof/2
		variables = np.arange(-i,i+1,dtype=int)
		remove = list(variables).index(0)
		del variables[remove]
	if dof % 2 != 0:
		i = (dof -1)/2
		variables = np.arange(-i,i+1,dtype=int)
	return variables

# basis functions below are the trigonometric basis fcns
#  used in Van de Walle, ICET, & others: these are implemented for
#  arbitrary numbers of components, but are less smooth compared
#  to chebychev basis  (more derivative discontinuities)

#/-------------------------------------------------------/
def phi_trig(int n,np.ndarray sigma,int dof):
	cdef np.ndarray phi
	if n ==0:
		phi=np.ones(np.shape(sigma)[0])

	if n % 2 !=0:
		phi = -np.cos( ( np.pi * (n+1) * sigma )  /dof)

	if n % 2 ==0:
		phi = -np.sin( ( np.pi *  n  * sigma )  /dof)

	return np.prod(phi)
#/-------------------------------------------------------/



# basis functions below are the discrete chebychev basis
# improved smoothness over trig basis, but only implemented
# up to 4 degrees of freedom
#/-------------------------------------------------------/
def phi_0(np.ndarray sigma,int dof):
	return np.prod(np.ones(np.shape(sigma)))

def phi_1(np.ndarray sigma,int dof):
	if dof ==1:
		raise ValueError('Site degrees of freedom must be >= 2')
	if dof ==2:
		return np.prod(sigma)
	if dof ==3:
		return np.prod(np.multiply(np.sqrt(3/2),sigma ))
	if dof ==4:
		return np.prod(np.multiply(np.sqrt(2/5),sigma ))


def phi_2(np.ndarray sigma,int dof):
	if dof < 3:
		raise ValueError('Site degrees of freedom must be >= 3')
	if dof ==3:
		return np.prod(np.sqrt(3) - np.multiply(np.sqrt(3), np.multiply(sigma,sigma) ))
	if dof ==4:
		return np.prod( (np.sqrt(34)/3) - np.multiply(((5/3)*np.sqrt(2/17)), np.multiply(sigma,sigma) ) )


def phi_3(np.ndarray sigma,int dof):
	if dof <4:
		raise ValueError('Site degrees of freedom must be >= 4')
	if dof ==4:
		return np.prod(np.multiply((-17/3)*np.sqrt(1/10),sigma) + np.multiply((1/3)*np.sqrt(5/2),(np.multiply(sigma,np.multiply(sigma,sigma)))))
#/--------------------------------------------------------/
