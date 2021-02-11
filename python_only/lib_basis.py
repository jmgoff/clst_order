import numpy as np

def get_variables(dof):
	if dof %2 ==0:
		i = dof/2
		i = int(i)
		variables = np.arange(-i,i+1).tolist()
		remove = variables.index(0)
		variables.pop(remove)
	if dof %2 != 0:
		i = (dof -1)/2
		i = int(i)
		variables = np.arange(-i,i+1).tolist()

	return variables

# basis functions below are the trigonometric basis fcns
#  used in Van de Walle, ICET, & others: these are implemented for
#  arbitrary numbers of components, but are less smooth compared
#  to chebychev basis  (more derivative discontinuities)

#/-------------------------------------------------------/
def phi_trig(n,sigma,dof):
    sigma = np.array(sigma)
    if n ==0:
        phi=np.ones(len(sigma))

    if n % 2 !=0:
        phi = -np.cos( ( np.pi * ( n +1) * sigma)  /dof)

    if n % 2 ==0:
        phi = -np.sin( ( np.pi *  n  * sigma)  /dof)

    return np.prod(phi)
#/-------------------------------------------------------/



# basis functions below are the discrete chebychev basis
# improved smoothness over trig basis, but only implemented
# up to 4 degrees of freedom
#/-------------------------------------------------------/
def phi_0(sigma,dof):
	return np.ones(np.shape(sigma))

def phi_1(sigma,dof):
	if dof ==1:
		raise ValueError('Site degrees of freedom must be >= 2')
	if dof ==2:
		return np.prod(sigma)
	if dof ==3:
		return np.prod(np.multiply(np.sqrt(3/2),sigma ))
	if dof ==4:
		return np.prod(np.multiply(np.sqrt(2/5),sigma ))


def phi_2(sigma,dof):
	if dof < 3:
		raise ValueError('Site degrees of freedom must be >= 3')
	if dof ==3:
		return np.prod(np.sqrt(3) - np.multiply(np.sqrt(3), np.multiply(sigma,sigma) ))
	if dof ==4:
		return np.prod( (np.sqrt(34)/3) - np.multiply(((5/3)*np.sqrt(2/17)), np.multiply(sigma,sigma) ) )


def phi_3(sigma,dof):
	if dof <4:
		raise ValueError('Site degrees of freedom must be >= 4')
	if dof ==4:
		return np.prod(np.multiply((-17/3)*np.sqrt(1/10),sigma) + np.multiply((1/3)*np.sqrt(5/2),(np.multiply(sigma,np.multiply(sigma,sigma)))))
#/--------------------------------------------------------/
