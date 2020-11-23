import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import norm
mpl.rcParams['font.family'] = 'Avenir Next'
mpl.rcParams['font.size'] = 18
mpl.rcParams['font.weight'] = 'medium'


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

# basis functions below are the discrete chebychev basis
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

# basis functions can be plotted - useful for comparing
#    smoothness to other site basis function types
def plot_basis():
	x = get_variables(4)

	y1 = phi_1(x,4)
	y2 = phi_2(x,4)
	y3 = phi_3(x,4)
	x = [int(i) for i in x]
	fig, ax = plt.subplots()
	plt.xlabel('Spin variable $\sigma$')
	plt.ylabel('$\phi_m$ ($\sigma$)')
	plt.ylim([-1.5,1.5])
	plt.xlim([-2,2])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.spines['bottom'].set_position(('outward',10))
	ax.spines['left'].set_position(('outward',5))
	plt.tight_layout()


	legend = []
	plt.plot(x, [1]*len(x), linestyle='dotted')
	legend.append(r'$\phi_0$')
	plt.plot(x,y1,linestyle='dashed')
	legend.append(r'$\phi_1$')
	plt.plot(x,y2,linestyle='dashdot')
	legend.append(r'$\phi_2$')
	plt.plot(x,y3,linestyle='solid')
	legend.append(r'$\phi_3$')
	plt.legend(legend,loc='best',prop={"size":10})
	
	fig.savefig('plot_basis.png')#, transparent=True)

#plot_basis()
