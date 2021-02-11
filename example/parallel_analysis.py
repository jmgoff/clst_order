import sys
import numpy as np
import pandas as pd
import json
import glob
from multiprocessing import pool

dat_label = sys.argv[1]

#target concentration
comp_labels = [0.080]

#actual concentrations in finite supercell 
#  (double check in your trajectory)
lst = [0.08 ]

#---------------------------------------------------------------
#Some windowing/autocorrelation procedures adapted from MCHammer

# Automated windowing procedure from Sokal (1989)
def auto_window(taus, c):
	m = np.arange(len(taus)) < c * taus
	if np.any(m):
		return np.argmin(m)
	return len(taus) - 1

# From Goodman & Weare (2010)
def acl(autocorr,c=5.0):
	f = autocorr.copy()
	taus = 2.0 * np.cumsum(f) - 1.0
	window = auto_window(taus, c)
	return taus[window]

def acf(data, max_lag):
	if max_lag is None:
		max_lag = len(data) - 1
	if 1 > max_lag >= len(data):
		raise ValueError('max_lag should be between 1 and len(data)-1.')
	series = pd.Series(data)
	acf = [series.autocorr(lag) for lag in range(0, max_lag)]
	return np.array(acf)
#---------------------------------------------------------------

# 4000 samples here = 100 MC passes for burn in
burn_in = 4000

def single_comp(args):
	comp = args['comp']
	comp_label = comp_labels[lst.index(comp)]
	tmp = {}
	auto = {}
	try:
		for d in range(100):
			f = '%s_%1.3f_%d.json'% (dat_label,comp_label,d)
			tmp[d] = []
			auto[d] = []
			dat = []
			with open(f, 'r') as readin:
				data = json.load(readin)
			ks = [int(m) for m in data.keys()]
			for key in sorted(ks):
				key = str(key)
				if np.average(data[key]) <= comp: # Maximum prob is comp/degen
					dat.append(1- ( np.average(data[key]) / ((comp)**4) ))
				else:
					#TODO address occasional symmetry detection error here
					pass
			tmp[d] = dat
			autocorr = acf(dat,len(dat) - 1)
			auto[d] = np.nan_to_num(autocorr)

	except (FileNotFoundError, ValueError):
		print('Missing JSON file',f)
		pass
	sub_vrs = []
	sub_avg = []
	sub_acl = []

	for i in tmp.keys():
		if len(tmp[i]) != 0:
			sub_vrs.append(np.var(tmp[i][burn_in:]))
			sub_avg.append(np.average(tmp[i][burn_in:]))
			ac = acl(auto[i])
			ac = np.nan_to_num(ac)
			sub_acl.append(ac)
		else:
			print ('No probability data for parallel mc calculation:',i)
			pass

	#effective sample size 
	n_eff = len(tmp[0][burn_in:])/( (2*np.average(sub_acl)) +1)
	#number of subsamples
	N = len(sub_avg)
	#effective variance (assumed to come from the same distribution 
	#  for all subsamples : we average to get this distribution
	var = np.average(sub_vrs)
	std_err = np.sqrt( (var**2/(N*n_eff) ) + ( np.var(sub_avg) / N ) )
	return {'comp':comp, 'result': np.average(sub_avg),'error': std_err}


args = []
for comp in lst: #np.arange(0,0.2, 0.008)[1:]:
	args.append({'comp':comp})

pl = pool.Pool(processes=20)
results = pl.map(single_comp, args)

res_dict = {'comps': [], 'results':[], 'errors':[]}

for result in results:
	key = result['comp']
	res_dict['comps'].append(key)
	res_dict['results'].append(result['result'])
	res_dict['errors'].append(result['error'])

with open('%s_results.json' %dat_label, 'w') as writejson:
	json.dump(res_dict, writejson, sort_keys=True, indent=2)
