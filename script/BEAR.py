#!/usr/bin/env python
import sys, os
import re, itertools, math
import traceback
from argparse import Namespace
import yaml
import pandas as pd, numpy as np, numexpr as ne
from ctypes import c_int, c_double, CDLL, POINTER # for Bmonreg
import warnings
warnings.filterwarnings("ignore")

# ----------------------------- basic functions ----------------------------- #
def mkdir(path):
	path = re.sub('/$', '', path)
	if os.access(os.path.dirname(path), os.W_OK):
		os.makedirs(path, exist_ok=True)
	else:
		print(f'Permission Denied: {path}')

def search_sorted_subsets_field(df, field, keywords):
	subsets = []
	for keyword in keywords:
		left_bound = df[field].searchsorted(keyword,'left')[0]
		right_bound = df[field].searchsorted(keyword,'right')[0]
		if right_bound-left_bound > 0:
			subsets.append(df[left_bound:right_bound])
	return subsets

def df_query(df, searches, return_col = None):
	if not isinstance(df, pd.DataFrame):
		print('`df` is not a pandas DataFrame !')
		return
	if isinstance(searches, list): 
		# query for multiple columns: input as list
		subsets = [df]
		for field, values in searches:
			subsets = list(itertools.chain(*[search_sorted_subsets_field(subset, field, values) for subset in subsets]))
		result = pd.DataFrame(np.vstack(subsets) if len(subsets)>0 else None,columns=df.columns).convert_objects()
	elif isinstance(searches, tuple): 
		# query for sigle column using numexpr: input as tuple
		col = searches[0]; val = searches[1]
		comm = f'{col} == {val}'
		exec(f'{col} = df.{col}.values')
		result = df[ne.evaluate(comm)]
	else:
		print('`searches`: list of conditions or tuple allowed !')
		return
	if return_col:
		if isinstance(return_col, str):
			return result[return_col].values
		else:
			return result.loc[:,return_col]
	else:
		return result

def set_args(mydict, arg = None):
	if arg is None:
		return Namespace(**mydict)
	else:
		for k, v in mydict.items():
			setattr(arg, k, v)
		return arg

# ------------------------------ preprocess -------------------------------- #

def init(config_file_path = None):
	'''
		set config and load input data
	'''
	args = set_config(config_file_path)
	# set output directory	
	if args.out_dir is None:
		mkdir('./result_set')
		args.out_dir = f'./result_set/{args.title}'
	mkdir(args.out_dir)
	args.input = pd.read_pickle(args.input_file_path)
	args = load_train_flag(args)
	return args	

def set_config(config_file_path):
	'''
		set config from config file
	'''
	if os.path.isfile(config_file_path):
		with open(config_file_path) as f:
			args = set_args(yaml.load(f, Loader=yaml.FullLoader))			
	else:
		print('##  ERROR  ## `set_config`: config file cannot be loaded !')
		return
	return args

def load_train_flag(args):
	'''
		load pd.DataFrame(*.pkl or *.dat) for train_flag (train_flag)  
	'''
	args.train_flag = pd.read_pickle(args.train_flag_file_path)
	args.train_eaids = args.train_flag.eaid.unique().tolist()
	return args

# ------------------------------ BEAR functions --------------------------------- #

# -- 1) assay score

def HEA(q, a, h):
	'''
		q: query compounds
		a: total compounds in the given assay
		h: hit compounds in the given assay
	'''	
	ecid2reg_score = None
	nQ = len(q)
	score, nA, nH, nHQ = get_hea_score(q, a, h)
	bin_n = math.floor(nA/args.min_comp_per_bin)
	max_bin_n = max(nQ, math.floor(nA/min(nH, args.max_comp_per_bin)))
	if bin_n >= args.min_bin_n:	
		comp_per_bin = args.min_comp_per_bin if bin_n <= max_bin_n else min(math.floor(nA/max_bin_n), args.max_comp_per_bin)
		tmp = []; bn = 1; low_limit = 0
		while bn <= min(bin_n, max_bin_n) and low_limit <= args.max_n_low_limit:	
			start_i = (bn-1)*comp_per_bin
			end_i = bn*comp_per_bin-1
			if bn == min(bin_n, max_bin_n):
				end_i = nA-1
			B = a[start_i:end_i]
			B_LR, _, _, _ = get_hea_score( q, a, B )
			if start_i > nH-1 and round(B_LR) <= 1 :
				low_limit = low_limit+1
			tmp.append(pd.DataFrame({'ecid':B, 'score':B_LR}))
			bn = bn + 1
		bin_LR_df = pd.concat(tmp)
		fit_LR = Bmonreg(np.array(range(len(bin_LR_df))), bin_LR_df.score)
		if fit_LR is not None:
			fit_LR = list(filter(lambda x: x and x > args.logscale**args.min_llr_score, fit_LR))
			if len(fit_LR) != 0:
				fit_LR = list(map(math.log2, fit_LR))
				ecid2reg_score = dict(zip(a[0:len(fit_LR)-1], fit_LR))
	return {'e_hea': score, 'lenQH_hea': nHQ, 'ecid2reg_score': ecid2reg_score}

def get_hea_score(q, a, h):
	'''
		q: query compounds (list)
		a: all compounds in the given assay (list)
		h: hit compounds in the given assay (list)
	'''
	score = args.min_hea_score
	# HEA parameters
	nA = max(args.min_a_size, len(a))
	nH = len(h)
	nQ = len(set(q).intersection(a))
	nHQ = len(set(q).intersection(h))
	nHQc = nH - nHQ; nQc = nA - nQ
	if nH != 0:
		HQ_ratio = (nHQ + args.laplacian_h) / (nHQc + args.laplacian_h) # escape 0 in (nH-nHQ).
		AQ_ratio = (nQ + 1) / (nQc + 1) # escape 0 in (nA-nQ).
		score = HQ_ratio/AQ_ratio
	return score, nA, nH, nHQ

def Bmonreg ( x, y, hd=0.3, hr=0.3, Kd="epanech", Kr="epanech", degree=0, inverse=0, monotonie="antinton" ):
	'''
		x, y: numpy array
	'''
	# Load the shared object file
	mon = CDLL('./script/monreg.so')

	# Casting Python numbers to ctype pointer
	c_int_p = POINTER(c_int)
	c_double_p = POINTER(c_double)

	def c_cast_int(pyobj):
		if isinstance(pyobj, int):
			cobj = np.array([pyobj], dtype = np.int).ctypes.data_as(c_int_p)
		elif isinstance(pyobj, list):
			cobj = np.array(pyobj, dtype = np.int).ctypes.data_as(c_int_p)
		elif isinstance(pyobj, np.ndarray):
			cobj = pyobj.ctypes.data_as(c_int_p)
		elif isinstance(pyobj, pd.Series):
			cobj = pyobj.to_numpy(dtype = np.int).ctypes.data_as(c_int_p)
		else:
			print('No valid python object is entered for `c_int_p` ')
			return
		return cobj

	def c_cast_double(pyobj):
		if isinstance(pyobj, float) or isinstance(pyobj, int):
			cobj = np.array([pyobj], dtype = np.double).ctypes.data_as(c_double_p)
		elif isinstance(pyobj, list):
			cobj = np.array(pyobj, dtype = np.double).ctypes.data_as(c_double_p)
		elif isinstance(pyobj, np.ndarray):
			cobj = pyobj.ctypes.data_as(c_double_p)
		elif isinstance(pyobj, pd.Series):
			cobj = pyobj.to_numpy(dtype = np.double).ctypes.data_as(c_double_p)
		else:
			print('No valid python object is entered for `c_double_p` ')
			return
		return cobj	

	# Declare input & output type for the function
	mon.mdach_a_inv.argtypes = [c_double_p, c_double_p, c_double_p, c_double_p, c_int_p, c_int_p, c_int_p, c_int_p, c_double_p, c_double_p, c_int_p, c_int_p, c_double_p, c_int_p, c_int_p, c_double_p]
	mon.mdach_a_inv.restype = None

	mon.mdach_i_inv.argtypes = [c_double_p, c_double_p, c_double_p, c_double_p, c_int_p, c_int_p, c_int_p, c_int_p, c_double_p, c_double_p, c_int_p, c_int_p, c_double_p, c_int_p, c_int_p, c_double_p]
	mon.mdach_i_inv.restype = None

	# Regression start
	a = min(x)
	b = max(x)
	N = len(x.shape)
	t = len(x)
	ldeg = len(degree) if isinstance(degree, list) else 1
	if ldeg > 1 and ldeg != N:
		print("lens of N and degree differ")
		return None

	xs = (x-a)/(b-a)

	if N == 1 and ldeg == 1:
		z = np.linspace(min(xs),max(xs),num=t)
	else:
		z = x.shape[0]

	if N == 1:
		lt = t
		t = np.repeat(0, t)
		tflag = 1
	else:
		lt = len(t)
		tflag = 0

	Kern1 = {'epanech':1, 'rectangle':2, 'biweight':3, 'triweight':4, 'triangle':5, 'cosine':6}.get(str(Kd), 1)
	Kern2 = {'epanech':1, 'rectangle':2, 'biweight':3, 'triweight':4, 'triangle':5, 'cosine':6}.get(str(Kr), 1)

	lx = len(x)
	lz = len(z)
	erg = np.linspace(1, 1, num=lt)

	result = None
	try:
		if len(x) != len(y):
			print("x and y must have the same len")
		else:
			if inverse != 0 and inverse != 1:
				print ("inverse must have value 0 or 1")
			else:
				if monotonie == "antiton":
					mon.mdach_a_inv(
						c_cast_double(xs), c_cast_double(y), c_cast_double(z), c_cast_double(t),
						c_cast_int(lx), c_cast_int(lz), c_cast_int(lt), c_cast_int(tflag), 
						c_cast_double(hd), c_cast_double(hr),
						c_cast_int(Kern1), c_cast_int(Kern2), c_cast_double(degree), 
						c_cast_int(ldeg), c_cast_int(inverse),
						c_cast_double(erg)
					)
				else:
					mon.mdach_i_inv(
						c_cast_double(xs), c_cast_double(y), c_cast_double(z), c_cast_double(t),
						c_cast_int(lx), c_cast_int(lz), c_cast_int(lt), c_cast_int(tflag), 
						c_cast_double(hd), c_cast_double(hr),
						c_cast_int(Kern1), c_cast_int(Kern2), c_cast_double(degree), 
						c_cast_int(ldeg), c_cast_int(inverse),
						c_cast_double(erg)
					)
	except Exception as e:
		pass
	# else:
	if inverse == 0:
		result = erg.copy()
	else:
		result = (erg*(b-a)+a).copy()
	return result


# -- 2) compound score

def BEAR(eaid2score):
	global args
	''''''
	hea_result = eaid2score.loc[(eaid2score.lenQH_hea >= args.min_oc_in_a) & (eaid2score.e_hea > args.logscale**args.min_hea_score), ['eaid','e_hea','ecid2reg_score','lenQH_hea']]
	if len(hea_result) > 0:
		ecid2total_score = pd.concat([get_bear_score(row) for _, row in hea_result.iterrows()])
		ecid2bear_score = (
			ecid2total_score
				.groupby('ecid', as_index = False, sort = False).bear_score.sum()
				.sort_values('bear_score', ascending = False)
				.reset_index(drop = True)
		)
		return ecid2bear_score
	else:
		return None

def get_bear_score(erow):
	global args
	'''
		erow: a certain row of args.eaid2score for the given q_eaid
		required columns: ['eaid','e_hea','ecid2reg_score']
	'''
	q_eaid, e_hea = erow.eaid, erow.e_hea
	H, NH = df_query(args.train_flag, ('eaid', q_eaid), ['h','nh']).iloc[0]
	es_hea = np.log(e_hea)/np.log(args.logscale)

	if erow.ecid2reg_score is None or pd.isnull(erow.ecid2reg_score):
		ecid2score = pd.concat([
			pd.DataFrame({'ecid': NH, 'bear_score': [0]*len(NH)}), 
			pd.DataFrame({'ecid': H, 'bear_score': [es_hea]*len(H)})
		])
	else:
		ecid2score = pd.DataFrame({'ecid': H + NH})
		ecid2score['bear_score'] = ecid2score.ecid.map(erow.ecid2reg_score)
		ecid2score.bear_score.fillna(0, inplace = True)
	return ecid2score

def get_bioassay_score(qrow):
	global args
	query_target = qrow.query_target
	try:
		eaid2score = get_eaid2score(qrow)
		if eaid2score is not None:
			mkdir(args.out_dir + '/eaid2score')	
			eaid2score.reset_index(drop = True, inplace = True)	
			eaid2score.to_pickle(f'{args.out_dir}/eaid2score/{str(query_target)}.pkl')

			print(f'## Finished ## `get_bioassay_score` for {str(query_target)} RESULT: {args.out_dir}/eaid2score/{str(query_target)}.pkl')	

			return eaid2score
	except Exception as e:
		traceback.print_tb(e.__traceback__)	

def get_eaid2score(qrow):
	global args
	'''
	<input>
		qrow: Target ligands (or hits) [type: dict or pandas.Series]
	<required vartiabales>
		train_flag : DataFrame of {'h': list of int64, 'nh': list of int64}}
			- h: Hit ecids from the bioassay
			- nh: Non-hit ecids from the bioassay
		param
	<return>
		DataFrame of HEA score by eaid
	'''
	eaid2score = None
	result = []
	for eaid in args.train_eaids:
		result.append(calc_scores(qrow, eaid))
	if len(result) > 0:
		eaid2score = pd.concat(result)
	return eaid2score

def calc_scores(qrow, eaid): 
	global args
	'''
	<input>
		qrow: {'query_target': string, 'list': list of int64, 'negative_list': list of int64}
			- query_target: title of query target
			- list: query compounds related to the query_target of interest
			- negative_list: negative compounds (not related to the query_target of interest)
	<return>
		Score : a single score of the bioassay for the query_target
	'''
	# -- calculate scores by assay
	assay_outcomes = args.train_flag[args.train_flag.eaid == eaid]
	e_res = assay_outcomes.apply(lambda row: pd.Series(HEA(qrow.list, row.h + row.nh, row.h)), axis=1)
	e_res = pd.concat([assay_outcomes['eaid'], e_res], axis=1) # column bind
	return e_res

def get_compound_score(row, eaid2score = None):	
	global args
	query_target = row.query_target
	mkdir(args.out_dir + '/ecid2score/')	

	# -- Get_score_by_compound
	Q = args.input[args.input.query_target == query_target].iloc[0].list
	N = args.input[args.input.query_target == query_target].iloc[0].negative_list

	ecid2bear_score = BEAR(eaid2score)
	if ecid2bear_score is None:
		print(f'## WARNING ## {query_target}: No result for BEAR !')
	else:
		ecid2score = ecid2bear_score[~ecid2bear_score.ecid.isin(Q)][~ecid2bear_score.ecid.isin(N)]
		ecid2score.bear_score.fillna(0, inplace = True)	
		ecid2score.to_pickle(f'{args.out_dir}/ecid2score/{query_target}.pkl')

if __name__ == "__main__":
	if len(sys.argv) > 1:
		config_f = sys.argv[1]
	else:
		currentdir = os.path.dirname(os.path.realpath(__file__))
		config_f = f'{currentdir}/config.yml'
	
	os.chdir(os.path.dirname(config_f) + '/../')
	
	# Initialization
	args = init(config_f)
	'''
		1) `get_bioassay_score` for each query_target input
			- `get_eaid2score`
				- `calc_scores` for each eaid
		2) `get_compound_score` for each query_target
	'''
	for inp_row in args.input.itertuples():
		eaid2score = get_bioassay_score(inp_row)
		if eaid2score is not None:
			get_compound_score(inp_row, eaid2score)
