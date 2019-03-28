import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import pickle
from sklearn.preprocessing import MinMaxScaler
from scripts.omics_to_rdf import *

def make_id_map(text, db_name='CCLE'):
	id_txt, db_txt = None, None
	idx = 0
	df = pd.DataFrame([], columns=['ID', db_name])
	for txt in text:
		tmp = txt.strip('\n').split('   ')
		if tmp[0] == 'ID':
			id_txt = tmp[1]
		elif tmp[0] == 'DR':
			tmp[1] = tmp[1].split('; ')
			if tmp[1][0] == db_name:
				db_txt = tmp[1][1]

		if (tmp[0] == '//') & (db_txt is not None):
			idx += 1
			df.loc[str(idx)] = [id_txt, db_txt]
			id_txt, db_txt = None, None
	return df


def check_raw_file(is_norm):
	if not os.path.isfile(get_param(DbName.CCLE, 'in_ccle_snp_file')[0]):
		print("CCLE SNPファイルをpickle化")
		df = pd.read_csv(get_param(DbName.PRE, 'in_ccle_snp_file')[0], sep='\t',low_memory=False)
		with open(get_param(DbName.CCLE, 'in_ccle_snp_file')[0], 'wb') as f:
			pickle.dump(df, f)
		del df

	if not os.path.isfile(get_param(DbName.CCLE, 'in_ccle_exp_file')[0]):
		print("CCLE EXPファイルをpickle化")
		df = pd.read_csv(get_param(DbName.PRE, 'in_ccle_exp_file')[0], skiprows=2, sep='\t', low_memory=False)
		df = df.rename(columns=lambda s: s.split(' ')[0])
		if is_norm:
			print("CCLE EXPファイル正規化")
			df = df.T
			df.loc[2:, :] = MinMaxScaler().fit_transform(df.iloc[2:, :].values)
			df = df.T
		with open(get_param(DbName.CCLE, 'in_ccle_exp_file')[0], 'wb') as f:
			pickle.dump(df, f)
		del df

	if not os.path.isfile(get_param(DbName.CCLE, 'in_ccle_sample_file')[0]):
		print("CCLE SAMPLEファイルをpickle化")
		df = pd.read_csv(get_param(DbName.PRE, 'in_ccle_sample_file')[0], sep='\t', low_memory=False)
		with open(get_param(DbName.CCLE, 'in_ccle_sample_file')[0], 'wb') as f:
			pickle.dump(df, f)
		del df

	if not os.path.isfile(get_param(DbName.COMMON, 'in_cellosaurus_ccle_file')[0]):
		print("cellosaurusファイルをccleについてpickle化")
		with open(get_param(DbName.PRE, 'in_cellosaurus_file')[0]) as f:
			text = f.readlines()
		df = make_id_map(text, 'CCLE')
		with open(get_param(DbName.COMMON, 'in_cellosaurus_ccle_file')[0], 'wb') as f:
			pickle.dump(df, f)
		del df

check_raw_file(args.normalize)
