import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import pandas as pd

from time import sleep
from SPARQLWrapper import SPARQLWrapper
from enum import Enum
import argparse
from joblib import Parallel
import time
import warnings
from tqdm import tqdm

class DbName(Enum):
	COMMON = 1
	PRE = 2
	GCN = 3
	CCLE = 4
	GDSC = 5

EXT = '.txt'
CCLE_PREFIX = 'result_ccle_'
GDSC_PREFIX = 'result_gdsc_'


parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gcn_path', help='GCN入力用のパラメータファイルパス', default='../config/prm_rdf_to_gcn.tsv')
parser.add_argument('-gp', '--gdsc_path', help='GDSC/turtleのパラメータファイルパス',
                    default='../config/prm_gdsc_to_rdf.tsv')
parser.add_argument('-cp', '--ccle_path', help='CCLE/turtleのパラメータファイルパス',
                    default='../config/prm_ccle_to_rdf.tsv')
parser.add_argument('-pp', '--precessing_path', help='前処理パラメータファイルパス', default='../config/prm_preprocessing.tsv')
parser.add_argument('-cm', '--common_path', help='共通パラメータファイルパス', default='../config/prm_common.tsv')
parser.add_argument("-n", "--normalize", help="各オミックスデータ正規化の有無", action="store_true")
args = parser.parse_args()
warnings.filterwarnings('ignore')

def get_param(category, prm_name):
	try:
		if category == DbName.GDSC:
			sweep = pd.read_table(args.gdsc_path, index_col=0).loc[prm_name, 'value']
		elif category == DbName.CCLE:
			sweep = pd.read_table(args.ccle_path, index_col=0).loc[prm_name, 'value']
		elif category == DbName.GCN:
			sweep = pd.read_table(args.gcn_path, index_col=0).loc[prm_name, 'value']
		elif category == DbName.PRE:
			sweep = pd.read_table(args.precessing_path, index_col=0).loc[prm_name, 'value']
		else:
			sweep = pd.read_table(args.common_path, index_col=0).loc[prm_name, 'value']
	except:
		sweep = ''
		print('パラメータファイル、または指定のパラメータ名がありません。')
	return sweep.split(',')


def text_progessbar(seq, total=None):
    step = 1
    tick = time.time()
    while True:
     time_diff = time.time()-tick
     avg_speed = time_diff/step
     total_str = 'of %n' % total if total else ''
     print('step', step, '%.2f' % time_diff,
       'avg: %.2f iter/sec' % avg_speed, total_str)
     step += 1
     yield next(seq)

all_bar_funcs = {
    'tqdm': lambda args: lambda x: tqdm(x, **args),
    'txt': lambda args: lambda x: text_progessbar(x, **args),
    'False': lambda args: iter,
    'None': lambda args: iter,
}

def ParallelExecutor(use_bar='tqdm', **joblib_args):
    def aprun(bar=use_bar, **tq_args):
     def tmp(op_iter):
      if str(bar) in all_bar_funcs.keys():
       bar_func = all_bar_funcs[str(bar)](tq_args)
      else:
       raise ValueError("Value %s not supported as bar type"%bar)
      return Parallel(**joblib_args)(bar_func(op_iter))
     return tmp
    return aprun

class OmicsToRDF(object):
	def __init__(self):
		self.sparql = SPARQLWrapper(endpoint='http://sparql.hegroup.org/sparql/', returnFormat='tsv')

	def get_bto_id(self, cell_id):
		self.sparql.setQuery("""                                                                                                                                                                      
				    PREFIX obo: <http://www.geneontology.org/formats/oboInOwl#>
				    SELECT DISTINCT ?x
				    from <http://purl.obolibrary.org/obo/merged/BTO>
				    WHERE
				    {
				    {?x rdfs:label  "%s cell"^^<http://www.w3.org/2001/XMLSchema#string>.}
				    UNION { ?x obo:hasRelatedSynonym  "%s cell"^^<http://www.w3.org/2001/XMLSchema#string>. }
				    }                                                                                                                    
				""" % (cell_id, cell_id))

		results = {}
		start_idx, end_idx = -1, -1
		for retry in range(5):
			try:
				results = self.sparql.query().convert().decode()
				start_idx, end_idx = results.find('BTO_'), results.rfind('"')
				break

			except:
				import traceback
				traceback.print_exc()
				sleep(1)

		return results[start_idx:end_idx] if start_idx is not -1 else 'NO ID'

	def clean_cell_line_name(self, name):
		name = name.split('@')[0].replace(' ', '')
		return name
