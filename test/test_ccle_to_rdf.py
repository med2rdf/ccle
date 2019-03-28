from unittest import TestCase
from scripts.ccle_to_rdf import CcleToRDF
from scripts.omics_to_rdf import *

from rdflib import Literal, Namespace
from rdflib.namespace import RDFS
import pickle

class TestCcleToRDF(TestCase):

	def setUp(self):
		# 入力データ
		with open(get_param(DbName.CCLE, 'in_ccle_snp_file')[0], 'rb') as f:
			self.snp = pickle.load(f)
		with open(get_param(DbName.CCLE, 'in_ccle_exp_file')[0], 'rb') as f:
			self.exp = pickle.load(f)
		with open(get_param(DbName.CCLE, 'in_ccle_sample_file')[0], 'rb') as f:
			self.sample = pickle.load(f)
		# 出力先
		self.output_path = get_param(DbName.CCLE, 'out_folder')[0] + 'debug_rdf.txt'
		# Namespace
		self.__ccle_ns = Namespace("https://portals.broadinstitute.org/ccle/")
		self.__skos_ns = Namespace("http://www.w3.org/2004/02/skos/core#")
		self.__cellline_ns = Namespace("https://portals.broadinstitute.org/links/cell_line/")
		# SPARQL
		self.sparql = SPARQLWrapper(endpoint='http://sparql.hegroup.org/sparql/', returnFormat='tsv')

		self.graph = None

	def tearDown(self):
		if self.graph is not None:
			self.graph.serialize(destination=self.output_path, format="turtle")
			with open(self.output_path) as f:
				result = f.read()
			print(result)


	def test_should_create_DMS53_cell_line_graph(self):
		snp_row = self.snp[self.snp['Tumor_Sample_Barcode'] == 'DMS53_LUNG']
		exp_column = self.exp.loc[:, ['Name', 'Description','DMS53_LUNG']]
		smpl_row = self.sample[self.sample['CCLE name'] == 'DMS53_LUNG']

		ctr = CcleToRDF(snp_row, exp_column, smpl_row)
		self.graph = ctr.create_ccle_single_cell_turtle(True)

		result = None
		for s1, p1, o1 in self.graph.triples((None, RDFS['label'], None)):
			if self.__cellline_ns['DMS53'] == s1:
				result = o1

		expected = Literal('DMS53_LUNG')
		self.assertEqual(result, expected)

	def test_should_create_22rv1_cell_line_graph(self):
		snp_row = self.snp[self.snp['Tumor_Sample_Barcode'] == '22RV1_PROSTATE']
		exp_column = self.exp.loc[:, ['Name', 'Description','22RV1_PROSTATE']]
		smpl_row = self.sample[self.sample['CCLE name'] == '22RV1_PROSTATE']

		ctr = CcleToRDF(snp_row, exp_column, smpl_row)
		self.graph = ctr.create_ccle_single_cell_turtle(True)

		result = None
		for s1, p1, o1 in self.graph.triples((None, RDFS['label'], None)):
			if self.__cellline_ns['22Rv1'] == s1:
				result = o1

		expected = Literal('22RV1_PROSTATE')
		self.assertEqual(result, expected)

	def test_should_create_rch_acv_cell_line_graph(self):
		snp_row = self.snp[self.snp['Tumor_Sample_Barcode'] == 'RCHACV_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE']
		exp_column = self.exp.loc[:, ['Name', 'Description','RCHACV_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE']]
		smpl_row = self.sample[self.sample['CCLE name'] == 'RCHACV_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE']

		ctr = CcleToRDF(snp_row, exp_column, smpl_row)
		self.graph = ctr.create_ccle_single_cell_turtle(True)

		result = None
		for s1, p1, o1 in self.graph.triples((None, RDFS['label'], None)):
			if self.__cellline_ns['RCH-ACV'] == s1:
				result = o1

		expected = Literal('RCHACV_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE')
		self.assertEqual(result, expected)

