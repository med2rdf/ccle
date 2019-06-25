import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from rdflib import Literal, Namespace, Graph, BNode
from rdflib.namespace import RDF, RDFS, DCTERMS, XSD

from sklearn.preprocessing import MinMaxScaler
from joblib import delayed

import re
import pickle
import uuid

from scripts.omics_to_rdf import *

EXT = ".ttl"

class CcleToRDF(OmicsToRDF):
	def __init__(self, snp_row, exp_column, sample_row):
		super().__init__()
		self.__result_path = get_param(DbName.CCLE, 'out_folder')[0]

		self.__id_map = None
		self.__sample = sample_row
		self.__exp = exp_column
		self.__snp = snp_row
		self.__is_normalize = False

		self.__ccle_ns = Namespace("http://med2rdf.org/ontology/ccle#")
		self.__assay_ns = Namespace("http://med2rdf.org/ccle/assay/")
		self.__cellline_ns = Namespace("https://portals.broadinstitute.org/ccle/page?cell_line=")
		self.__gene_ns = Namespace("https://portals.broadinstitute.org/ccle/page?gene=")
		self.__ensembl_ns = Namespace("https://www.ensembl.org/id/")
		self.__variants_ns = Namespace("http://med2rdf.org/ccle/variant/")
		self.__m2r_ns = Namespace("http://med2rdf.org/ontology/med2rdf#")
		self.__cellosaurus_ns = Namespace("https://web.expasy.org/cgi-bin/cellosaurus/")
		self.__hco_ns = Namespace("http://identifiers.org/hco/")
		self.__faldo_ns = Namespace("http://biohackathon.org/resource/faldo#")
		self.__sio_ns = Namespace("http://semanticscience.org/resource/")
		self.__uo_ns = Namespace("http://purl.obolibrary.org/obo/")
		self.__skos_ns = Namespace("http://www.w3.org/2004/02/skos/core#")

		self.g = None
		self.init_graph()

		if os.path.isfile(self.__result_path):
			os.remove(self.__result_path)


	@property
	def id_map(self):
		if self.__id_map is None:
			with open(get_param(DbName.COMMON, 'in_cellosaurus_ccle_file')[0], 'rb') as f:
				self.__id_map = pickle.load(f)
		return self.__id_map

	@id_map.setter
	def id_map(self, value):
		self.__id_map = value


	@property
	def snp(self):
		if self.__snp is None:
			with open(get_param(DbName.CCLE, 'in_ccle_snp_file')[0], 'rb') as f:
				self.__snp = pickle.load(f)
		return self.__snp

	@snp.setter
	def snp(self, value):
		self.__snp = value

	@property
	def exp(self):
		if self.__exp is None:
			with open(get_param(DbName.CCLE, 'in_ccle_exp_file')[0], 'rb') as f:
				self.__exp = pickle.load(f)

			self.__exp = self.__exp.rename(columns=lambda s: s.split(' ')[0])
			if self.__is_normalize:
				self.__exp[self.__exp.columns[2:]] =  MinMaxScaler().fit_transform(self.__exp[self.__exp.columns[2:]].values)
		return self.__exp

	@exp.setter
	def exp(self, value):
		self.__exp = value

	@property
	def sample(self):
		if self.__sample is None:
			self.__sample = pd.read_table(get_param(DbName.CCLE, 'in_ccle_exp_file')[0], low_memory=False)

		return self.__sample

	@sample.setter
	def sample(self, value):
		self.__sample = value


	def init_graph(self):
		if self.g is not None:
			self.g.remove((None, None, None))
			self.g.close()

		self.g = Graph()
		self.g.bind("ccle", self.__ccle_ns)
		self.g.bind("cellosaurus", self.__cellosaurus_ns)
		self.g.bind("cclec", self.__cellline_ns)
		self.g.bind("ccleg", self.__gene_ns)
		self.g.bind("cclev", self.__variants_ns)
		self.g.bind("hco", self.__hco_ns)
		self.g.bind("m2r", self.__m2r_ns)
		self.g.bind("faldo", self.__faldo_ns)
		self.g.bind("dct", DCTERMS)
		self.g.bind("sio", self.__sio_ns)
		self.g.bind("obo", self.__uo_ns)
		self.g.bind("skos", self.__skos_ns)

	def save_turtle(self, name):
		name = re.sub(r'\W', "_", name)
		out_path = self.__result_path + CCLE_PREFIX + name + EXT
		idx = 2
		while (os.path.isfile(out_path)):
			out_path = self.__result_path + CCLE_PREFIX + name + str(idx) + EXT
			print('change name:{}'.format(out_path))
			idx += 1

		self.g.serialize(destination=out_path, format="turtle")
		self.init_graph()


	def create_ccle_single_cell_turtle(self, is_debug=False):
		if is_debug:
			target_gene = ["ENSG00000141510.11", "ENSG00000143631.10", 'ENSG00000101384.7', 'ENSG00000270408.1']#"ENSG00000142615.7", "ENSG00000020633.14",
			self.exp = self.exp[self.exp['Name'].isin(target_gene)]

		self.exp = self.exp.assign(Name=self.exp.loc[:, 'Name'].apply(lambda x: x.split('.')[0]))
		cell_name = self.sample.loc[:, 'CCLE name'].values[0]
		id_name = self.get_cell_line_id(cell_name)
		gender = self.sample[self.sample["CCLE name"] == cell_name]["Gender"].values[0]
		primary = self.sample[self.sample["CCLE name"] == cell_name]["Site Primary"].values[0]

		if id_name == 'NO_ID':
			return self.g

		self.create_cell_line_turtle(cell_name, id_name, gender, primary)

		for idx in range(self.exp.shape[0]):
			self.create_gene_turtle(idx, id_name, cell_name)

		if not is_debug:
			self.save_turtle(id_name)

		return self.g

	def create_cell_line_turtle(self, cell_name, cell_id, gender, site_primary):
		bot_id = self.get_bto_id(cell_id)
		if not bot_id == 'NO ID':
			self.g.add((self.__cellline_ns[cell_name], RDF["type"], self.__uo_ns[bot_id]))

		self.g.add((self.__cellline_ns[cell_name], RDF["type"], self.__m2r_ns["CellLine"]))
		self.g.add((self.__cellline_ns[cell_name], RDFS["label"], Literal(cell_name)))
		self.g.add((self.__cellline_ns[cell_name], self.__ccle_ns["gender"], Literal(gender)))
		self.g.add((self.__cellline_ns[cell_name], self.__ccle_ns["site_primary"], Literal(site_primary)))

		self.create_snp_turtle(cell_name)

	def create_gene_turtle(self, gene_idx, cell_id, cell_name):
		ensembl = self.exp.loc[:, 'Name'].values[gene_idx]
		hgnc = self.exp.loc[:, 'Description'].values[gene_idx]


		self.g.add((self.__gene_ns[ensembl], RDF["type"], self.__m2r_ns["Gene"]))
		self.g.add((self.__gene_ns[ensembl], RDFS["label"], Literal(hgnc)))
		self.g.add((self.__gene_ns[ensembl], RDFS["seeAlso"], self.__ensembl_ns[ensembl]))

		self.create_assay_turtle(gene_idx, cell_id, cell_name, ensembl)

	def create_assay_turtle(self, gene_idx, cell_id, cell_name, ensembl):
		assay = ensembl + "@" + cell_name
		exp_val = self.exp.loc[:, cell_name].values[gene_idx]

		self.g.add((self.__assay_ns[assay], self.__ccle_ns["celline"], self.__cellline_ns[cell_name]))
		self.g.add((self.__assay_ns[assay], self.__ccle_ns["gene"], self.__gene_ns[ensembl]))

		val_node = str(uuid.uuid4())

		self.g.add((self.__assay_ns[assay], self.__sio_ns["SIO_000216"], BNode(val_node)))
		self.g.add((BNode(val_node), self.__sio_ns['SIO_000300'], Literal(str(exp_val), datatype=XSD.decimal)))
		self.g.add((BNode(val_node), RDF["type"], self.__uo_ns["STATO_0000206"]))
		self.g.add((self.__uo_ns["STATO_0000206"], RDFS['label'], Literal('RPKM')))

	def create_snp_turtle(self, cell_name):
		snp_df = self.snp[self.snp['Tumor_Sample_Barcode'] == cell_name]
		ref_alle = snp_df['Reference_Allele'].values
		tumor_alle = snp_df['Tumor_Seq_Allele1'].values
		anno_trans = snp_df['Annotation_Transcript'].values
		ncbi_build = snp_df['NCBI_Build'].values
		start_pos = snp_df['Start_position'].values
		chromo = snp_df['Chromosome'].values
		end_pos = snp_df['End_position'].values
		var_class= snp_df['Variant_Classification'].values
		var_type = snp_df['Variant_Type'].values

		for snp_idx in range(snp_df.shape[0]):
			var_node = str(uuid.uuid4())

			self.g.add((self.__cellline_ns[cell_name], self.__ccle_ns["has_variant"], BNode(var_node)))
			self.g.add((BNode(var_node), self.__m2r_ns["gene"], self.__m2r_ns["Variation"]))
			self.g.add((BNode(var_node), RDF["type"], self.__m2r_ns["Variation"]))
			self.g.add((BNode(var_node), self.__m2r_ns["reference_allele"], Literal(ref_alle[snp_idx])))
			self.g.add((BNode(var_node), self.__m2r_ns["alternative_allele"], Literal(tumor_alle[snp_idx])))
			self.g.add((BNode(var_node), self.__ccle_ns["annotation_transcript"], Literal(anno_trans[snp_idx])))
			self.g.add((BNode(var_node), self.__ccle_ns["hcobuild"], self.__hco_ns["GRCh{}".format(ncbi_build[snp_idx])]))
			self.g.add((BNode(var_node), self.__ccle_ns["variant_classification"], Literal(var_class[snp_idx])))
			self.g.add((BNode(var_node), self.__ccle_ns["variant_type"], Literal(var_type[snp_idx])))

			self.g.add((BNode(var_node), self.__faldo_ns["location"], BNode(var_node + '_p')))
			self.g.add((BNode(var_node + '_p'), RDF["type"], self.__faldo_ns["Region"]))
			self.g.add((BNode(var_node + '_p'), self.__faldo_ns["begin"], BNode(var_node + '_p' + '_b')))
			self.g.add((BNode(var_node + '_p'), self.__faldo_ns["end"], BNode(var_node + '_p' + '_e')))
			self.g.add((BNode(var_node + '_p' + "_b"), RDF["type"], self.__faldo_ns["ExactPosition"]))
			self.g.add((BNode(var_node + '_p' + "_b"), self.__faldo_ns["position"], Literal(start_pos[snp_idx])))
			self.g.add((BNode(var_node + '_p' + "_b"), self.__faldo_ns["reference"], self.__hco_ns[chromo[snp_idx]]))
			self.g.add((BNode(var_node + '_p' + "_e"), RDF["type"], self.__faldo_ns["ExactPosition"]))
			self.g.add((BNode(var_node + '_p' + "_e"), self.__faldo_ns["position"], Literal(end_pos[snp_idx])))
			self.g.add((BNode(var_node + '_p' + "_e"), self.__faldo_ns["reference"], self.__hco_ns[chromo[snp_idx]]))

	def get_cell_line_id(self, name):
		rtrn = self.id_map[self.id_map["CCLE"] == name]["ID"]
		if len(rtrn) == 0:
			return 'NO_ID'
		return rtrn.values[0].replace(' ', '')


def process(snp_row, exp_column, sample_row):
	ctr = CcleToRDF(snp_row, exp_column, sample_row)
	ctr.create_ccle_single_cell_turtle()
	del ctr


if __name__ == '__main__':
	with open(get_param(DbName.CCLE, 'in_ccle_snp_file')[0], 'rb') as f:
		snp = pickle.load(f)

	with open(get_param(DbName.CCLE, 'in_ccle_exp_file')[0], 'rb') as f:
		exp = pickle.load(f)

	with open(get_param(DbName.CCLE, 'in_ccle_sample_file')[0], 'rb') as f:
		sample = pickle.load(f)

	jobs = int(get_param(DbName.COMMON, 'n_jobs')[0])
	# SNPファイル、EXPファイル、SAMPLEファイルで全てに存在する細胞株名を取得
	cell_line_list = list(set(snp["Tumor_Sample_Barcode"]) & set(exp.columns[2:]) & set(sample["CCLE name"]))
	print("CCLEデータからRDF/turtleファイルに変換")
	aprun = ParallelExecutor(n_jobs=jobs)
	aprun(total=len(cell_line_list))(delayed(process)(snp[snp['Tumor_Sample_Barcode'] == name], exp.loc[:, ['Name', 'Description', name]], sample[sample['CCLE name'] == name]) for name in cell_line_list)

