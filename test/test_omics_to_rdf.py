from unittest import TestCase
from scripts.omics_to_rdf import OmicsToRDF

class TestOmicsToRDF(TestCase):
	def test_should_convert_cell_id_to_bto_id(self):
		cell_id = 'DMS 53'
		bto_id = OmicsToRDF().get_bto_id(cell_id)
		self.assertTrue('BTO_' in bto_id)

	def test_should_convert_no_cell_id_to_statement_NO_ID(self):
		cell_id = 'DMS 53333333333'
		bto_id = OmicsToRDF().get_bto_id(cell_id)
		self.assertTrue('NO ID' in bto_id)

















