from utils.conformation_utils import *
from utils.design_utils import *

def test_chainbreak():
	is_chainbreak, secstruct_new = parse_out_chainbreak('(( ))')
	assert(is_chainbreak == [0,1,0,0])
	assert(secstruct_new == '(())')

def test_convert_structure_to_bps_1():
	bps = convert_structure_to_bps('((.))')
	assert(bps == [[1,3],[0,4]])

def test_convert_structure_to_bps_2():
	bps = convert_structure_to_bps('((.aa))aa')
	assert(bps == [[3, 8], [4, 7], [1, 5], [0, 6]])

def test_secstruct_to_partner():
	partner = secstruct_to_partner('((.))')
	assert(all(partner == [4., 3., -1., 1., 0.]))

def test_parse_stems_from_bps():

	stems = parse_stems_from_bps(convert_structure_to_bps('(.).(.)'))
	assert(stems == [[[0,2]], [[4,6]]])

def test_figure_out_stem_assignment():
	stem_assignment = figure_out_stem_assignment('(.).(.)')
	assert(all(stem_assignment == [1., 0., 1., 0., 2., 0., 2.]))

def test_get_nts():
	assert( get_nts_for_symbol( 'A' )==['A'] );
	assert(  get_nts_for_symbol( 'N' )==['A','C','G','U'] ); 
	assert( get_nts_for_symbol( 'W' )==['A','U'] );
	assert( get_sequences_for_pattern( 'AANUU' )==['AAAUU','AACUU','AAGUU','AAUUU'] );

if __name__=='__main__':
	test_chainbreak()
	test_convert_structure_to_bps_1()
	test_convert_structure_to_bps_2()
	test_parse_stems_from_bps()
	test_secstruct_to_partner()
	test_figure_out_stem_assignment()
	test_get_nts()
