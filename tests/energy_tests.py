from utils.energy_utils import *
from conformation import Conformation

def test_score_motif():

	counts = score_motif(['((.))','.....'],'((.))')
	assert(all(counts == [1,0]))

	counts = score_motif(['((.))','.....'],'((x))')
	assert(all(counts == [1,0]))

	counts = score_motif(['((.))','.....'],'(( ))')
	assert(all(counts == [1,0]))

	counts = score_motif(['((.((.))','((......'],'(( ))')
	assert(all(counts == [1,0]))

def test_ener_0():
	conf = Conformation(secstruct='((.))', params = Parameters(sigma=0, epsilon=0, delta=0, 
	                                                           gnm=False, motifs=[Motif('((x))', dG=-1)]))
	conf.run()
	assert(conf.energies[0] == -1)

if __name__ == '__main__':
	test_score_motif()
	test_ener_0()