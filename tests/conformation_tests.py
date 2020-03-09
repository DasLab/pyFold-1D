from conformation import Conformation

def test_conf_0():

	secstruct_test_list = ['.....','((.))','(...)','.((...))',
	'(( ))', '.( )', '(. )', '( .)',
	'..((.((..((...)))..)..(((...)))....))' ]

	n_conf_list = [2**4, 1, 3, 6, 1, 2, 2, 2, 4*126*3*3*3*3]

	for i in range(len(secstruct_test_list)):
		mdl = Conformation(secstruct=secstruct_test_list[i])
		mdl.run()
		assert (mdl.n_conf == n_conf_list[i])

def test_conf_1():
	seq_list = ['NNN','AAA','UAUUA']
	n_conf_list = [5,4,28]

	for i in range(len(seq_list)):
		mdl = Conformation(sequence=seq_list[i])
		mdl.run()
		assert (mdl.n_conf == n_conf_list[i])

if __name__ == '__main__':
	test_conf_0()
	test_conf_1()

