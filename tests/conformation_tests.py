from conformation import Conformation

def test_conf_0():

	secstruct_list = ['.....','((.))','(...)','.((...))',
	'(( ))', '.( )', '(. )', '( .)']
	#'..((.((..((...)))..)..(((...)))....))' ]

	n_conf_list = [2**4, 1, 3, 6, 1, 2, 2, 2,] # 4*126*3*3*3*3]

	for i in range(len(secstruct_list)):
		print(secstruct_list[i])
		mdl = Conformation(secstruct=secstruct_list[i])
		mdl.run()
		assert (mdl.n_conf == n_conf_list[i])

def test_conf_1():

	seq_list = ['NNN', 'AAA', 'UAUUA']
	n_conf_list = [5,4,28]

	for i in range(len(seq_list)):
		print(seq_list[i])
		mdl = Conformation(sequence=seq_list[i])
		mdl.run()
		assert (mdl.n_conf == n_conf_list[i])

def test_conf_2():

	secstruct_list = ['_(.)_','(___)']
	n_conf_list = [5,5]

	for i in range(len(secstruct_list)):
		print(secstruct_list[i])
		mdl = Conformation(secstruct=secstruct_list[i])
		mdl.run()
		assert (mdl.n_conf == n_conf_list[i])

if __name__ == '__main__':
	test_conf_0()
	print('test 0 done')
	test_conf_1()
	print('test 1 done')
	test_conf_2()
	print('test 2 done')
