from design import Design

def test_switch_0():
	#These two models should design the same way!

	mdl0 = Design(secstruct=['...(.)', '......'], n_switch_states=2, score_function = (lambda z1, z2: z1/z2),
             params=[Parameters(gnm=False)]*2)
	mdl0.run()
	assert(len(mdl0.sequences) == 16)

	mdl1 = Design(secstruct='...(.)', params=Parameters(gnm=False))
	mdl1.run()

	assert(len(mdl1.sequences) == 16)

	assert(mdl0.score[0] == mdl1.score[0])

if __name__ == '__main__':
	test_switch_0()

