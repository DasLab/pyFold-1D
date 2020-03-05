from utils.design_utils import *
import matplotlib.pyplot as plt

def enumerative_design( secstruct, pattern=None, params=None)
#   Enumerate over all sequences that could form target secondary
#    structure (and conform to optional design pattern) and see which ones
#    fold best.
#
# Inputs
# secstruct = Target secondary structure in dot-parens notation, e.g.
#                 '((.))'. Give [] or '' if you want
#                  the secondary structures to be enumerated
# pattern = [Optional] sequence like AANUU, where characters like 
#                N,R,Y,W, and S are wild cards, and A,U,C, and G are preserved.
#                If not specified, code will use ANNNNNN... (note that
#                first base can be set to A without loss of generality in
#                current energy model)
#  params = Energy parameter values for delta, epsilon, etc. [MATLAB struct]
#
#
# Outputs
# sequences = all sequences that match secstruct and pattern, ordered
#              with 'best' design (by p_target) first
# p_target  = For each sequence, fraction of conformations with target 
#                secondary structure.
# x = For each sequence, all sets of conformations.
#        If there are no base pairs specified, should get
#        2^(Nbeads-1). First position is always 0.   
# d = For each sequence, input directions (array of +/-1's)
# p = For each sequence, partners  (0 if bead is unpaired,
#        otherwise index of partner from 1,... Nbeads )

	if params is None:
		params = Parameters()

	if pattern is None:
		pattern = ''.join(['A']+['N']*(len(secstruct)-1))

	assert(len(pattern) == len(secstruct))

	sequences = get_sequences_for_pattern(pattern)
	sequences = filter_by_secstruct(sequences, secstruct)

	p_target, x, d, p = [], [], [], []

	for sequence in sequences:
		p_target_i, x_i, d_i, p_i = test_design(sequence, secstruct, params)
		p_target.append(p_target_i)
		x.append(x_i)
		d.append(d_i)
		p.append(p_i)

	idx = np.argsort(-1*p_target)
	p_target, x, d, p, sequences = p_target[idx], x[idx], d[idx], p[idx], sequences[idx]

	subplot(1,2,1)
	imshow(get_bpp(x[0], d[0], p[0], params))
	title('Best design\n %s \n p = %.3f' % (sequences[0], p[0]))
	subplot(1,2,2)
	
	imshow(get_bpp(x[-1], d[-1], p[-1], params ))
	title('Best design\n %s \n p = %.3f' % (sequences[-1], p[-1]))

	return sequences, p_target, x, d, p