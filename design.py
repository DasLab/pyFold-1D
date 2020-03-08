from utils.design_utils import *
from utils.conformation_utils import *
from utils.energy_utils import *
import matplotlib.pyplot as plt
from conformation import Conformation

class Design(object):

	def __init__(self, secstruct = None, pattern = None, params = None):

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

		if params:
			self.params = params
		else:
			self.params = Parameters() # default parameter values

		if secstruct is None:
			raise RuntimeError('Must provide target secondary structure.')
		else:
			self.is_chainbreak, self.secstruct = parse_out_chainbreak(secstruct) # just doing to get breaks

		if pattern is None:
			patt = ['A']
			
			self.pattern = ''.join(['A']+['N']*(len(self.secstruct)-1))
		else:
			self.pattern = pattern

		assert(len(self.pattern) == len(self.secstruct))

		self.sequences = []
		self.conformations = []
		self.target_probabilities = []

	def run(self):

		#set sequences to search over
		self.filter_sequences()
		self.test()
		if len(self.sequences) > 0:
			self.plot_best_and_worst()

	def filter_sequences(self):
		sequences = get_sequences_for_pattern(self.pattern)
		bps = convert_structure_to_bps(self.secstruct)

		for seq in sequences:
			ok=True
			for [i, j] in bps:
				if not check_pair(seq[i], seq[j]):
					ok = False
					break
			if ok:
				self.sequences.append(seq)

	def test(self):

		print('Found %d sequences to test' % len(self.sequences))

		for sequence in self.sequences:
			p_target_i, base_mdl = self.test_design_(sequence)
			self.target_probabilities.append(p_target_i)
			self.conformations.append(base_mdl)

		idx = np.argsort([-1*x for x in self.target_probabilities]) # to sort highest to lowest
		self.target_probabilities = [self.target_probabilities[x] for x in idx]
		self.sequences = [self.sequences[x] for x in idx]
		self.conformations = [self.conformations[x] for x in idx]

	def test_design_(self, sequence):
		# Main script for testing if a sequence folds well into target 
		#   secondary structure.
		#
		# Inputs
		#  sequence  = sequence like 'AAACCCGGA'
		#  secstruct = target secondary structure in dot parens notation
		#  params = Energy parameter values for delta, epsilon, etc. [MATLAB struct]
		
		# Output
		#  p_target = Fraction of conformations with target 
		#                secondary structure. (Higher is better.)
		#  x = [Nbeads x Nconformations] all sets of conformations.
		#        If there are no base pairs specified, should get
		#        2^(Nbeads-1). First position is always 0.   
		#  d = [Nbeads x Nconformations] input directions (array of +/-1's)
		#  p = [Nbeads x Nconformations] partners  (0 if bead is unpaired,
		#        otherwise index of partner from 1,... Nbeads )

		base_mdl = Conformation(secstruct= None, sequence = sequence, params = self.params)
		target_mdl = Conformation(secstruct= self.secstruct, sequence = sequence, params = self.params)

		base_mdl.run()
		target_mdl.run()

		p_target = target_mdl.Z / base_mdl.Z

		return p_target, base_mdl

	def plot_best_and_worst(self):
		plt.subplot(1,2,1)
		plt.imshow(self.conformations[0].bpps,cmap='gist_heat_r')
		plt.title('Best design\n %s \n prob = %.3f' % (self.sequences[0], self.target_probabilities[0]))

		plt.subplot(1,2,2)
		plt.imshow(self.conformations[-1].bpps,cmap='gist_heat_r')
		plt.title('Worst design\n %s \n prob = %.3f' % (self.sequences[-1], self.target_probabilities[-1]))



