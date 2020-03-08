from utils.design_utils import *
from utils.conformation_utils import *
from utils.energy_utils import *
import matplotlib.pyplot as plt
from conformation import Conformation

class Design(object):

	def __init__(self, secstruct = None, pattern = None, params = None):
		'''
	    Enumerate over all sequences that could form target secondary
	    structure (and conform to optional design pattern) and see which ones
	    fold best.
		
		Inputs
		secstruct (str): Target secondary structure in dot-parens notation, e.g.
		                '((.))'. If a list of structures, will treat the structure as a switch [coming].
		
		pattern (str) (optional): Sequence like AANUU, where characters like 
		               N,R,Y,W, and S are wild cards, and A,U,C, and G are preserved.
		               If not specified, code will use ANNNNNN... (note that
		               first base can be set to A without loss of generality in
		               current energy model)
		
		params (Parameters class): custom parameter values.
			If designing switches! pass a list of params classes that describe the energy models used 
			for each switch state. [coming]

		To evaluate model: call run()
		
		Attributes
		sequences (list): All sequences that match secstruct and pattern, ordered
		             with 'best' design (by p_target) first
		target_probabilities  (list): For each sequence, fraction of conformations with target 
		               secondary structure.
		conformations (list of Conformation objects): Conformation objects for each sequence.

		Functions

		run()
		filter_sequence()
		score()
		test()
		test_design_()
		plot_best_and_worst()
		'''

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
		'''
		Determine probability that sequence folds into target secondary structure.

		Inputs: sequence  = sequence like 'AAACCCGGA'
		'''

		base_mdl = Conformation(secstruct= None, sequence = sequence, params = self.params)
		target_mdl = Conformation(secstruct= self.secstruct, sequence = sequence, params = self.params)

		base_mdl.run()
		target_mdl.run()

		p_target = target_mdl.Z / base_mdl.Z

		return p_target, base_mdl

	def plot_best_and_worst(self):
		'''
		Plot base pair probability matrix of best and worst found design.
		'''

		plt.subplot(1,2,1)
		plt.imshow(self.conformations[0].bpps,cmap='gist_heat_r')
		plt.title('Best design\n %s \n prob = %.3f' % (self.sequences[0], self.target_probabilities[0]))

		plt.subplot(1,2,2)
		plt.imshow(self.conformations[-1].bpps,cmap='gist_heat_r')
		plt.title('Worst design\n %s \n prob = %.3f' % (self.sequences[-1], self.target_probabilities[-1]))



