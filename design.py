from utils.design_utils import *
from utils.conformation_utils import *
from utils.energy_utils import *
import matplotlib.pyplot as plt
from conformation import Conformation

class DesignSwitch(object):

	def __init__(self, input_motifs=None, output_motifs=None, truth_table=None):
		'''
		Design a switch. Serves as a wrapper to the Design class. 
		Enumerates all the secstruct states, parameter classes, and creates
		the function to evaluate switch performance based on a couple options.

		Example: Classic Eterna FMN-MS2 switches

		`fmn_motif = Motif('(x(&)x)',dG=-1)
		ms2_motif = Motif('((x))',dG=-1)

		#Truth tables
		tt_off = {[0]:[1], [1]:[0]}
		tt_on = {[0]:[0], [1]:[1]}
		off_switch_mdl = DesignSwitch(input_motifs=[fmn_motif], output_motifs = [ms2_motif], truth_table=tt_off)
		on_switch_mdl = DesignSwitch(input_motifs=[fmn_motif], output_motifs = [ms2_motif], truth_table=tt_on)
		off_switch_mdl.run()
		on_switch_mdl.run()`

		Example: two-oligo XOR gate

		`oligo_A = Motif('(((&)))',dG=-2)
		oligo_B = Motif('(((&)))',dG=-2)
		output_stem = Motif('(((......)))',dG=-2)

		tt = {[0,0]:[0], [0,1]: [1], [1,0]:[1], [1,1]:[0]}
		xor_mdl = DesignSwitch(input_motifs=[oligo_A, oligo_B], output_motifs = [output_stem], truth_table=tt)
		xor_mdl.run()`
		'''

		# TODO: implement rest


class Design(object):

	def __init__(self, secstruct = None, pattern = None, params = None, n_switch_states=0):
		'''
	    Enumerate over all sequences that could form target secondary
	    structure (and conform to optional design pattern) and see which ones
	    fold best.

	    If designing switch: secstruct and params inputs should be a list of states (see below).
		
		Inputs
		secstruct (str): Target secondary structure in dot-parens notation, e.g.
		                '((.))'. If a list of structures, will treat the structure as a switch.
		
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
		scores  (list): For desigining single state: for each sequence, fraction of conformations with target 
		               secondary structure.
		               For switch design: output evaluating switch performance [coming].
		conformations (list of Conformation objects): Conformation objects for each sequence.

		Functions

		run()
		filter_sequence()
		score()
		test()
		test_design_()
		plot_best_and_worst()
		'''

		self.n_switch_states = n_switch_states

		if self.n_switch_states > 0:
			# switch design mode!
			if params:
				if len(params) != self.n_switch_states:
					raise RuntimeError('Error: Length of parameter list does not match n_switch_states.')
				self.params = params
			else:
				self.params = [Parameters()]*self.n_switch_states # same default parameters for each state.
																  # probably not what you want.
		
		else:
			if params:
				self.params = params
			else:
				self.params = Parameters() # default parameter values

		if secstruct is None:
			raise RuntimeError('Must provide target secondary structure.')
		else:
			if isinstance(secstruct,list):

				if len(params) != self.n_switch_states:
					raise RuntimeError('Error: length of secstruct list does not match n_switch_states.')

				# switch design mode!
				self.secstruct=[]

				for s in secstruct:
					is_chainbreak_tmp, secstruct_tmp = parse_out_chainbreak(s)
					self.secstruct.append(secstruct_tmp)
					self.is_chainbreak = is_chainbreak_tmp
			else:
				self.is_chainbreak, self.secstruct = parse_out_chainbreak(secstruct)

		if pattern is None:
			patt = ['A']
			
			self.pattern = ''.join(['A']+['N']*(len(self.secstruct)-1))
		else:
			self.pattern = pattern

		if isinstance(self.secstruct, list):
			assert(len(self.pattern) == len(self.secstruct[0]))

		else:
			assert(len(self.pattern) == len(self.secstruct))

		self.sequences = []
		self.conformations = []
		self.scores = []

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

		if self.n_switch_states == 0:
			# not a switch
			for sequence in self.sequences:
				p_target_i, base_mdl = self.test_design_nonswitch_(sequence)
				self.scores.append(p_target_i)
				self.conformations.append(base_mdl)

		else:
			# score switches
			for sequence in self.sequences:
				switch_score_i = self.test_design_switch_(sequence)
				self.scores.append(switch_score_i)

		idx = np.argsort([-1*x for x in self.scores]) # to sort highest to lowest
		self.scores = [self.scores[x] for x in idx]
		self.sequences = [self.sequences[x] for x in idx]
		self.conformations = [self.conformations[x] for x in idx]

	def test_design_nonswitch_(self, sequence):
		'''
		Determine probability that sequence folds into single target secondary structure.

		Inputs: sequence  = sequence like 'AAACCCGGA'
		'''
		base_mdl = Conformation(secstruct= None, sequence = sequence, params = self.params)
		target_mdl = Conformation(secstruct= self.secstruct, sequence = sequence, params = self.params)

		base_mdl.run()
		target_mdl.run()

		p_target = target_mdl.Z / base_mdl.Z

		return p_target, base_mdl

	def test_design_switch_(self, sequence):
		'''
		Determine performance of switch.
		'''

		Z_list = []
		for s_ind in range(self.n_switch_states):
			mdl = Conformation(secstruct= self.secstruct[s_ind], 
				sequence = sequence, params = self.params[s_ind])
			mdl.run()
			state_mdls.append(mdl)
			Z_list.append(mdl.Z)

			# TODO: score switch based on Z values in Z_list

		return score

	def plot_best_and_worst(self):
		'''
		Plot base pair probability matrix of best and worst found design.
		'''

		plt.subplot(1,2,1)
		plt.imshow(self.conformations[0].bpps,cmap='gist_heat_r')
		plt.title('Best design\n %s \n prob = %.3f' % (self.sequences[0], self.scores[0]))

		plt.subplot(1,2,2)
		plt.imshow(self.conformations[-1].bpps,cmap='gist_heat_r')
		plt.title('Worst design\n %s \n prob = %.3f' % (self.sequences[-1], self.scores[-1]))



