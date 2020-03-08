import numpy as np
from conformation import Conformation
from design import Design
from utils.energy_utils import Parameters

class NascentChain(object):

	def __init__(self, secstruct=None, sequence=None, params = None):

		if params:
			self.params = params
		else:
			self.params = Parameters() # default parameter values

		if sequence is None:
			raise RuntimeError('Cotranscriptional class requires sequence.')
		else:
			self.sequence = sequence
		self.N = len(sequence)
		self.models = []
		self.evaluated=False

		self.MFE_folding_path = []

	def run(self):
		print('getting conformations for each step of nascent chain ...')
		for i in range(1,self.N+1):
		    conf = Conformation(sequence=self.sequence[:i], params = self.params)
		    conf.run()
		    self.models.append(conf)
		    self.evaluated = True
		    print(self.sequence[:i])
		    print(conf.dbn_strings[0])


	def get_MFE_folding_path(self):

		if not self.evaluated: raise RuntimeError('Not run yet, call .run()')

		for i in range(self.N):
			self.MFE_folding_path.append(self.models[i].dbn_strings[0])


