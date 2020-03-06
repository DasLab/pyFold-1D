import numpy as np
from utils.conformation_utils import *
from utils.design_utils import *
from utils.energy_utils import *
from copy import copy

class Conformation(object):
	'''Figure out all conformations (bead positions) that are consistent
	with a secondary structure. 

	WLOG, conformations assumed to start at origin (0)
	and initially point in positive direction.

	Conformations are returned so that first configurations
	have the most pairs and least bends.

	Input:
		secstruct (str): Secondary structure in dot-parens notation, e.g. '((.))'.
		sequence (str) (optional): Input sequence (array of 0's and 1's) with 'colors'.
			Required if you want secstruct's to be enumerated.
		params (array-like): Energy parameter values for delta, epsilon, etc.

	Attributes:
		x = [Nbeads x Nconformations] all sets of conformations. If no base pairs specified,
		should get 2^(Nbeads - 1). First position is always 0.

		d = [Nbeads x Nconformations] Input directions (array of +/-1's)

		p = [Nbeads x Nconformations] partners (-1 if bead is unpaired,
		 otherwise index of partner from 0, ... Nbeads-1)

		E = [Nconformations] Energies for each conformation.
	'''

	def __init__(self, secstruct=None, sequence=None, params=None):

		if params:
			self.params = params
		else:
			self.params = Parameters() # default parameter values

		if secstruct is not None:

			self.is_chainbreak, self.secstruct = parse_out_chainbreak(secstruct)

			self.partner = secstruct_to_partner(self.secstruct)
			self.stem_assignment = figure_out_stem_assignment(self.secstruct)

			self.N = len(self.secstruct)

			if sequence is not None:
				is_chainbreak_sequence, self.sequence = parse_out_chainbreak(sequence)

				assert (is_chainbreak_sequence == self.is_chainbreak) #breaks in design for multisequence things?

				assert(len(self.sequence) == len(self.secstruct))

		else: # no secstruct provided

			if sequence is None:
				raise RuntimeError("Must input at a minimum a secstruct or a sequence")

			else: #sequence, no secstruct
			
				self.is_chainbreak, self.sequence = parse_out_chainbreak(sequence)

				self.N = len(self.sequence)
				self.partner = -1*np.ones([self.N])
				self.stem_assignment = np.zeros([self.N])

		self.x = []
		self.d = []
		self.p = []

		self.Z = None
		self.energies = []
		self.conf_probabilities = []
		self.bpps = np.zeros([self.N, self.N])
		self.connectivity_matrices = []

	def run(self):
		self.get_conformations() # enumerate conformations
		self.score() # calculate energies, Z, probabilities

		idx = np.argsort(self.energies)
		self.x = self.x[:,idx]
		self.d = self.d[:,idx]
		self.p = self.p[:,idx]
		self.energies = self.energies[idx]
		self.conf_probabilities = self.conf_probabilities[idx]

		self.get_bpp()
		self.get_connectivity_matrices()

	def get_conformations(self):

		# doing x, d, p as local vars and setting to class objs at the end for -laziness- legibility
		x = np.zeros([self.N,1])
		d = np.zeros([self.N,1])
		p = np.zeros([self.N,1])
		x[:] = np.NaN
		d[:] = np.NaN
		p[:] = np.NaN

		#initialize first bead location
		x[0,0] = 0
		d[0,0] = 1
		p[0,0] = self.partner[0]

		for i in range(1,self.N):

			if not self.is_chainbreak[i-1]:
				if self.stem_assignment[i] > 0 and self.stem_assignment[i]==self.stem_assignment[i-1]:
					#continuing a stem. go in same direction

					d[i] = d[i-1]
				else:
					q = x.shape[1] # number of conformations enumerated so far

					# Two choices for next move: forward or backward.
					# Forward:

					d[i,:q] = 1

					#Backward:

					#But first, expand second axis for arrays (implicit in matlab)
					x = np.hstack([x,np.zeros([self.N,q])])
					d = np.hstack([d,np.zeros([self.N,q])])
					p = np.hstack([p,-1*np.ones([self.N,q])])

					newblock = q + np.arange(q)

					x[:,newblock] = copy(x[:,:q])
					d[:,newblock] = copy(d[:,:q])
					p[:,newblock] = copy(p[:,:q])
					d[i, newblock] = -1
				x[i] = x[i-1]+d[i]
				p[i] = self.partner[i]

			else:
				q = x.shape[1]

				# New strand!
				for xx in np.arange(-self.N-1, self.N+1):

				# There's one more extra conformation
				# introduced here than needed, does it get filtered?
				# How to pick bounds? Affects chemical potential of second strand

					for dd in [1,-1]:
						newblock = x.shape[1] + np.arange(q)

						#expand second axis
						x = np.hstack([x,np.zeros([self.N,q])])
						d = np.hstack([d,np.zeros([self.N,q])])
						p = np.hstack([p,-1*np.ones([self.N,q])])

						x[:,newblock] = copy(x[:,:q])
						d[:,newblock] = copy(d[:,:q])
						p[:,newblock] = copy(p[:,:q])		
						
						x[i,newblock] = xx
						d[i,newblock] = dd
						p[i,newblock] = self.partner[i]

			if np.sum(self.stem_assignment) > 0:
				# Filter trajectories that obey user-inputted pairs: positions at same level
				#going in opposite directions
				if self.partner[i] == -1:
					continue
				elif self.partner[i] > i:
					continue
				else:
					gp_x = np.where(x[i] == x[int(self.partner[i])])
					gp_d = np.where(d[i] == -1*d[int(self.partner[i])])
					gp = np.intersect1d(gp_x, gp_d)

					x = x[:,gp]
					d = d[:,gp]
					p = p[:,gp]

			else:
				# Enumerate secondary structures by figuring out possible partners of the new bead
				# with any previously positioned beads
				# Just using for loop now

				# This all is really ugly in python translation, vectorizing will be nicer

				assert(all(p[i]==-1))

				xnew = [] # hacky initialization, see below

				for m in range(x.shape[1]):

					xm = copy(x[:,m])
					dm = copy(d[:,m])
					pm = copy(p[:,m])

					#no partner for bead i

					if len(xnew)==0: # xnew actually initialized here
						xnew = xm.reshape(-1,1)
						dnew = dm.reshape(-1,1)
						pnew = pm.reshape(-1,1)

					else:
						xnew = np.hstack([xnew, xm.reshape(-1,1)])
						dnew = np.hstack([dnew, dm.reshape(-1,1)])
						pnew = np.hstack([pnew, pm.reshape(-1,1)])

					# look for partners based on matching position

					if self.sequence[i] == 'N':
						sequence_match = [x for x in range(i)]
					elif self.sequence[i] == 'A':
						sequence_match = find_all(self.sequence[:i],['N','U'])
					elif self.sequence[i] == 'U':
						sequence_match = find_all(self.sequence[:i],['N','A'])
					elif self.sequence[i] == 'G':
						sequence_match = find_all(self.sequence[:i],['N','C'])
					elif self.sequence[i] == 'C':
						sequence_match = find_all(self.sequence[:i],['N','G'])
					elif self.sequence[i].islower():
						sequence_match = find_all(self.sequence[:i], self.sequence[i])

					partners=[]
					for s in sequence_match:
						if xm[s] == xm[i] and dm[s] == -1*dm[i] and pm[s] == -1:
							partners.append(s)

					partners = list(set(partners))

					for pp in partners:
						pm = copy(p[:,m])
						pm[i] = pp
						pm[pp] = i

						xnew = np.hstack([xnew, xm.reshape(-1,1)])
						dnew = np.hstack([dnew, dm.reshape(-1,1)])
						pnew = np.hstack([pnew, pm.reshape(-1,1)])

				x, d, p = xnew, dnew, pnew

		self.x, self.d, self.p = x, d, p

	def score(self):
		# Energy in toyfold model
		#
		# Inputs
		#  d = [Nbeads x Nconformations] input directions (array of +/-1's)
		#  p = [Nbeads x Nconformations] partners  (-1 if bead is unpaired,
		#        otherwise index of partner from 0,... Nbeads-1 )
		# params = Energy parameter values for delta, epsilon, etc. [Params class]
		#
		# Output
		# E = [Nconformations] energies for each conformation 

		num_bends = score_bends(self.d)
		num_pairs = score_pairs(self.p)

		self.energies = self.params.delta * num_bends + self.params.epsilon * num_pairs
		self.Z = np.sum(np.exp(-1*self.energies))
		self.conf_probabilities = np.exp(-self.energies)/self.Z

	def get_bpp(self):
	# Get base pair probability matrix from get_conformations() output of
	#   ensemble of conformations (positions,directions,pairings).
	#
	# INPUT
	#  x = [Nbeads x Nconformations] all sets of conformations.
	#        If there are no base pairs specified, should get
	#        2^(Nbeads-1). First position is always 0.   
	#  d = [Nbeads x Nconformations] input directions (array of +/-1's)
	#  p = [Nbeads x Nconformations] partners  (0 if bead is unpaired,
	#        otherwise index of partner from 1,... Nbeads )
	#  params = Energy parameter values for delta, epsilon, etc. [MATLAB struct]
	#
	# Output
	# bpp = [Nbeads x Nbeads] matrix of probabilities (from 0 to 1) that  
	#                bead i is paired to bead j in the ensemble.

		for q in range(self.x.shape[1]):
			for m in range(self.N):
				if self.p[m,q] >= 0:
					self.bpps[m, int(self.p[m,q])] += self.conf_probabilities[q]

	def get_connectivity_matrices(self):

		'''For each conformation, get matrix that describes connectivity:
		both backbone links and base pairs.'''

		for q in range(self.x.shape[1]):

			new_mat = np.zeros([self.N, self.N])

			for m in range(self.N-1):
				if self.is_chainbreak[m]!=1:
					new_mat[m,m+1] += 1
					new_mat[m+1,m] += 1

			for m in range(self.N):
				if self.p[m,q] >= 0:
					new_mat[m, int(self.p[m,q])] += 1

		self.connectivity_matrices.append(new_mat)



