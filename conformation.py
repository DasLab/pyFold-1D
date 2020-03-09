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
		secstruct (str): Secondary structure in dot-parens notation, e.g. '((.))'.  For unspecified nts, use '_'.
		sequence (str) (optional): Input sequence (array of 0's and 1's) with 'colors'.
		params (Parameters class): input custom parameters model (see Parameters doc for how this works.)

	Attributes:
		x: [Nbeads x Nconformations] all sets of conformations. If no base pairs specified,
		should get 2^(Nbeads - 1). First position is always 0.

		d: [Nbeads x Nconformations] Input directions (array of +/-1's)

		p: [Nbeads x Nconformations] partners (-1 if bead is unpaired,
		 otherwise index of partner from 0, ... Nbeads-1)

		energies: [Nconformations] Energies for each conformation.
		bpps: [Nbeads x Nbeads array] equilibrium-average base pair probability matrix.
		connectivity_matrices: list of arrays [Nbeads x Nbeads]: connectivity matrices used for GNM energy calculation
			(can be used down the line for visualization)
		evaluated (bool): True if enumeration sequence has been performed, False otherwise
		n_conf (int): number of conformations found
		dbn_strings: list of dot-bracket strings for each conformation found
		stems_list: list of (list of bps) (aka stems) for each conformation

	'''

	def __init__(self, secstruct=None, sequence=None, params=None):

		if params:
			self.params = params
		else:
			self.params = Parameters() # default parameter values

		if secstruct is not None:

			self.is_chainbreak, self.secstruct = parse_out_chainbreak(secstruct)
			self.N = len(self.secstruct)

			if sequence is not None:
				is_chainbreak_sequence, self.sequence = parse_out_chainbreak(sequence)

				assert (is_chainbreak_sequence == self.is_chainbreak) #breaks in design for multisequence things?

				assert(len(self.sequence) == len(self.secstruct))
			else:
				self.sequence = 'N'*self.N

		else: # no secstruct provided

			if sequence is None:
				raise RuntimeError("Must input at a minimum a secstruct or a sequence")

			else: #sequence, no secstruct
				self.is_chainbreak, self.sequence = parse_out_chainbreak(sequence)
				self.N = len(self.sequence)
				self.secstruct = '_'*self.N

		#writes self.constraints, self.stem_assignment
		self.get_starting_constraints()

		self.x = []
		self.d = []
		self.p = []

		self.Z = None
		self.energies = []
		self.conf_probabilities = []
		self.bpps = np.zeros([self.N, self.N])
		self.connectivity_matrices = []
		self.evaluated = False
		self.n_conf = 0
		self.dbn_strings = []
		self.stems_list = []

	def get_starting_constraints(self):
		''' 
		Writes self.constraints and self.stem_assignment.
	    -1 if forced to be unpaired, Nan if unspecified.
	    Otherwise self.constraints: assigned stem 1 to max(N_stems)
		Note sem_assignments not switched to zero-indexing for python version, unlike constraints syntax.
		'''

		self.constraint = -1*np.ones([self.N]) 
		self.stem_assignment = -1*np.ones([self.N])

		bps = convert_structure_to_bps(self.secstruct)
		unspecified = get_unspecified_spots(self.secstruct)
		stems = parse_stems_from_bps(bps)

		for (i,j) in bps:
			self.constraint[i] = j
			self.constraint[j] = i

		for i, stem in enumerate(stems):
			for bp in stem:
				self.stem_assignment[bp[0]] = i+1
				self.stem_assignment[bp[1]] = i+1

		for u in unspecified:
			if self.constraint[u] != -1 or self.stem_assignment[u] != -1:
				raise RuntimeError('conflicting base pairs and unspecified spot.')
			else:
				self.constraint[u] = np.NaN
				self.stem_assignment[u] = np.NaN

	def run(self):
		self.get_conformations() # enumerate conformations
		self.n_conf = self.x.shape[1]
		self.parse_conformations()
		self.get_connectivity_matrices()
		self.score() # calculate energies, Z, probabilities
		self.get_bpp()

		#sort all the things
		idx = np.argsort(self.energies)
		self.x = self.x[:,idx]
		self.d = self.d[:,idx]
		self.p = self.p[:,idx]
		self.energies = self.energies[idx]
		self.conf_probabilities = self.conf_probabilities[idx]
		self.dbn_strings = [self.dbn_strings[x] for x in idx]
		self.stems_list = [self.stems_list[x] for x in idx]
		self.connectivity_matrices = [self.connectivity_matrices[x] for x in idx]

	def get_conformations(self, debug=False):

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
		p[0,0] = -1 #self.constraint[0]

		for i in range(1,self.N):

			if debug:
				print('i:', i)
				print('x\n',x)
				print('d\n',d)
				print('p\n',p)

			# This just puts all of them out there unless continuing a stem in the constraints.
			if not self.is_chainbreak[i-1]:

				if self.stem_assignment[i] > -1 and self.stem_assignment[i]==self.stem_assignment[i-1]:
					#continuing a stem. go in same direction

					d[i] = d[i-1]
				else: # case where the next base are either unpaired or unspecified. Will be filtered later

					# this base is unpaired

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
				if np.isnan(self.constraint[i]): # assume its unpaired for now
					p[i] = -1

				else:
					p[i] = self.constraint[i]
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

						if np.isnan(self.constraint[i]): # assume its unpaired for now
							p[i,newblock] = -1

						else:
							p[i,newblock] = self.constraint[i]

			# Now we filter based on constraints either in sequence or lack thereof in constraints.
			# Filter trajectories that obey user-inputted pairs: positions at same level
			# going in opposite directions

			if debug: print(self.constraint[i])

			if self.constraint[i] == -1:
				continue

			elif self.constraint[i] > -1:

				if self.constraint[i] > i:
					continue
				else:
					# we have a constraint
					gp_x = np.where(x[i] == x[int(self.constraint[i])])
					gp_d = np.where(d[i] == -1*d[int(self.constraint[i])])
					gp = np.intersect1d(gp_x, gp_d)

					x = x[:,gp]
					d = d[:,gp]
					p = p[:,gp]

			elif np.isnan(self.constraint[i]):
				# Enumerate secondary structures by figuring out possible partners of the new bead
				# with any previously positioned beads
				# Just using for loop now

				# This all is really ugly in python translation, vectorizing will be nicer

				#assert(all(p[i]==-1))

				xnew = [] # hacky initialization, see below

				for m in range(x.shape[1]):

					if debug: print('m',m)
					xm = copy(x[:,m])
					dm = copy(d[:,m])
					pm = copy(p[:,m])

					if debug: print(pm)

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
					if debug: print('partners', partners)

					for pp in partners:
						pm = copy(p[:,m])
						pm[i] = pp
						pm[pp] = i

						xnew = np.hstack([xnew, xm.reshape(-1,1)])
						dnew = np.hstack([dnew, dm.reshape(-1,1)])
						pnew = np.hstack([pnew, pm.reshape(-1,1)])

				x, d, p = xnew, dnew, pnew

		self.x, self.d, self.p = x, d, p
		self.evaluated = True

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
		if not self.evaluated:
			raise RuntimeError('Not run yet, call .run() for pipeline or .get_conformations() for first step')

		num_bends = score_bends(self.d)
		num_pairs = score_pairs(self.p)
		num_stacks = score_stacks(self.stems_list)

		self.energies = self.params.delta * num_bends + self.params.epsilon * num_pairs + self.params.sigma * num_stacks

		if self.params.gnm:
			print('Scoring includes GNM elastic energy')
			self.energies += score_fluctuations(self.connectivity_matrices)

		if len(self.params.motifs) > 0:
			for motif in self.params.motifs:
				self.energies += motif.dG*score_motif(self.dbn_strings, motif.secstruct)

		self.Z = np.sum(np.exp(-1*self.energies))
		self.conf_probabilities = np.exp(-self.energies)/self.Z

	def parse_conformations(self):
		# To be run after `get_conformations`. Given filled-out p matrix, gets 
		# list [length n_conformation] of list [length variable] of stems.

		if not self.evaluated:
			raise RuntimeError('Not run yet, call .run() for pipeline or .get_conformations() for first step')
		for m in range(self.n_conf):
			pm = copy(self.p[:,m])
			bp_list = partner_to_bp_list(pm)
			stems = parse_stems_from_bps(bp_list)
			self.dbn_strings.append(write_dbn_from_partner(pm))
			self.stems_list.append(stems)

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

		if not self.evaluated:
			raise RuntimeError('Not run yet, call .run() for pipeline or .get_conformations() for first step')

		for q in range(self.n_conf):
			for m in range(self.N):
				if self.p[m,q] >= 0:
					self.bpps[m, int(self.p[m,q])] += self.conf_probabilities[q]

	def get_connectivity_matrices(self):

		'''For each conformation, get matrix that describes connectivity:
		both backbone links and base pairs.'''
		
		if not self.evaluated:
			raise RuntimeError('Not run yet, call .run() for pipeline or .get_conformations() for first step')

		for q in range(self.n_conf):

			new_mat = np.zeros([self.N, self.N])

			for m in range(self.N-1):
				if self.is_chainbreak[m]!=1:
					new_mat[m,m+1] += 1
					new_mat[m+1,m] += 1

			for m in range(self.N):
				if self.p[m,q] >= 0:
					new_mat[m, int(self.p[m,q])] += 1

			self.connectivity_matrices.append(new_mat)



