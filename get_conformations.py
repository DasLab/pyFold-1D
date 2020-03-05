import numpy as np
from utils import *
from copy import copy

class Parameters(object):

	def __init__(self, epsilon=-2, delta=1):

		self.epsilon = epsilon
		self.delta = delta

def score_bends(d):

#	Take a bunch of conformations and score the number of bends in each
# conformation. For use in determining bending energy penalty ('Delta' in
# toyfold notes).
#
# INPUT
#
# d = [Nbeads x Nconformations] input directions
#
# OUTPUT
#
# num_bends = [num_conformations] number of bends
	return np.sum((d[:-1] != d[1:]),axis=0)

def score_pairs(p):
# Take a bunch of pairing patterns and return number of pairs,
#  needed for scoring.
#
# INPUT
#
#  p = [Nbeads x Nconformations] partners  (-1 if bead is unpaired,
#        otherwise index of partner from 0,... Nbeads-1 )
#
# OUTPUT
#
# num_bends = [num_conformations] number of pairs
	return np.sum((p>=0),axis=0)/2

def get_energy(d, p, params=None):
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

	if params is None:
		params = Parameters()

	num_bends = score_bends(d)
	num_pairs = score_pairs(p)

	return params.delta * num_bends + params.epsilon * num_pairs

def get_Z(x,d,p,E, params):

# Partition function calculation.
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
# OUTPUT
# Z: partition function of system.
# conf_prob: probabilities of enumerated conformations.

	Z = np.sum(np.exp(-E))
	conf_prob = np.exp(-E)/Z

	return Z, conf_prob

def get_conformations(secstruct, sequence=None, params=None):

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

	Returns:
	x = [Nbeads x Nconformations] all sets of conformations. If no base pairs specified,
	should get 2^(Nbeads - 1). First position is always 0.

	d = [Nbeads x Nconformations] Input directions (array of +/-1's)

	p = [Nbeads x Nconformations] partners (-1 if bead is unpaired,
	 otherwise index of partner from 0, ... Nbeads-1)

	E = [Nconformations] Energies for each conformation.
	'''

	if params is None:
		params = Parameters() # default values

	is_chainbreak, secstruct_new = parse_out_chainbreak(secstruct)

	if sequence is not None:
		is_chainbreak_sequence, sequence = parse_out_chainbreak(sequence)

	if len(secstruct_new) == 0:
		if sequence is None:
			raise RuntimeError("Must input at a minimum a secstruct or a sequence")
		N = len(sequence)
		partner = -1*np.ones([N])
		stem_assignment = np.zeros([N])
		is_chainbreak = is_chainbreak_sequence

	else:
		partner = secstruct_to_partner(secstruct_new)
		stem_assignment = figure_out_stem_assignment(secstruct_new)
		N = len(secstruct_new)

	#initialize arrays with one conformation
	x = np.zeros([N,1])
	d = np.zeros([N,1])
	p = np.zeros([N,1])

	x[:] = np.NaN
	d[:] = np.NaN
	p[:] = np.NaN

	#initialize first bead location
	x[0,0] = 0 #already set
	d[0,0] = 1
	p[0,0] = partner[0]

	for i in range(1,N):

		if not is_chainbreak[i-1]:
			if stem_assignment[i] > 0 and stem_assignment[i]==stem_assignment[i-1]:
				#continuing a stem. go in same direction

				d[i] = d[i-1]
			else:
				q = x.shape[1] # number of conformations enumerated so far

				# Two choices for next move: forward or backward.
				# Forward:

				d[i,:q] = 1

				#Backward:

				#But first, expand second axis for arrays (implicit in matlab)
				x = np.hstack([x,np.zeros([N,q])])
				d = np.hstack([d,np.zeros([N,q])])
				p = np.hstack([p,-1*np.ones([N,q])])

				newblock = q + np.arange(q)

				x[:,newblock] = copy(x[:,:q])
				d[:,newblock] = copy(d[:,:q])
				p[:,newblock] = copy(p[:,:q])
				d[i, newblock] = -1
			x[i] = x[i-1]+d[i]
			p[i] = partner[i]

		else:
			q = x.shape[1]

			# New strand!
			for xx in np.arange(-2*N-1,2*N):

			# There's one more extra conformation
			# introduced here than needed, does it get filtered?
			# How to pick 2N bounds? Affects chemical potential of second strand

				for dd in [1,-1]:
					#print(xx,dd)
					newblock = x.shape[1] + np.arange(q)

					#expand second axis
					x = np.hstack([x,np.zeros([N,q])])
					d = np.hstack([d,np.zeros([N,q])])
					p = np.hstack([p,-1*np.ones([N,q])])

					x[:,newblock] = copy(x[:,:q])
					d[:,newblock] = copy(d[:,:q])
					p[:,newblock] = copy(p[:,:q])		
					
					x[i,newblock] = xx
					d[i,newblock] = dd
					p[i,newblock] = partner[i]

					# print('x:', x)
					# print('d:', d)
					# print('p:', p)

		if np.sum(stem_assignment) > 0:
			# Filter trajectories that obey user-inputted pairs: positions at same level
			#going in opposite directions
			if partner[i] == -1:
				continue
			elif partner[i] > i:
				continue
			else:
				# print(x)
				# print(d)
				# print(x[i], x[int(partner[i]-1)])
				# print(d[i], d[int(partner[i]-1)])

				gp_x = np.where(x[i] == x[int(partner[i])])
				gp_d = np.where(d[i] == -1*d[int(partner[i])])
				gp = np.intersect1d(gp_x, gp_d)
				# print(gp)

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

				if sequence[i] == 'N':
					sequence_match = [x for x in range(i)]
				elif sequence[i] == 'A':
					sequence_match = find_all(sequence[:i],['N','U'])
				elif sequence[i] == 'U':
					sequence_match = find_all(sequence[:i],['N','A'])
				elif sequence[i] == 'G':
					sequence_match = find_all(sequence[:i],['N','C'])
				elif sequence[i] == 'C':
					sequence_match = find_all(sequence[:i],['N','G'])
				elif sequence[i].islower():
					sequence_match = find_all(sequence[:i], sequence[i])

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

			x = xnew
			d = dnew
			p = pnew

	E = get_energy(d,p,params)

	idx = np.argsort(E)

	x = x[:,idx]
	d = d[:,idx]
	p = p[:,idx]
	E = E[idx]

	return x, d, p, E

