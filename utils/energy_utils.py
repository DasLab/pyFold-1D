import numpy as np

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
