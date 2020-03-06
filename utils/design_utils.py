from collections import Counter
import numpy as np
from copy import copy

def complement_to_(string):
	'''from arnie'''
	base_pairing_dct = {'a':'u', 'u':'a', 'g':'c', 'c':'g','t':'a'}
	return ''.join(base_pairing_dct[x.lower()] for x in string[::-1])

def get_nts_for_symbol(s):
# Returns, e.g., 'ACGU' for 'N'

	if s in ['A','G','C','U',',','+',' ','&']:
		return [s]
	elif s=='R':
		return ['A','G']
	elif s=='Y':
		return ['C','U']
	elif s=='N':
		return ['A','C','G','U']
	elif s=='W':
		return ['A','U']
	elif s=='S':
		return ['C','G']

def check_pair(s1, s2):
# Returns 1 if s1, s2 are OK for Watson-Crick pairing.
#  E.g., 
#    check_pair( 'A', 'U' ) --> 1
#    check_pair( 'U', 'U' ) --> 0
#
#  Note that we don't have G*U in current toyfold model:
#    check_pair( 'G', 'U' ) --> 0
#
# Inputs
#  s1 = character (A,C,G,U)
#  s2 = character (A,C,G,U)
#
# Output
#  ok = are they pairable in Toyfold model?
	return (''.join([s1,s2]) in ['AU','UA','CG','GC'])

def get_sequences_for_pattern(pattern):
# get list of all possible sequences consistent with design pattern.
# E.g., pattern of AANUU  will lead to AAAUU, AACUU, AAUUU, AAGUU.
#
# Input
#   pattern = sequence like AANUU, where characters like N,R,Y,W, and S are
#                wild cards, and A, U, C, and G are preserved.
# Output
#   sequences = list of sequence strings.

	N = len(pattern)
	num_nts = []
	for j in range(N):
		num_nts.append(len(get_nts_for_symbol(pattern[j])))

	num_sequences = np.cumprod(num_nts)

	sequences=[]

	for q in range(num_sequences[-1]):
		sequence=[]
		prev_num_sequences = 1
		for j in range(N):
			nts = get_nts_for_symbol(pattern[j])
			sequence.append(nts[int(np.ceil(q/prev_num_sequences) % len(nts))])
			prev_num_sequences = num_sequences[j]
		sequences.append(''.join(sequence))

	assert( len(sequences) == num_sequences[-1]) # check we found them all
	assert( len(sequences) == len(set(sequences))) # check all are unique

	return sequences

def filter_by_secstruct(sequences, secstruct):
	# filter list of sequences to ones that are consistent with target
	# secstruct
	#
	# Inputs:
	# sequences = list of input sequences
	# secstruct = target secondary structure (dot-parens notation)
	#
	# Output:
	# ok_sequences = list of sequences, filtered for those that are consistent
	#                  with secstruct.

	bps = convert_structure_to_bps(secstruct)

	ok_sequences=[]
	ok = True

	for sequence in sequences:
		for [i, j] in bps:
			if not check_pair(sequence[i], sequence[j]):
				ok = False
				break
		if ok:
			ok_sequences.append(sequence)
	return ok_sequences

