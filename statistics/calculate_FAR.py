#!/usr/bin/env python

import numpy as np
from math import floor
import matplotlib.pyplot as plt

from ../grouping import external_grouping

# Read all_sorted.cand 40 lines at a time
# Push each group through a simple filter/grouping/grouping with R^2 filtering
# Convert the mjd times post-filtering into integer second values
# For every seconds value occuring, the number of candidates post-grouping is the number of false activations per minute
# Use seconds instead of e.g. hours because I don't know how long each observation ran for and therefore how much of each hour was actually covered
# Seconds are small enough where I can be reasonably sure that almost all seconds are 100% covered

# Index values for various candidate fields
sn 	= 0		# S/N
t 	= 1		# sampno
s 	= 2		# secs from file start
w 	= 3		# boxcar width
idt = 4		# number of samples difference between top and bottom of frequency range
dm  = 5		# DM
bno = 6		# Beam number
mjd = 7		# Modified Julian Date

lbl = 8		# Label
ning= 9 	# Number of candidates in the group this candidate represens
rsql= 10	# R^2 for cands with time and dm less than brightest in group
rsqm= 11	# R^2 for cands with time and dm more than brightest in group

# Number of candidates to take at a time
n_cands = 40

# Simple filter parameters
dmmin = 100
wmax = 10
#snmin = 10

def _main():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Calculate False Activation Rate (FAR) for simple filtering and grouping algorithms based on all_sorted.cand')
	parser.add_argument('sigma', type=float, help='Sigma threshold [Required]')
	parser.add_argument('rsqmmin', type=float, help='R^2 (more) threshold [Required]')
	args = parser.parse_args()

	snmin = args.sigma
	rsqmmin = args.rsqmmin

	with open('all.cand', 'r') as f:
		cands = []
		seconds_sfil = dict()
		seconds_grouped = dict()
		count = 0
		while True:
			for i in range(n_cands):
				line = f.readline()
				if len(line) > 0 and line[0] != '#' and len(line) > 5:
					new_cand = map(float, line.split('\t')[0:8])
					new_cand[mjd] = int(floor(new_cand[mjd]*86400))	# Convert mjd to seconds
					cands.append(new_cand)

			if len(cands) < 10:
				break
			#else:
			#	print(len(cands))

			# Simple filtering
			sfil_cands = simple_filter(cands, snmin)

			# Grouping (R on)
			#grouped_cands = external_grouping(cands, dmmin, wmax, snmin, True, rsqmmin)
			
			for cand in sfil_cands:
				if cand[mjd] in seconds_sfil:
					seconds_sfil[cand[mjd]] += 1
				else:
					seconds_sfil[cand[mjd]] = 1
			"""
			for cand in grouped_cands:
				if cand[mjd] in seconds_grouped:
					seconds_grouped[cand[mjd]] += 1
				else:
					seconds_grouped[cand[mjd]] = 1
			"""
			count += 1
			if count % 10 == 0:
				print(count)

			cands = []

			#for i in range(40000):
			#	line = f.readline()
		
		sfil_fname = 'FARs/FAR_simple_' + str(snmin) + '.csv'
		sfil_fwrite = open(sfil_fname, 'w')

		seconds_sfil_arr = np.array([])
		for key in seconds_sfil:
			seconds_sfil_arr = np.append(seconds_sfil_arr, float(seconds_sfil[key]))
			sfil_fwrite.write(str(key) + ',' + str(seconds_sfil[key]) + '\n')

		sfil_fwrite.close()
		"""
		grouped_fname = 'FARs/FAR_grouped_' + str(snmin) + '_' + str(rsqmmin) + '.csv'
		grouped_fwrite = open(grouped_fname, 'w')

		seconds_grouped_arr = np.array([])
		for key in seconds_grouped:
			seconds_grouped_arr = np.append(seconds_grouped_arr, float(seconds_grouped[key]))
			grouped_fwrite.write(str(key) + ',' + str(seconds_grouped[key]) + '\n')

		grouped_fwrite.close()
		"""
		print('FAR (per second, simple filter) = ' + str(np.mean(seconds_sfil_arr)))
#		print('FAR (per second, grouping)      = ' + str(np.mean(seconds_grouped_arr)))


def simple_filter(cands, snmin):
	sfil_cands = []
	for cand in cands:
		if cand[dm] >= dmmin and cand[w] <= wmax and cand[sn] >= snmin:
			sfil_cands.append(cand)

	return sfil_cands

if __name__ == '__main__':
	_main()
