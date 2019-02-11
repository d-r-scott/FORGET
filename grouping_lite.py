#!/usr/bin/env python
"""
grouping_lite.py

AUTHOR: David Scott [david.r.scott@graduate.curtin.edu.au]

PURPOSE: Groups candidate events in given files by simply comparing all pairs to find coincident events
         Adapted from the clustering algorithm used in Heimdall by Ben Barsdell and Andrew Jameson
             https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/label_candidate_clusters.cu
         The structure of the algorithm draws inspiration from friends-of-friends used in FREDDA by Keith Bannister
"""

import numpy as np

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

__author__ = "David Scott <david.r.scott@graduate.curtin.edu.au>"

def _main():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Group candidate events. Input is any file containing candidates, output is to <in_filename>.grouped.\n\nTo use in another python script, simply import external_grouping from this script, and call with the syntax:\nnew_cands = external_grouping(cands, dmmin, wmax, snmin, rsqmin)', formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--ttol', type=int, help='Time tolerance - how many time samples apart are coincident events?', default=3)
	parser.add_argument('-d', '--dmtol', type=float, help='DM tolerance - how many DM units (in pc cm^-3) apart are coincident events?', default=2.)
	parser.add_argument('-w', '--wtol', type=int, help='Width tolerance - how close do the widths (in number of time samples) have to be for events to be coincident?', default=2)
	parser.add_argument('--dmmin', type=float, help='Minimum DM to return (pc/cm3)', default=0.)
	parser.add_argument('--wmax', type=int, help='Maximum width to return (time samples)', default=20)
	parser.add_argument('--snmin', type=float, help='Minimum S/N to return', default=0.)
	parser.add_argument('--rsqmin', type=float, help='Minimum R^2 to return', default=-100.)

	parser.add_argument(dest='files', nargs='+')
	args = parser.parse_args()

	for fname in args.files:
		old_cands = open_file(fname, args)		# Open file and return all valid cands

		if len(old_cands) > 0:
			new_cands = group(old_cands, args)		# Group candidates together and return a list of the groups
		else:
			new_cands = []

		print "Reduced number of candidates in %s from %d to %d" % (fname, len(old_cands), len(new_cands))

		if len(new_cands) > 0:
			write_cands(fname, new_cands)

# Open the file with given filename and return an array of candidates
# 	The candidates are structured as an array with the following fields:
# 	[ S/N, sampno, secs from file start, width, idt, DM, beam number, MJD, label ]
#	These are the fields that are assumed to be in the file, with the addition of the label field
def open_file(fname):
	cands = []
	with open(fname, 'r') as f:
		for i, line in enumerate(f):
			if line[0] != '#' and len(line) > 5:
				# In case the file has more columns than we need, trim the extras off the end
				new_cand = map(float, line.split()[0:mjd])

				# Sometimes there's no mjd field by default, so we need to add it
				while len(new_cand) <= mjd:
					new_cand.append(0.0)

				# Add a label field and number of candidates in group field
				new_cand.append(0)
				new_cand.append(1)
				cands.append(new_cand)

	# Sort by S/N (descending)
	cands.sort(key=lambda x: x[sn])

	# Give each candidate a label corresponding to its rank when ordered by S/N
	for i, cand in enumerate(cands):
		cand[lbl] = i

	return cands

# Iterate over the cands and group them with other cands that are nearby in time-DM-width space
# The returned list of cands will contain only those in each group with the highest S/N
def group(cand_list, args):
	# For every cand, compare it to each other cand
	# If the cands are within the time, DM, and width tolerances of each other, give them both the label of the cand with the highest S/N
	for i in range(len(cand_list)):
		for j in range(len(cand_list)):
			if i != j:
				t_diff  = abs(cand_list[i][t]  - cand_list[j][t])
				dm_diff = abs(cand_list[i][dm] - cand_list[j][dm])
				w_diff  = abs(cand_list[i][w]  - cand_list[j][w])

				if t_diff <= args.ttol and dm_diff <= args.dmtol and w_diff <= args.wtol:
					# Give both cands the label of whichever is brightest
					if cand_list[i][sn] >= cand_list [j][sn]:
						cand_list[j][lbl] = cand_list[i][lbl]
					else:	
						# j brighter than i
						cand_list[i][lbl] = cand_list[j][lbl]

	# We've found all the groups, now trace the label chains to give each group the same label
	for cand in cand_list:
		trace_equivalency_chain(cand, cand_list)

	# Count how many cands there are in each group
	# Get the labels of each group by finding the unique labels in cand_list
	labels = [ cand[lbl] for cand in cand_list ]
	unique_labels = list(set(labels))
	for label in unique_labels:
		cand_list[label][ning] = labels.count(label)

	# Determine the correlation coefficient of each group
	det_corr_coefs(cand_list, unique_labels)

	# Return just a list of the cands that are at the top of each chain
	# Effectively, the cands whose label have NOT changed, i.e. cand_list[i][lbl] == i
	new_cands = []
	for i in range(len(cand_list)):
		if cand_list[i][lbl] == i:
			new_cands.append(cand_list[i])

	# Filter out undesired candidates
	# Even if the rsqmin values are left as default, this will filter out those with NaN as R^2 values
	# cand[rsql] >= -100 because the lower quadrant can't be relied on to always be positive, but there will always be cands there
	#	This check is effectively asking "Are there ANY cands in the lower quadrant?"
	new_cands = [ cand for cand in new_cands if cand[rsqm] >= args.rsqmin
											and cand[rsql] >= -100
	                                        and cand[sn] >= args.snmin 
	                                        and cand[dm] >= args.dmmin
	                                        and cand[w] <= args.wmax ]

	return new_cands

# Move the cand along the label chain and set the imported cand's label to the label at the end of it
# Return a list of the candidates in the chain
def trace_equivalency_chain(cand, cand_list):
	chain = []
	while True:
		chain.append(cand)
		old_lbl = cand[lbl]
		old_sn = cand[sn]
		cand[ning] = cand_list[old_lbl][ning]
		cand = cand_list[old_lbl]

		if cand[lbl] == old_lbl and cand[sn] == old_sn:		# We've reached the end of the chain
			# Set all labels to be the same
			for link in chain:
				link[lbl] = old_lbl
			break


"""
Determine the correlation coefficients for every group
There will be two lines used to calculate these coefficients. Both will start at the brightest cand in
	the group. One goes to the point at the minimum values of time and DM, and the other to the
	maximum of these values. The correlation coefficients of each line are determined and added to the
	brightest cand in each group.
These coefficients could be used to determine if the group corresponds to a real FRB.
"""
def det_corr_coefs(cand_list, unique_labels):

	for label in unique_labels:
		# Each label is the index of the brightest candidate in each group
		# Get a list containing just the cands in the current group
		group = [ cand for cand in cand_list if cand[lbl] == label ]

		# Calculate R^2 for cands with start time and DM less than group's brightest cand
		less_group = [ cand for cand in group if cand[t] <= cand_list[label][t] and cand[dm] <= cand_list[label][dm] ]
		less_ts = [ cand[t] for cand in less_group ]
		less_dms = [ cand[dm] for cand in less_group ]
		less_r_squared = calc_r_squared(np.array(less_ts), np.array(less_dms))

		# Calculate R^2 for cands with start time and DM more than group's brightest cand
		more_group = [ cand for cand in group if cand[t] >= cand_list[label][t] and cand[dm] >= cand_list[label][dm] ]
		more_ts = [ cand[t] for cand in more_group ]
		more_dms = [ cand[dm] for cand in more_group ]
		more_r_squared = calc_r_squared(np.array(more_ts), np.array(more_dms))

		# Give brightest cand the R^2 values
		cand_list[label].append(less_r_squared)
		cand_list[label].append(more_r_squared)


# Calculate the correlation coefficient (R^2) of given data
def calc_r_squared(x, y):
	with np.errstate(divide='ignore', invalid='ignore'):	# Hide warning messages about zero division
		x_min = min(x)
		x_max = max(x)
		y_min = min(y)
		y_max = max(y)

		# Calculate gradient of line going from (x_min, y_min) to (x_max, y_max) and use it as a model
		m = (y_max - y_min)/(x_max - x_min)
		f = lambda i : m*(i - x_min) + y_min

		model_y = np.array([ f(i) for i in x ])

		# Sums of squares
		ss_res = np.sum(np.square(y - model_y))	# Residuals
		ss_tot = np.sum(np.square(y - np.average(y)))	# Total

		r_squared = 1 - ss_res/ss_tot

		return r_squared

def write_cands(fname, cands):
	header = 'S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd, label, number in group, R^2'
	intf = '%d'
	floatf = '%0.2f'
	formats = (floatf, intf, floatf, intf, intf, floatf, intf, '%0.15f', intf, intf, floatf, floatf)
	npcands = np.array(cands)
	np.savetxt(fname+'.grouped', npcands, fmt=formats, header=header)

# Allows use outside of this file
def external_grouping(cands, ext_dmmin, ext_wmax, ext_snmin, ext_rsqmin):
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Group candidate events. Input is any file containing candidates, output is to <in_filename>.grouped', formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--ttol', type=int, help='Time tolerance - how many time samples apart are coincident events?', default=3)
	parser.add_argument('-d', '--dmtol', type=float, help='DM tolerance - how many DM units (in pc cm^-3) apart are coincident events?', default=2.)
	parser.add_argument('-w', '--wtol', type=int, help='Width tolerance - how close do the widths (in number of time samples) have to be for events to be coincident?', default=2)
	parser.add_argument('--dmmin', type=float, help='Minimum DM to return (pc/cm3)', default=0.)
	parser.add_argument('--wmax', type=int, help='Maximum width to return (time samples)', default=20)
	parser.add_argument('--snmin', type=float, help='Minimum S/N to return', default=0.)
	parser.add_argument('--rsqmin', type=float, help='Minimum R^2 to return', default=-100.)

	parser.add_argument(dest='files', nargs='+')

	parse_str = '--dmmin ' + str(ext_dmmin) + ' --wmax ' + str(ext_wmax) + ' --snmin ' + str(ext_snmin)
	parse_str = parse_str + ' --rsqmin ' + str(ext_rsqmin)
	parse_str = parse_str + ' not_a_real_file'

	ext_args = parser.parse_args(parse_str.split())

	cands.sort(key=lambda x: x[sn])

	# Give each candidate a label corresponding to its initial position in the list
	#for i, cand in enumerate(cands):
		#cand.append(0)
		#cand.append(1)
		#cand[lbl] = i

	return group(cands, ext_args)

if __name__ == '__main__':
	_main()
