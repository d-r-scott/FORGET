#!/usr/bin/env python
"""
grouping.py

AUTHOR: David Scott [david.r.scott@graduate.curtin.edu.au]

PURPOSE: Groups candidate events in given files by simply comparing all pairs to find coincident events
         Adapted from the clustering algorithm used in Heimdall by Ben Barsdell and Andrew Jameson
             https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/label_candidate_clusters.cu
         The structure of the algorithm draws inspiration from friends-of-friends used in FREDDA by Keith Bannister
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.time import Time

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
	# For consistency many arguments are the same as for friends of friends
	# Notable exception: -w is the width tolerance, not the maximum width considered
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Group candidate events. Input is any file containing candidates, output is to <in_filename>.grouped', formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--ttol', type=int, help='Time tolerance - how many time samples apart are coincident events?', default=3)
	parser.add_argument('-d', '--dmtol', type=float, help='DM tolerance - how many DM units (in pc cm^-3) apart are coincident events?', default=2.)
	parser.add_argument('-w', '--wtol', type=int, help='Width tolerance - how close do the widths (in number of time samples) have to be for events to be coincident?', default=2)
	parser.add_argument('--tmin', type=int, help='Earliest time sample to consider', default=0)
	parser.add_argument('--tmax', type=int, help='Latest time sample to consider', default=9999999999)
	parser.add_argument('--dmmin', type=float, help='Minimum DM to consider (pc/cm3)', default=0.)
	parser.add_argument('--dmmax', type=float, help='Maximum DM to consider (pc/cm3)', default=10000.)
	parser.add_argument('--wmax', type=int, help='Maximum width to consider (time samples)', default=20)
	parser.add_argument('--snmin', type=float, help='Minimum S/N to consider', default=0.)
	parser.add_argument('-r', '--rsq', action='store_true', help='Enable calculation of and filtering by correlation coefficients', default=False)
	parser.add_argument('--rsqlmin', type=float, help='Minimum R^2 (less) to consider', default=-100.)
	parser.add_argument('--rsqmmin', type=float, help='Minimum R^2 (more) to consider', default=-100.)
	parser.add_argument('-p', '--plot', action='store_true', help='Create plots', default=False)
	parser.add_argument('-l', '--latency', action='store_true', help='Measure and output latency', default=False)
	parser.add_argument('-s', '--stats', action='store_true', help='Generate and print extra stats with a more concise output', default=False)

	parser.add_argument(dest='files', nargs='+')
	args = parser.parse_args()

	for fname in args.files:
		# Global variables for extra statistics
		global t_first, t_last
		t_first = 999999999		# Time index of earliest cand in file (not necessarily valid cand)
		t_last = 0				# Time index of latest cand in file (not necessarily valid cand)
		old_cands = open_file(fname, args)		# Open file and return all valid cands
		if args.latency:
			start_time = Time.now()
		if len(old_cands) > 0:
			new_cands = group(old_cands, args)		# Group candidates together and return a list of the groups
			new_cands = [ cand for cand in new_cands if cand[sn] >= args.snmin ]
		else:
			new_cands = []

		if args.latency:
			end_time = Time.now()

		if args.plot and len(old_cands) > 0:
			plot_cands(old_cands, new_cands, args)

		if args.stats:
			print "#file n_bef n_aft t_first t_last"
			print "%s %d %d %d %d" % (fname, len(old_cands), len(new_cands), t_first, t_last)
		else:
			print "Reduced number of candidates in %s from %d to %d" % (fname, len(old_cands), len(new_cands))
		if args.latency:
			print "Latency: {} ms".format((end_time.mjd-start_time.mjd)*86400.0*1e3)

		#if len(new_cands) > 0:
			#write_cands(fname, new_cands)			# Write the grouped candidates to file in the same format

# Open the file with given filename and return an array of candidates
# 	The candidates are structured as an array with the following fields:
# 	[ S/N, sampno, secs from file start, width, idt, DM, beam number, MJD, label ]
#	These are the fields that are assumed to be in the file, with the addition of the label field
def open_file(fname, args):
	cands = []
	with open(fname, 'r') as f:
		for i, line in enumerate(f):
			if line[0] != '#' and len(line) > 5:
				# In case the file has more columns than we need
				new_cand = map(float, line.split()[0:7])

				# Calculate extra stats
				if args.stats:
					global t_first, t_last
					t_first = new_cand[t] if new_cand[t] < t_first else t_first
					t_last = new_cand[t] if new_cand[t] > t_last else t_last

				# Filter out candidates with values we want to exclude
				if new_cand[t] >= args.tmin and new_cand[t] <= args.tmax and new_cand[w]  <= args.wmax and new_cand[dm] >= args.dmmin and new_cand[dm] <= args.dmmax:
					# Sometimes there's no mjd field by default, so we need to add it
					while len(new_cand) <= mjd:
						new_cand.append(0.0)
					# Add a label field and number candidates in group field
					new_cand.append(0)
					new_cand.append(1)
					cands.append(new_cand)

	cands.sort(key=lambda x: x[sn])

	# Give each candidate a label corresponding to its initial position in the list
	for i, cand in enumerate(cands):
		cand[lbl] = i

	return cands

# Iterate over the cands and group them with other cands that are nearby in time-DM-width space
# The returned list of cands will contain only those in each group with the highest S/N
def group(cand_list, args):
	#cand_list.sort(key=lambda x : x[t])
	# For every cand, compare it to each other cand
	# If the cands are within the time, DM, and width tolerances of each other, give them both the label of the cand with the highest S/N
	for i in range(len(cand_list)):
		for j in range(len(cand_list)):
			if i != j:
				t_diff = abs(cand_list[i][t] - cand_list[j][t])
				dm_diff = abs(cand_list[i][dm] - cand_list[j][dm])
				w_diff = abs(cand_list[i][w] - cand_list[j][w])

				if t_diff <= args.ttol and dm_diff <= args.dmtol and w_diff <= args.wtol:
					# Give both cands the label of whichever is brightest
					if cand_list[i][sn] >= cand_list [j][sn]:
						cand_list[j][lbl] = cand_list[i][lbl]
					else:	# j brighter than i
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

	# Determine the correlation coefficients of each group
	if args.rsq:
		det_corr_coefs(cand_list, unique_labels)
	else:
		for cand in cand_list:
			# Not calculating the R^2 values, but need something to put in the column
			cand.append(0.)
			cand.append(0.)

	# Return just a list of the cands that are at the top of each chain
	# Effectively, the cands whose label have NOT changed, i.e. cand_list[i][lbl] == i
	new_cands = []
	for i in range(len(cand_list)):
		if cand_list[i][lbl] == i:
			new_cands.append(cand_list[i])

	# Filter by the R^2 values
	# Even if the rsqXmin values are left as default, this will filter out those with NaN as R^2 values
	new_cands = [ cand for cand in new_cands if cand[rsql] >= args.rsqlmin and cand[rsqm] >= args.rsqmmin ] 

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

def plot_cands(old_cands, new_cands, args):
	plt.style.use('dark_background')
	colourmap='autumn'
	marker_size = 5

	fig = plt.figure(figsize=(18, 24))

	old_t = [ cand[t] for cand in old_cands ]
	old_dm = [ cand[dm] for cand in old_cands ]
	old_w = [ cand[w] for cand in old_cands ]
	old_sn = [ cand[sn] for cand in old_cands ]

	new_t = [ cand[t] for cand in new_cands ]
	new_dm = [ cand[dm] for cand in new_cands ]
	new_w = [ cand[w] for cand in new_cands ]
	new_sn = [ cand[sn] for cand in new_cands ]

	# Before grouping
	ax = fig.add_subplot(431, projection='3d')

	ax.scatter(old_t, old_dm, old_w, c=old_sn, cmap=colourmap)
	ax.scatter(new_t, new_dm, new_w, marker='o', facecolors='none', edgecolors='white', s=marker_size*10)

	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("DM (pc cm^-3)")
	ax.set_zlabel("Width (number of samples)")

	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_dm)-args.dmtol, max(old_dm)+args.dmtol])
	ax.set_zlim([min(old_w)-args.wtol, max(old_w)+args.wtol])

	# Before grouping with error bars
	ax = fig.add_subplot(432, projection='3d')

	ax.scatter(old_t, old_dm, old_w, c=old_sn, cmap=colourmap)

	# Error bars
	for i in range(len(old_t)):
		ax.plot([old_t[i]+args.ttol, old_t[i]-args.ttol], [old_dm[i], old_dm[i]], [old_w[i], old_w[i]], marker="_", c='white', alpha=0.1)
		ax.plot([old_t[i], old_t[i]], [old_dm[i]+args.dmtol, old_dm[i]-args.dmtol], [old_w[i], old_w[i]], marker="_", c='white', alpha=0.1)
		ax.plot([old_t[i], old_t[i]], [old_dm[i], old_dm[i]], [old_w[i]+args.wtol, old_w[i]-args.wtol], marker="_", c='white', alpha=0.1)

	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("DM (pc cm^-3)")
	ax.set_zlabel("Width (number of samples)")

	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_dm)-args.dmtol, max(old_dm)+args.dmtol])
	ax.set_zlim([min(old_w)-args.wtol, max(old_w)+args.wtol])

	# After grouping
	ax = fig.add_subplot(433, projection='3d')

	ax.scatter(new_t, new_dm, new_w, c=new_sn, cmap=colourmap)

	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("DM (pc cm^-3)")
	ax.set_zlabel("Width (number of samples)")
	
	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_dm)-args.dmtol, max(old_dm)+args.dmtol])
	ax.set_zlim([min(old_w)-args.wtol, max(old_w)+args.wtol])

	# 2D plots
	# Before grouping
	ax = fig.add_subplot(434)
	ax.scatter(old_t, old_dm, c=old_sn, cmap=colourmap, s=marker_size)
	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("DM (pc cm^-3)")
	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_dm)-args.dmtol, max(old_dm)+args.dmtol])
	ax.scatter(new_t, new_dm, marker='o', facecolors='none', edgecolors='white', s=marker_size*10)

	ax = fig.add_subplot(437)
	ax.scatter(old_t, old_w, c=old_sn, cmap=colourmap, s=marker_size)
	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("Width (number of samples)")
	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_w)-args.dmtol, max(old_w)+args.dmtol])
	ax.scatter(new_t, new_w, marker='o', facecolors='none', edgecolors='white', s=marker_size*10)
	
	ax = fig.add_subplot(4,3,10)
	ax.scatter(old_dm, old_w, c=old_sn, cmap=colourmap, s=marker_size)
	ax.set_xlabel("DM (pc cm^-3)")
	ax.set_ylabel("Width (number of samples)")
	ax.set_xlim([min(old_dm)-args.ttol, max(old_dm)+args.ttol])
	ax.set_ylim([min(old_w)-args.dmtol, max(old_w)+args.dmtol])
	ax.scatter(new_dm, new_w, marker='o', facecolors='none', edgecolors='white', s=marker_size*10)

	# Before grouping with error bars
	ax = fig.add_subplot(435)
	ax.errorbar(old_t, old_dm, xerr=args.ttol, yerr=args.dmtol, fmt='o', c='white', ms=2, alpha=0.25)
	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("DM (pc cm^-3)")
	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_dm)-args.dmtol, max(old_dm)+args.dmtol])

	ax = fig.add_subplot(438)
	ax.errorbar(old_t, old_w, xerr=args.ttol, yerr=args.wtol, fmt='o', c='white', ms=2, alpha=0.25)
	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("Width (number of samples)")
	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_w)-args.dmtol, max(old_w)+args.dmtol])
	
	ax = fig.add_subplot(4,3,11)
	ax.errorbar(old_dm, old_w, xerr=args.dmtol, yerr=args.wtol, fmt='o', c='white', ms=2, alpha=0.25)
	ax.set_xlabel("DM (pc cm^-3)")
	ax.set_ylabel("Width (number of samples)")
	ax.set_xlim([min(old_dm)-args.ttol, max(old_dm)+args.ttol])
	ax.set_ylim([min(old_w)-args.dmtol, max(old_w)+args.dmtol])

	# After grouping
	ax = fig.add_subplot(436)
	ax.scatter(new_t, new_dm, c=new_sn, cmap=colourmap, s=marker_size)
	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("DM (pc cm^-3)")
	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_dm)-args.dmtol, max(old_dm)+args.dmtol])

	ax = fig.add_subplot(439)
	ax.scatter(new_t, new_w, c=new_sn, cmap=colourmap, s=marker_size)
	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("Width (number of samples)")
	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_w)-args.dmtol, max(old_w)+args.dmtol])
	
	ax = fig.add_subplot(4,3,12)
	ax.scatter(new_dm, new_w, c=new_sn, cmap=colourmap, s=marker_size)
	ax.set_xlabel("DM (pc cm^-3)")
	ax.set_ylabel("Width (number of samples)")
	ax.set_xlim([min(old_dm)-args.ttol, max(old_dm)+args.ttol])
	ax.set_ylim([min(old_w)-args.dmtol, max(old_w)+args.dmtol])

	plt.savefig("grouping.png")
	plt.tight_layout()
	plt.show()

def write_cands(fname, cands):
	header = 'S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd, label, number in group, R^2 (less), R^2 (more)'
	intf = '%d'
	floatf = '%0.2f'
	formats = (floatf, intf, floatf, intf, intf, floatf, intf, '%0.15f', intf, intf, floatf, floatf)
	npcands = np.array(cands)
	np.savetxt(fname+'.grouped', npcands, fmt=formats, header=header)

# Allows use outside of this file
def external_grouping(cands, ext_dmmin, ext_wmax, ext_snmin, r_on, ext_rsqmmin):
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Group candidate events. Input is any file containing candidates, output is to <in_filename>.grouped', formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--ttol', type=int, help='Time tolerance - how many time samples apart are coincident events?', default=3)
	parser.add_argument('-d', '--dmtol', type=float, help='DM tolerance - how many DM units (in pc cm^-3) apart are coincident events?', default=2.)
	parser.add_argument('-w', '--wtol', type=int, help='Width tolerance - how close do the widths (in number of time samples) have to be for events to be coincident?', default=2)
	parser.add_argument('--tmin', type=int, help='Earliest time sample to consider', default=0)
	parser.add_argument('--tmax', type=int, help='Latest time sample to consider', default=9999999999)
	parser.add_argument('--dmmin', type=float, help='Minimum DM to consider (pc/cm3)', default=0.)
	parser.add_argument('--dmmax', type=float, help='Maximum DM to consider (pc/cm3)', default=10000.)
	parser.add_argument('--wmax', type=int, help='Maximum width to consider (time samples)', default=20)
	parser.add_argument('--snmin', type=float, help='Minimum S/N to consider', default=0.)
	parser.add_argument('-r', '--rsq', action='store_true', help='Enable calculation of and filtering by correlation coefficients', default=False)
	parser.add_argument('--rsqlmin', type=float, help='Minimum R^2 (less) to consider', default=-100.)
	parser.add_argument('--rsqmmin', type=float, help='Minimum R^2 (more) to consider', default=-100.)
	parser.add_argument('-p', '--plot', action='store_true', help='Create plots', default=False)
	parser.add_argument('-l', '--latency', action='store_true', help='Measure and output latency', default=False)
	parser.add_argument('-s', '--stats', action='store_true', help='Generate and print extra stats with a more concise output', default=False)

	parser.add_argument(dest='files', nargs='+')

	parse_str = '--dmmin ' + str(ext_dmmin) + ' --wmax ' + str(ext_wmax) + ' --snmin ' + str(ext_snmin)
	if r_on:
		parse_str = parse_str + ' -r --rsqmmin ' + str(ext_rsqmmin)

	parse_str = parse_str + ' grouping.py'

	ext_args = parser.parse_args(parse_str.split())

	cands.sort(key=lambda x: x[sn])

	# Give each candidate a label corresponding to its initial position in the list
	for i, cand in enumerate(cands):
		cand.append(0)
		cand.append(1)
		cand[lbl] = i

	return group(cands, ext_args)

if __name__ == '__main__':
	_main()
