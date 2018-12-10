#!/usr/bin/env python
"""
grouping.py

AUTHOR: David Scott [david.r.scott@student.curtin.edu.au]

PURPOSE: Groups candidate events in given files by simply comparing all pairs to find coincident events
         Adapted from the clustering algorithm used in Heimdall by Ben Barsdell and Andrew Jameson
             https://sourceforge.net/p/heimdall-astro/code/ci/master/tree/Pipeline/label_candidate_clusters.cu
         The structure of the algorithm draws inspiration from friends-of-friends used in FREDDA by Keith Bannister
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

__author__ = "David Scott <david.r.scott@student.curtin.edu.au>"

def _main():
	# For consistency many arguments are the same as for friends of friends
	# Notable exception: -w is the width tolerance, not the maximum width considered
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Group candidate events. Input is any file containing candidates, output is to <in_filename>.grouped', formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--ttol', type=int, help='Time tolerance - how many time samples apart are coincident events?', default=32)
	parser.add_argument('-d', '--dmtol', type=float, help='DM tolerance - how many DM units (in pc cm^-3) apart are coincident events?', default=20)
	parser.add_argument('-w', '--wtol', type=int, help='Width tolerance - how close do the widths (in number of time samples) have to be for events to be coincident?', default=5)
	parser.add_argument('--tmin', type=int, help='Earliest time sample to consider', default=0)
	parser.add_argument('--dmmin', type=float, help='Minimum DM to consider (pc/cm3)', default=0.)
	parser.add_argument('--wmax', type=int, help='Maximum width to consider (time samples)', default=32)
	parser.add_argument('--snmin', type=float, help='Minimum S/N to consider', default=0.)
	parser.add_argument('-p', '--plot', action='store_true', help='Create plots', default=False)

	parser.add_argument(dest='files', nargs='+')
	args = parser.parse_args()

	for fname in args.files:
		old_cands = open_file(fname, args)		# Open file and return all valid cands
		new_cands = group(old_cands, args)		# Group candidates together and return a list of the groups

		if args.plot:
			plot_cands(old_cands, new_cands, args)

		print "Reduced number of candidates in %s from %d to %d" % (fname, len(old_cands), len(new_cands))
		write_cands(fname, new_cands)			# Write the grouped candidates to file in the same format

# Open the file with given filename and return an array of candidates
# 	The candidates are structured as an array with the following fields:
# 	[ S/N, sampno, secs from file start, width, idt, DM, beam number, MJD, label ]
#	These are the fields that are assumed to be in the file, with the addition of the label field
def open_file(fname, args):
	cands = []
	with open(fname, 'r') as f:
		for i, line in enumerate(f):
			if line[0] != '#' and line.strip() != '':
				# Not a comment and not empty
				new_cand = map(float, line.split())

				# Filter out candidates with values we want to exclude
				if new_cand[sn] >= args.snmin and new_cand[t]  >= args.tmin and new_cand[w]  <= args.wmax and new_cand[dm] >= args.dmmin:
					# Add a label field and then add the cand to the list
					new_cand.append(0)
					cands.append(new_cand)

	# Give each candidate a label corresponding to its initial position in the list
	# Sorting is no longer required
	#cands.sort(key=lambda x: x[sn])
	for i, cand in enumerate(cands):
		cand[lbl] = i

	return cands

# Iterate over the cands and group them with other cands that are nearby in time-DM-width space
# The returned list of cands will contain only those in each group with the highest S/N
def group(cand_list, args):
	# For every cand, compare it to each other cand
	# If the cands are within the time, DM, and width tolerances of each other, give them both the label of the cand with the highest S/N
	for i in range(len(cand_list)):
		for j in range(i+1, len(cand_list)):	# Start at i+1 to avoid doing things twice
			t_diff = abs(cand_list[i][t] - cand_list[j][t])
			dm_diff = abs(cand_list[i][dm] - cand_list[j][dm])
			w_diff = abs(cand_list[i][w] - cand_list[j][w])

			# Check tolerances all in one boolean to allow for short circuit evaluation if any fail
			if t_diff <= args.ttol and dm_diff <= args.dmtol and w_diff <= args.wtol:
				# Give both cands the label of whichever is brightest
				if cand_list[i][sn] >= cand_list [j][sn]:
					cand_list[j][lbl] = cand_list[i][lbl]
				else:	# j brighter than i
					cand_list[i][lbl] = cand_list[j][lbl]

	# We've found all the groups, now trace the label chains to give each group the same label
	for cand in cand_list:
		trace_equivalency_chain(cand, cand_list)

	# Return just a list of the cands that are at the top of each chain
	# Effectively, the cands whose label have NOT changed, i.e. cand_list[i][lbl] == i
	new_cands = []
	for i in range(len(cand_list)):
		if cand_list[i][lbl] == i:
			new_cands.append(cand_list[i])

	return new_cands

# Move the cand along the label chain and set the imported cand's label to the label at the end of it
def trace_equivalency_chain(cand, cand_list):
	while True:
		old_lbl = cand[lbl]
		old_sn = cand[sn]

		cand = cand_list[old_lbl]

		if cand[lbl] == old_lbl and cand[sn] == old_sn:		# We've reached the end of the chain
			break

def plot_cands(old_cands, new_cands, args):
	fig = plt.figure(figsize=(18,5))

	# Before grouping
	ax = fig.add_subplot(131, projection='3d')

	old_t = [ cand[t] for cand in old_cands ]
	old_dm = [ cand[dm] for cand in old_cands ]
	old_w = [ cand[w] for cand in old_cands ]
	old_sn = [ cand[sn] for cand in old_cands ]

	ax.scatter(old_t, old_dm, old_w, c=old_sn, cmap='winter')

	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("DM (pc cm^-3)")
	ax.set_zlabel("Width (number of samples)")

	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_dm)-args.dmtol, max(old_dm)+args.dmtol])
	ax.set_zlim([min(old_w)-args.wtol, max(old_w)+args.wtol])

	# Before grouping with error bars
	ax = fig.add_subplot(132, projection='3d')

	ax.scatter(old_t, old_dm, old_w, c=old_sn, cmap='winter')

	# Error bars
	for i in range(len(old_t)):
		ax.plot([old_t[i]+args.ttol, old_t[i]-args.ttol], [old_dm[i], old_dm[i]], [old_w[i], old_w[i]], marker="_", c="black", alpha=0.1)
		ax.plot([old_t[i], old_t[i]], [old_dm[i]+args.dmtol, old_dm[i]-args.dmtol], [old_w[i], old_w[i]], marker="_", c="black", alpha=0.1)
		ax.plot([old_t[i], old_t[i]], [old_dm[i], old_dm[i]], [old_w[i]+args.wtol, old_w[i]-args.wtol], marker="_", c="black", alpha=0.1)

	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("DM (pc cm^-3)")
	ax.set_zlabel("Width (number of samples)")

	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_dm)-args.dmtol, max(old_dm)+args.dmtol])
	ax.set_zlim([min(old_w)-args.wtol, max(old_w)+args.wtol])

	# After grouping
	ax = fig.add_subplot(133, projection='3d')

	new_t = [ cand[t] for cand in new_cands ]
	new_dm = [ cand[dm] for cand in new_cands ]
	new_w = [ cand[w] for cand in new_cands ]
	new_sn = [ cand[sn] for cand in new_cands ]

	ax.scatter(new_t, new_dm, new_w, c=new_sn, cmap='winter')

	ax.set_xlabel("Time (number of samples)")
	ax.set_ylabel("DM (pc cm^-3)")
	ax.set_zlabel("Width (number of samples)")
	
	ax.set_xlim([min(old_t)-args.ttol, max(old_t)+args.ttol])
	ax.set_ylim([min(old_dm)-args.dmtol, max(old_dm)+args.dmtol])
	ax.set_zlim([min(old_w)-args.wtol, max(old_w)+args.wtol])

	plt.show()

def write_cands(fname, cands):
	header = '# S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd, label'
	intf = '%d'
	floatf = '%0.3f'
	formats = (floatf, intf, floatf, intf, intf, floatf, intf, '%0.15f', intf)
	npcands = np.array(cands)
	np.savetxt(fname+'.grouped', npcands, fmt=formats, header=header)

if __name__ == '__main__':
	_main()
