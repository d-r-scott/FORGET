#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Index values for various candidate fields
sn  = 0     # S/N
t   = 1     # sampno
s   = 2     # secs from file start
w   = 3     # boxcar width
idt = 4     # number of samples difference between top and bottom of frequency range
dm  = 5     # DM
bno = 6     # Beam number
mjd = 7     # Modified Julian Date

def _main():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Plot .cand files, intended for the real FRB events')
	
	parser.add_argument(dest='files', nargs='+')
	args = parser.parse_args()

	cands_list_list = []

	for fname in args.files:
		cands_list_list.append(open_file(fname, args))

	plot_all(cands_list_list, args)

def open_file(fname, args):
	cands = []
	with open(fname, 'r') as f:
		for i, line in enumerate(f):
			if line[0] != '#' and len(line) > 5:
				# In case the file has more columns than we need
				new_cand = line.split()[0:7]
				cands.append(new_cand)


	return cands

def plot_all(cands_list_list, args):
	plt.style.use('dark_background')
	colourmap='autumn'
	marker_size=3
	line_width=1

	# For now 14 files assumed --> 4 columns, 4 rows
	n_col = 4
	n_row = 4

	fig = plt.figure(figsize=(24,24))

	for i, cands_list in enumerate(cands_list_list):
		ax = fig.add_subplot(n_row, n_col, i+1)

		ts = [ int(cand[t]) for cand in cands_list ]
		dms = [ float(cand[dm]) for cand in cands_list ]
		sns = [ float(cand[sn]) for cand in cands_list ]

		brightest = max(cands_list, key=lambda cand: float(cand[sn]))

		more = [cand for cand in cands_list if cand[t] >= brightest[t] and cand[dm] >= brightest[dm]]
		more_ts = [ int(cand[t]) for cand in more ]
		more_dms = [ float(cand[dm]) for cand in more ]

		less = [cand for cand in cands_list if cand[t] <= brightest[t] and cand[dm] <= brightest[dm]]
		less_ts = [ int(cand[t]) for cand in less ]
		less_dms = [ float(cand[dm]) for cand in less ]

		# Candidates
		ax.scatter(ts, dms, c=sns, cmap=colourmap)

		# Quadrants
		# More
		ax.plot([int(brightest[t]),int(brightest[t])], [float(brightest[dm]), max(more_dms)], c='white', linewidth=line_width)
		ax.plot([int(brightest[t]), max(more_ts)], [max(more_dms), max(more_dms)], c='white', linewidth=line_width)
		ax.plot([int(brightest[t]), max(more_ts)], [float(brightest[dm]), float(brightest[dm])], c='white', linewidth=line_width)
		ax.plot([max(more_ts), max(more_ts)], [float(brightest[dm]), max(more_dms)], c='white', linewidth=line_width)
		

		ax.plot([int(brightest[t]),int(brightest[t])], [float(brightest[dm]), min(less_dms)], c='white', linewidth=line_width)
		ax.plot([int(brightest[t]), min(less_ts)], [min(less_dms), min(less_dms)], c='white', linewidth=line_width)
		ax.plot([int(brightest[t]), min(less_ts)], [float(brightest[dm]), float(brightest[dm])], c='white', linewidth=line_width)
		ax.plot([min(less_ts), min(less_ts)], [float(brightest[dm]), min(less_dms)], c='white', linewidth=line_width)


		ax.set_xlabel("Time (number of samples)")
		ax.set_ylabel("DM (pc/cm3)")
		ax.set_title(args.files[i])

	plt.savefig("real_FRB_cands.png")
	#plt.show()

if __name__ == "__main__":
	_main()
