#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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

	cands_list_list = sorted(cands_list_list, key=lambda cand_list: max(cand_list, key=lambda cand: float(cand[sn])))

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
	marker_size=10
	line_width=1

	more_gradients = []
	less_gradients = []
	more_rsqs = []
	less_rsqs = []

	# For now 16 files assumed --> 4 columns, 4 rows
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
		ax.scatter(ts, dms, c=sns, cmap=colourmap, s=marker_size)

		# Quadrants
		# More
		ax.plot([int(brightest[t]),int(brightest[t])], [float(brightest[dm]), max(more_dms)], c='white', linewidth=line_width)
		ax.plot([int(brightest[t]), max(more_ts)], [max(more_dms), max(more_dms)], c='white', linewidth=line_width)
		ax.plot([int(brightest[t]), max(more_ts)], [float(brightest[dm]), float(brightest[dm])], c='white', linewidth=line_width)
		ax.plot([max(more_ts), max(more_ts)], [float(brightest[dm]), max(more_dms)], c='white', linewidth=line_width)
		
		# Less
		ax.plot([int(brightest[t]),int(brightest[t])], [float(brightest[dm]), min(less_dms)], c='white', linewidth=line_width)
		ax.plot([int(brightest[t]), min(less_ts)], [min(less_dms), min(less_dms)], c='white', linewidth=line_width)
		ax.plot([int(brightest[t]), min(less_ts)], [float(brightest[dm]), float(brightest[dm])], c='white', linewidth=line_width)
		ax.plot([min(less_ts), min(less_ts)], [float(brightest[dm]), min(less_dms)], c='white', linewidth=line_width)

		# Fit line
		ax.plot([int(brightest[t]), max(more_ts)], [float(brightest[dm]), max(more_dms)], c='#0000FF', linewidth=2*line_width)
		ax.plot([int(brightest[t]), min(less_ts)], [float(brightest[dm]), min(less_dms)], c='#0000FF', linewidth=2*line_width)

		# Calculate values
		more_gradients.append((max(more_dms) - float(brightest[dm]))/(max(more_ts)-int(brightest[t])))
		less_gradients.append((min(less_dms) - float(brightest[dm]))/(min(less_ts)-int(brightest[t])))
		more_rsqs.append(calc_r_squared(np.array(more_ts), np.array(more_dms)))
		less_rsqs.append(calc_r_squared(np.array(less_ts), np.array(less_dms)))

		ax.set_xlabel("Time (number of samples)")
		ax.set_ylabel("DM (pc/cm3)")
		ax.set_title(args.files[i])

	plt.savefig("DM_t.png")
	#plt.show()

	df = pd.DataFrame({'Gradient (more)':more_gradients, 'Gradient (less)':less_gradients, 'R^2 (more)':more_rsqs, 'R^2 (less)':less_rsqs})
	print(df.describe())

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
		ss_res = np.sum(np.square(y - model_y)) # Residuals
		ss_tot = np.sum(np.square(y - np.average(y)))   # Total

		r_squared = 1 - ss_res/ss_tot

		return r_squared


if __name__ == "__main__":
	_main()
