#!/usr/bin/env python

import numpy as np

sn  = 0	 # S/N
t   = 1	 # sampno
s   = 2	 # secs from file start
w   = 3	 # boxcar width
idt = 4	 # number of samples difference between top and bottom of frequency range
dm  = 5	 # DM
bno = 6	 # Beam number
mjd = 7	 # Modified Julian Date


def _main():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Simple snoopy-style filtering', formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument('--tmin', type=int, help='Earliest time sample to consider', default=0)
	parser.add_argument('--tmax', type=int, help='Latest time sample to consider', default=9999999999)
	parser.add_argument('--dmmin', type=float, help='Minimum DM to consider (pc/cm3)', default=0.)
	parser.add_argument('--dmmax', type=float, help='Maximum DM to consider (pc/cm3)', default=10000.)
	parser.add_argument('--wmax', type=int, help='Maximum width to consider (time samples)', default=20)
	parser.add_argument('--snmin', type=float, help='Minimum S/N to consider', default=0.)
	parser.add_argument('-s', '--stats', action='store_true', help='Generate and print extra stats with a more concise output', default=False)

	parser.add_argument(dest='files', nargs='+')
	args = parser.parse_args()

	for fname in args.files:
		global t_first, t_last, n_bef, n_aft
		t_first = 99999999
		t_last = 0
		n_bef = 0
		n_aft = 0
		cands = open_file(fname, args)

		if args.stats:
			print "#file n_bef n_aft t_first t_last"
			print "%s %d %d %d %d" % (fname, n_bef, n_aft, t_first, t_last)

def open_file(fname, args):
	cands = []
	with open(fname, 'r') as f:
		for i, line in enumerate(f):
			if line[0] != '#' and len(line) > 5:
				# In case the file has more columns than we need
				new_cand = map(float, line.split()[0:7])

				# Calculate extra stats
				if args.stats:
					global n_bef, t_first, t_last
					n_bef += 1
					if new_cand[t] < t_first:
						t_first = new_cand[t]
					if new_cand[t] > t_last:
						t_last = new_cand[t]

				# Filter out candidates with values we want to exclude
				if new_cand[t] >= args.tmin and new_cand[t] <= args.tmax and new_cand[w]  <= args.wmax and new_cand[dm] >= args.dmmin and new_cand[dm] <= args.dmmax and new_cand[sn] >= args.snmin and int(new_cand[bno]) != 35:
					global n_aft
					n_aft += 1
					cands.append(new_cand)

	return cands

if __name__ == "__main__":
	_main()
