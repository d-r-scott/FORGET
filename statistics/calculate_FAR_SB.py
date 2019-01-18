#!/usr/bin/env python

import numpy as np
import csv
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from grouping import external_grouping

SB_path = '../../schedblocks/'
SB_info_path = SB_path + 'SB_info.csv'

smp_path = 'smp_FARs/'
grp_path = 'grp_FARs/'

min_path = 'min_FARs/'

#Format: list of sublists. sublist[0] is obs .cand path relative to SB_path, sublist[1] is duration
obs = []
opath = 0
odur = 1

# Index values for various candidate fields
sn  = 0	 # S/N
t   = 1	 # sampno
s   = 2	 # secs from file start
w   = 3	 # boxcar width
idt = 4	 # number of samples difference between top and bottom of frequency range
dm  = 5	 # DM
bno = 6	 # Beam number
mjd = 7	 # Modified Julian Date

lbl = 8	 # Label
ning= 9	 # Number of candidates in the group this candidate represens
rsql= 10	# R^2 for cands with time and dm less than brightest in group
rsqm= 11	# R^2 for cands with time and dm more than brightest in group

# Number of candidates to take at a time
n_cands = 40

# Simple filter parameters
dmmin = 100
wmax = 10

def _main():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Calculate False Activation Rate (FAR) for each observation inside ../../schedblocks')
	parser.add_argument('sigma', type=float, help='Sigma threshold [Required]')
	parser.add_argument('rsqmmin', type=float, help='R^2 (more) threshold [Required]')
	parser.add_argument('-m', '--minimal', action='store_true', help='Minimal output to a single output file inside ' + min_path)
	args = parser.parse_args()

	snmin = args.sigma
	rsqmmin = args.rsqmmin
	
	# Get observation paths and durations
	with open(SB_info_path, 'r') as f:
		reader = csv.reader(f, delimiter=',')
		for row in reader:
			if row[1] != '-1.0':	# A duration of -1 indicates inability to determine duration
				new_obs = [SB_path + row[0], float(row[1])]
				obs.append(new_obs)

	header = '# num, duration (s), minimum rate (s^-1), maximum rate (s^-1)\n'
	intf = '%d '
	floatf = '%0.5f '
	formats = intf + floatf + floatf + floatf + '\n'

	# Iterate over all observations, calculate FAR for both simple cutoff and grouping by processing
	# in groups of n_cands
	for o in obs:
		sbid_obsid = o[opath][len(SB_path):len(SB_path)+22]
		print(sbid_obsid)

		#len(SB_path):len(SB_path)+22 makes the outfile names simply SBID/OBSID....csv
		smp_outfile = smp_path + sbid_obsid + '_s' + str(snmin) +  '.csv'
		grp_outfile = grp_path + sbid_obsid + '_s' + str(snmin) + '_r' + str(rsqmmin) + '.csv'

		# Total number of candidates post-grouping
		n_smp = 0
		n_grp = 0

		# Time of last candidate (seconds after file start)
		t_last = 0.0

		with open(o[opath], 'r') as f:
			total_cands = file_len(o[opath]) - 1	#-1 for header
			cands = []
			count = 0
			ii = 0
			while True:
				for i in range(n_cands):
					line = f.readline()
					if len(line) > 1 and line[0] != '#':
						new_cand = map(float, line.split(' ')[0:8])
						t_last = new_cand[s] if new_cand[s] > t_last else t_last
						cands.append(new_cand)

				if len(cands) < 5:
					break

				# if not os.path.isfile(outfile): check means we don't redo something that's already been done
				# Simple filtering 
				smp_cands = []
				if not os.path.isfile(smp_outfile):
					smp_cands = simple_filter(cands, snmin)

				n_smp += len(smp_cands)

				# Grouping
				grp_cands = []
				if not os.path.isfile(grp_outfile):
					grp_cands = external_grouping(cands, dmmin, wmax, snmin, True, rsqmmin)

				n_grp += len(grp_cands)

				#Output counter
				count += len(cands)
				ii += 1
				if ii % 10000 == 0:
					print("#%d/%d - %0.1f%%" % (count, total_cands, float(count)/float(total_cands)*100))

				#Skip
				#for j in range(3960):
				#	line = f.readline()

				cands = []


		min_smp_rate = float(n_smp)/o[odur]
		max_smp_rate = float(n_smp)/t_last

		min_grp_rate = float(n_grp)/o[odur]
		max_grp_rate = float(n_grp)/t_last

		smp_write = formats % (n_smp, o[odur], min_smp_rate, max_smp_rate)
		grp_write = formats % (n_grp, o[odur], min_grp_rate, max_grp_rate)

		print(header)
		print(smp_write)
		print(grp_write)

		if args.minimal:
			min_outfile = '%ss%0.3f_r%0.3f' % (min_path, snmin, rsqmmin)
			with open(min_outfile, 'a') as mf:
				mf.write(sbid_obsid + '\n')
				mf.write(smp_write)
				mf.write(grp_write)
		else:
			if not os.path.isfile(smp_outfile):
				with open(smp_outfile, 'w') as sf:
					sf.write(header)
					sf.write(smp_write)
			else:
				print("# Skipped simple")

			if not os.path.isfile(grp_outfile):
				with open(grp_outfile, 'w') as gf:
					gf.write(header)
					gf.write(grp_write)
			else:
				print("# Skipped grouping")

def simple_filter(cands, snmin):
    sfil_cands = []
    for cand in cands:
        if cand[dm] >= dmmin and cand[w] <= wmax and cand[sn] >= snmin:
            sfil_cands.append(cand)

    return sfil_cands

def file_len(fname):
	with open(fname) as f:
		retval = sum(1 for _ in f)
	return retval

if __name__ == '__main__':
	_main()
