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

smp = 0
grp = 1

#Skips based on a run at snmin=7 rsqmmin=0. If an obs returned 0 candidates, all higher values will too
skips = [[7450, 190101212657, [False, False]],
[7451, 190101213644, [True, True]],
[7451, 190101213754, [True, True]],
[7451, 190101223627, [True, True]],
[7451, 190101233640, [True, True]],
[7451, 190102001400, [True, True]],
[7451, 190102003653, [True, True]],
[7451, 190102013706, [True, True]],
[7451, 190102023729, [True, True]],
[7453, 190102033451, [True, True]],
[7453, 190102034238, [True, True]],
[7453, 190102034837, [True, True]],
[7453, 190102043215, [True, True]],
[7453, 190102051443, [True, True]],
[7457, 190102070152, [False, False]],
[7457, 190102071030, [False, True]],
[7457, 190102071159, [False, False]],
[7457, 190102071319, [False, True]],
[7457, 190102071449, [False, False]],
[7471, 190102120527, [False, False]],
[7601, 190112105100, [True, True]],
[7601, 190112105210, [False, False]],
[7601, 190112115103, [False, False]],
[7601, 190112121656, [False, False]],
[7601, 190112122504, [False, False]],
[7601, 190112125116, [False, True]],
[7601, 190112212911, [False, False]],
[7601, 190112213948, [False, False]],
[7601, 190112214735, [False, False]],
[7601, 190112215314, [True, True]],
[7602, 190112232249, [False, False]],
[7602, 190113001534, [False, False]],
[7602, 190113010211, [False, False]],
[7602, 190113011547, [False, False]],
[7602, 190113020543, [True, True]],
[7603, 190113022637, [False, False]],
[7606, 190113041446, [True, True]],
[7606, 190113041825, [False, False]],
[7606, 190113042154, [False, False]],
[7606, 190113042623, [False, False]],
[7606, 190113043142, [False, False]],
[7606, 190113043730, [False, False]],
[7607, 190113044657, [False, False]],
[7607, 190113045315, [True, True]],
[7607, 190113045515, [False, True]],
[7607, 190113050004, [False, False]],
[7610, 190113053424, [False, False]],
[7610, 190113053803, [False, False]],
[7610, 190113060715, [False, False]],
[7610, 190113064454, [False, False]],
[7610, 190113085737, [False, False]],
[7610, 190113093456, [False, False]],
[7610, 190113104317, [False, False]],
[7610, 190113112645, [False, False]],
[7611, 190113114330, [False, False]],
[7611, 190113124104, [False, False]],
[7611, 190113133040, [False, False]],
[7611, 190113134107, [False, True]],
[7611, 190113144100, [False, False]],
[7611, 190113154104, [False, False]],
[7611, 190113164106, [False, False]],
[7611, 190113174109, [False, False]],
[7611, 190113184113, [False, True]],
[7611, 190113204004, [False, False]],
[7611, 190113204114, [False, False]],
[7612, 190113220909, [True, True]],
[7612, 190113222415, [True, True]],
[7612, 190113223611, [False, True]],
[7612, 190113230134, [False, False]],
[7612, 190113230932, [True, True]],
[7612, 190113235350, [False, False]],
[7612, 190114000945, [False, False]],
[7612, 190114005103, [True, True]],
[7612, 190114005223, [False, False]],
[7612, 190114010350, [False, True]],
[7612, 190114010510, [False, False]],
[7613, 190114011427, [False, False]],
[7613, 190114024631, [False, False]],
[7614, 190114032520, [False, False]],
[7615, 190114033557, [True, True]],
[7615, 190114033707, [False, False]],
[7615, 190114034016, [False, False]],
[7625, 190114120941, [True, True]],
[7625, 190114122456, [False, False]],
[7625, 190114125159, [True, True]],
[7625, 190114125318, [False, False]],
[7625, 190114130345, [False, False]],
[7625, 190114130705, [False, False]],
[7625, 190114140708, [False, False]],
[7625, 190114150711, [False, False]],
[7625, 190114160714, [False, False]],
[7625, 190114192029, [False, False]],
[7625, 190114200716, [False, False]],
[7632, 190115110751, [False, False]],
[7633, 190115112117, [True, True]],
[7633, 190115112914, [False, False]],
[7633, 190115122120, [False, False]],
[7633, 190115132123, [False, False]],
[7633, 190115142126, [False, False]],
[7633, 190115152129, [False, False]],
[7633, 190115162132, [False, False]],
[7633, 190115172125, [False, False]],
[7633, 190115182128, [False, False]],
[7633, 190115192141, [False, False]],
[7633, 190115202135, [False, False]],
[7633, 190115212138, [False, False]],
[7636, 190116003343, [False, False]],
[7636, 190116005557, [True, True]],
[7636, 190116005707, [True, True]],
[7636, 190116005817, [True, True]],
[7636, 190116005936, [False, False]],
[7637, 190116010714, [False, False]]]


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
		t_last = 0.0001

		obsid = int(sbid_obsid[-12:])
		skip_this = []
		for sk in skips:
			if sk[1] == obsid:
				skip_this = sk[2]

		# If we're skipping both simple filtering and grouping, then don't even bother opening the file!
		if len(skip_this) == 0 or (not skip_this[smp] and not skip_this[grp]):
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

					# If there's nothing in the group above our sigma limit, don't even bother trying
					if max([cand[sn] for cand in cands]) >= snmin:
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
