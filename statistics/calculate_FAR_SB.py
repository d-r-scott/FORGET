#!/usr/bin/env python

import numpy as np
import csv
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from grouping import external_grouping

SB_path = '../../cands/'
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
skips = [['../../cands/000.cand', [False, False]],
['../../cands/001.cand', [False, False]],
['../../cands/002.cand', [False, False]],
['../../cands/003.cand', [False, False]],
['../../cands/004.cand', [False, False]],
['../../cands/005.cand', [False, False]],
['../../cands/006.cand', [False, False]],
['../../cands/007.cand', [False, False]],
['../../cands/008.cand', [False, False]],
['../../cands/009.cand', [True, True]],
['../../cands/010.cand', [True, True]],
['../../cands/011.cand', [False, True]],
['../../cands/012.cand', [True, True]],
['../../cands/013.cand', [False, True]],
['../../cands/014.cand', [False, True]],
['../../cands/015.cand', [True, True]],
['../../cands/016.cand', [False, False]],
['../../cands/017.cand', [False, False]],
['../../cands/018.cand', [False, False]],
['../../cands/019.cand', [False, False]],
['../../cands/020.cand', [False, False]],
['../../cands/021.cand', [False, False]],
['../../cands/022.cand', [False, False]],
['../../cands/023.cand', [False, False]],
['../../cands/024.cand', [False, False]],
['../../cands/025.cand', [False, False]],
['../../cands/026.cand', [False, True]],
['../../cands/027.cand', [False, True]],
['../../cands/028.cand', [True, True]],
['../../cands/029.cand', [False, False]],
['../../cands/030.cand', [False, False]],
['../../cands/031.cand', [False, True]],
['../../cands/032.cand', [True, True]],
['../../cands/033.cand', [False, False]],
['../../cands/034.cand', [False, False]],
['../../cands/035.cand', [False, False]],
['../../cands/036.cand', [False, False]],
['../../cands/037.cand', [False, False]],
['../../cands/038.cand', [False, False]],
['../../cands/039.cand', [False, False]],
['../../cands/040.cand', [False, False]],
['../../cands/041.cand', [False, False]],
['../../cands/042.cand', [False, False]],
['../../cands/043.cand', [False, False]],
['../../cands/044.cand', [False, False]],
['../../cands/045.cand', [False, False]],
['../../cands/046.cand', [False, False]],
['../../cands/047.cand', [True, True]],
['../../cands/048.cand', [False, False]],
['../../cands/049.cand', [False, True]],
['../../cands/050.cand', [False, False]],
['../../cands/051.cand', [False, True]],
['../../cands/052.cand', [False, False]],
['../../cands/053.cand', [False, True]],
['../../cands/054.cand', [False, False]],
['../../cands/055.cand', [False, False]],
['../../cands/056.cand', [False, False]],
['../../cands/057.cand', [False, False]],
['../../cands/058.cand', [False, False]],
['../../cands/059.cand', [False, False]],
['../../cands/060.cand', [False, False]],
['../../cands/061.cand', [False, False]],
['../../cands/062.cand', [False, False]],
['../../cands/063.cand', [False, False]],
['../../cands/064.cand', [False, False]],
['../../cands/065.cand', [False, False]],
['../../cands/066.cand', [False, True]],
['../../cands/067.cand', [False, False]],
['../../cands/068.cand', [False, False]],
['../../cands/069.cand', [False, False]],
['../../cands/070.cand', [True, True]],
['../../cands/071.cand', [False, False]],
['../../cands/072.cand', [False, True]],
['../../cands/073.cand', [False, False]],
['../../cands/074.cand', [False, True]],
['../../cands/075.cand', [False, False]],
['../../cands/076.cand', [True, True]],
['../../cands/077.cand', [False, False]],
['../../cands/078.cand', [False, False]],
['../../cands/079.cand', [False, False]],
['../../cands/080.cand', [False, False]],
['../../cands/081.cand', [False, False]],
['../../cands/082.cand', [False, False]],
['../../cands/083.cand', [False, False]],
['../../cands/084.cand', [False, False]],
['../../cands/085.cand', [False, False]],
['../../cands/086.cand', [False, False]],
['../../cands/087.cand', [False, False]],
['../../cands/088.cand', [False, False]],
['../../cands/089.cand', [False, False]],
['../../cands/090.cand', [False, False]],
['../../cands/091.cand', [False, False]],
['../../cands/092.cand', [False, False]],
['../../cands/093.cand', [False, False]],
['../../cands/094.cand', [False, False]],
['../../cands/095.cand', [False, False]],
['../../cands/096.cand', [False, False]],
['../../cands/097.cand', [False, False]],
['../../cands/098.cand', [False, False]],
['../../cands/099.cand', [True, True]],
['../../cands/100.cand', [True, True]],
['../../cands/101.cand', [True, True]],
['../../cands/102.cand', [False, True]],
['../../cands/103.cand', [False, True]],
['../../cands/104.cand', [True, True]],
['../../cands/105.cand', [True, True]],
['../../cands/106.cand', [False, False]],
['../../cands/107.cand', [False, False]],
['../../cands/108.cand', [False, False]],
['../../cands/109.cand', [False, False]],
['../../cands/110.cand', [False, False]],
['../../cands/111.cand', [False, False]],
['../../cands/112.cand', [False, True]],
['../../cands/113.cand', [False, False]],
['../../cands/114.cand', [False, False]],
['../../cands/115.cand', [False, False]],
['../../cands/116.cand', [False, False]],
['../../cands/117.cand', [False, False]],
['../../cands/118.cand', [False, False]],
['../../cands/119.cand', [False, False]],
['../../cands/120.cand', [False, False]],
['../../cands/121.cand', [False, False]],
['../../cands/122.cand', [False, False]],
['../../cands/123.cand', [False, False]],
['../../cands/124.cand', [False, False]],
['../../cands/125.cand', [False, False]],
['../../cands/126.cand', [False, False]],
['../../cands/127.cand', [False, False]],
['../../cands/128.cand', [False, False]],
['../../cands/129.cand', [False, False]],
['../../cands/130.cand', [False, True]],
['../../cands/131.cand', [False, True]],
['../../cands/132.cand', [False, False]],
['../../cands/133.cand', [True, True]],
['../../cands/134.cand', [True, True]],
['../../cands/135.cand', [False, False]],
['../../cands/136.cand', [False, True]],
['../../cands/137.cand', [False, False]],
['../../cands/138.cand', [False, True]],
['../../cands/139.cand', [False, False]],
['../../cands/140.cand', [False, False]],
['../../cands/141.cand', [False, False]],
['../../cands/142.cand', [False, False]],
['../../cands/143.cand', [False, False]],
['../../cands/144.cand', [False, False]],
['../../cands/145.cand', [False, False]],
['../../cands/146.cand', [False, False]],
['../../cands/147.cand', [False, False]],
['../../cands/148.cand', [False, False]],
['../../cands/149.cand', [False, False]],
['../../cands/150.cand', [False, False]],
['../../cands/151.cand', [False, True]],
['../../cands/152.cand', [False, False]],
['../../cands/153.cand', [False, True]],
['../../cands/154.cand', [False, False]],
['../../cands/155.cand', [False, True]],
['../../cands/156.cand', [False, False]],
['../../cands/157.cand', [False, True]],
['../../cands/158.cand', [False, False]],
['../../cands/159.cand', [True, True]],
['../../cands/160.cand', [False, False]],
['../../cands/161.cand', [False, False]],
['../../cands/162.cand', [False, False]],
['../../cands/163.cand', [False, False]],
['../../cands/164.cand', [False, False]],
['../../cands/165.cand', [False, False]],
['../../cands/166.cand', [False, False]],
['../../cands/167.cand', [False, False]],
['../../cands/168.cand', [False, False]],
['../../cands/169.cand', [False, False]],
['../../cands/170.cand', [False, False]],
['../../cands/171.cand', [False, False]],
['../../cands/172.cand', [False, False]],
['../../cands/173.cand', [False, False]],
['../../cands/174.cand', [False, False]],
['../../cands/175.cand', [False, False]],
['../../cands/176.cand', [False, False]],
['../../cands/177.cand', [False, False]],
['../../cands/178.cand', [False, False]],
['../../cands/179.cand', [False, False]],
['../../cands/180.cand', [False, False]],
['../../cands/181.cand', [False, False]],
['../../cands/182.cand', [False, False]],
['../../cands/183.cand', [True, True]],
['../../cands/184.cand', [True, True]],
['../../cands/185.cand', [True, True]],
['../../cands/186.cand', [True, True]],
['../../cands/187.cand', [True, True]],
['../../cands/188.cand', [True, True]],
['../../cands/189.cand', [True, True]],
['../../cands/190.cand', [False, False]],
['../../cands/191.cand', [True, True]],
['../../cands/192.cand', [True, True]],
['../../cands/193.cand', [True, True]],
['../../cands/194.cand', [False, True]],
['../../cands/195.cand', [False, True]],
['../../cands/196.cand', [True, True]],
['../../cands/197.cand', [False, False]],
['../../cands/198.cand', [True, True]],
['../../cands/199.cand', [True, True]],
['../../cands/200.cand', [False, True]],
['../../cands/201.cand', [False, False]],
['../../cands/202.cand', [True, True]],
['../../cands/203.cand', [True, True]],
['../../cands/204.cand', [True, True]],
['../../cands/205.cand', [False, False]],
['../../cands/206.cand', [False, True]],
['../../cands/207.cand', [False, True]],
['../../cands/208.cand', [True, True]],
['../../cands/209.cand', [False, False]],
['../../cands/210.cand', [True, True]],
['../../cands/211.cand', [True, True]],
['../../cands/212.cand', [False, False]],
['../../cands/213.cand', [True, True]],
['../../cands/214.cand', [False, False]],
['../../cands/215.cand', [False, False]],
['../../cands/216.cand', [False, False]],
['../../cands/217.cand', [False, False]],
['../../cands/218.cand', [False, False]],
['../../cands/219.cand', [True, True]],
['../../cands/220.cand', [True, True]],
['../../cands/221.cand', [False, False]],
['../../cands/222.cand', [True, True]],
['../../cands/223.cand', [True, True]],
['../../cands/224.cand', [True, True]],
['../../cands/225.cand', [True, True]],
['../../cands/226.cand', [False, False]],
['../../cands/227.cand', [True, True]],
['../../cands/228.cand', [True, True]],
['../../cands/229.cand', [True, True]],
['../../cands/230.cand', [True, True]],
['../../cands/231.cand', [True, True]],
['../../cands/232.cand', [True, True]],
['../../cands/233.cand', [True, True]],
['../../cands/234.cand', [False, True]],
['../../cands/235.cand', [False, True]],
['../../cands/236.cand', [False, True]],
['../../cands/237.cand', [False, True]],
['../../cands/238.cand', [True, True]],
['../../cands/239.cand', [True, True]],
['../../cands/240.cand', [False, True]],
['../../cands/241.cand', [False, False]],
['../../cands/242.cand', [False, False]],
['../../cands/243.cand', [False, False]],
['../../cands/244.cand', [False, False]],
['../../cands/245.cand', [False, False]],
['../../cands/246.cand', [False, False]],
['../../cands/247.cand', [False, False]],
['../../cands/248.cand', [False, False]],
['../../cands/249.cand', [False, False]],
['../../cands/250.cand', [False, False]],
['../../cands/251.cand', [False, False]],
['../../cands/252.cand', [False, False]],
['../../cands/253.cand', [True, True]],
['../../cands/254.cand', [False, True]],
['../../cands/255.cand', [False, True]],
['../../cands/256.cand', [True, True]],
['../../cands/257.cand', [False, True]],
['../../cands/258.cand', [False, True]],
['../../cands/259.cand', [True, True]],
['../../cands/260.cand', [True, True]],
['../../cands/261.cand', [False, False]],
['../../cands/262.cand', [False, False]],
['../../cands/263.cand', [False, False]],
['../../cands/264.cand', [False, False]],
['../../cands/265.cand', [False, False]],
['../../cands/266.cand', [False, False]],
['../../cands/267.cand', [True, True]],
['../../cands/268.cand', [True, True]],
['../../cands/269.cand', [True, True]]]

def _main():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Calculate False Activation Rate (FAR) for each observation inside ../../cands')
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
			#if row[1] != '-1.0':	# A duration of -1 indicates inability to determine duration
			new_obs = [SB_path + row[0]]
			obs.append(new_obs)
	
	header = '# num, duration (s), minimum rate (s^-1), maximum rate (s^-1)\n'
	intf = '%d '
	floatf = '%0.5f '
	formats = intf + floatf + floatf + floatf + floatf + '\n'

	# Iterate over all observations, calculate FAR for both simple cutoff and grouping by processing
	# in groups of n_cands
	for o in obs:
		#sbid_obsid = o[opath][len(SB_path):len(SB_path)+22]
		print(o[opath])

		#len(SB_path):len(SB_path)+22 makes the outfile names simply SBID/OBSID....csv
		#smp_outfile = smp_path + sbid_obsid + '_s' + str(snmin) +  '.csv'
		#grp_outfile = grp_path + sbid_obsid + '_s' + str(snmin) + '_r' + str(rsqmmin) + '.csv'

		# Total number of candidates post-grouping
		n_smp = 0
		n_grp = 0

		# Time of last candidate (seconds after file start)
		t_last = 0.0001

		skip_this = []
		for sk in skips:
			if sk[0] == o[opath]:
				skip_this = sk[1]

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
						#if not os.path.isfile(smp_outfile):
						smp_cands = simple_filter(cands, snmin)
	
						n_smp += len(smp_cands)
	
						# Grouping
						grp_cands = []
						#if not os.path.isfile(grp_outfile):
						grp_cands = external_grouping(cands, dmmin, wmax, snmin, True, rsqmmin)
	
						n_grp += len(grp_cands)
	
					#Output counter
					count += len(cands)
					ii += 1
					if ii % 1000 == 0:
						print("%d/%d - %0.1f%%" % (count, total_cands, float(count)/float(total_cands)*100))

					#Skip
					#for j in range(3960):
					#	line = f.readline()

					cands = []


		#min_smp_rate = float(n_smp)/o[odur]
		max_smp_rate = float(n_smp)/t_last

		#min_grp_rate = float(n_grp)/o[odur]
		max_grp_rate = float(n_grp)/t_last

		smp_write = formats % (n_smp, 1, t_last, 0, max_smp_rate)
		grp_write = formats % (n_grp, 1, t_last, 0, max_grp_rate)

		print(header)
		print(smp_write)
		print(grp_write)

		if args.minimal:
			min_outfile = '%ss%0.3f_r%0.3f' % (min_path, snmin, rsqmmin)
			with open(min_outfile, 'a') as mf:
				mf.write(o[opath] + '\n')
				mf.write(smp_write)
				mf.write(grp_write)
		else:
			"""
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
			"""

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
