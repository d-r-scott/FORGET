#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# Indexes for file_arr fields
sbid  = 0	# Schedblock number
obsid = 1	# Unique observation id
nums  = 2	# List of numbers of candidates (0: simple, 1: grouped)
durs  = 3	# List of durations of observations
minrs = 4	# List of minimum FARs
maxrs = 5	# List of maximum FARs

# Indexes for the sublists
smp = 0
grp = 1

def _main():
	parser = ArgumentParser(description='Parse and plot given (minimal) FAR .out files')
	parser.add_argument(dest='files', nargs='+')
	args = parser.parse_args()

	file_arr_list = []

	for fname in args.files:
		file_arr = []

		# Files are structured like so:
		#	SBID/OBSID
		#	nums[smp] durs[smp] minrs[smp] maxrs[smp]
		#	nums[grp] durs[grp] minrs[grp] maxrs[grp]
		#	...
		with open(fname, 'r') as f:
			i = 0
			count = 0
			for line in f:
				if i == 0:
					file_arr.append([])
					# First line: SBID/OBSID
					new_sbid = int(line[2:7])	# Skip the leading 'SB' and just get the number
					new_obsid = int(line[9:])
					file_arr[count].append(new_sbid)
					file_arr[count].append(new_obsid)
					i += 1
				elif i == 1:
					# Second line: smp fields
					line = line.split(' ')
					new_nums = []
					new_durs = []
					new_minrs = []
					new_maxrs = []
					new_nums.append(int(line[0]))
					new_durs.append(float(line[1]))
					new_minrs.append(float(line[2]))
					new_maxrs.append(float(line[3]))
					file_arr[count].append(new_nums)
					file_arr[count].append(new_durs)
					file_arr[count].append(new_minrs)
					file_arr[count].append(new_maxrs)
					i += 1
				elif i == 2:
					# Third line
					line = line.split(' ')
					file_arr[count][nums].append(int(line[0]))
					file_arr[count][durs].append(float(line[1]))
					file_arr[count][minrs].append(float(line[2]))
					file_arr[count][maxrs].append(float(line[3]))
					print(file_arr[count])
					i = 0
					count += 1

		file_arr_list.append(file_arr)

	for file_arr in file_arr_list:
		total_num_smp = sum(row[nums][smp] for row in file_arr)
		total_num_grp = sum(row[nums][grp] for row in file_arr)

		total_dur_smp = sum(row[durs][smp] for row in file_arr)
		total_dur_grp = sum(row[durs][grp] for row in file_arr)

		total_minr_smp = total_num_smp/total_dur_smp
		total_minr_grp = total_num_grp/total_dur_grp

		print(total_minr_smp, total_minr_grp)



if __name__ == '__main__':
	_main()
