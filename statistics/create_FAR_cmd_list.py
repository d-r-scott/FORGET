#!/usr/bin/env python

import numpy as np

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Create list of commands to calculate False Activation Rate for ranges of sigma and R^2 thresholds', formatter_class=ArgumentDefaultsHelpFormatter)

parser.add_argument('-s', '--smin', type=float, help='Start of sigma threshold range', required=True)
parser.add_argument('-S', '--smax', type=float, help='End of sigma threshold range', required=True)
parser.add_argument('--snum', type=int, help='Number of sigma thresholds in range', required=True)

parser.add_argument('-r', '--rmin', type=float, help='Start of R^2 threshold range', required=True)
parser.add_argument('-R', '--rmax', type=float, help='End of R^2 threshold range', required=True)
parser.add_argument('--rnum', type=int, help='Number of R^2 thresholds in range', required=True)

args = parser.parse_args()

sigmas = np.linspace(args.smin, args.smax, num=args.snum)
rsqmmins = np.linspace(args.rmin, args.rmax, num=args.rnum)

print('sigmas=np.linspace( %f %f %d )' % (args.smin, args.smax, args.snum))
print('rsqmmins=np.linspace( %f %f %d )' % (args.rmin, args.rmax, args.rnum))

with open('FAR_cmd_list.txt', 'w') as f:
	for s in sigmas:
		for r in rsqmmins:
			write_str = './calculate_FAR.py ' + str(s) + ' ' + str(r) + ' > out_FARs/simple_' + str(s) + '_' + str(r) + '.out\n'
			f.write(write_str)
