#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys

def calculate_FAR(df):
	if df['t_range'] > 0:
		return df['n_aft']/df['t_range']
	else:
		return np.nan

secs_per_time_sample = 8.652e-4

fname = sys.argv[1]

fdf = pd.read_csv(fname, sep=' ', comment='#')
fdf.columns = ['File', 'n_bef', 'n_aft', 't_first', 't_last']

fdf['t_range'] = fdf['t_last'] - fdf['t_first'] + 1

fdf['FAR'] = fdf.apply(calculate_FAR, axis=1)

fdf['FAR (s^-1)'] = fdf['FAR']/secs_per_time_sample

fdf['FAR (h^-1)'] = fdf['FAR (s^-1)']*3600

#print(fdf.describe())

#print("MEAN FAR (h^-1)")
print("%f" % ((fdf['n_aft'].mean() / fdf['t_range'].mean())*3600/secs_per_time_sample))
