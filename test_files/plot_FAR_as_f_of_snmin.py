#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

simple_fn = "simple_filter_FAR.out"
r_off_fn = "r_off_FAR.out"
r_on_fn = "r_on_FAR.out"

simple = []
r_off = []
r_on = []

with open(simple_fn, 'r') as f:
	for i, line in enumerate(f):
		simple.append(float(line))

with open(r_off_fn, 'r') as f:
	for i, line in enumerate(f):
		r_off.append(float(line))

with open(r_on_fn, 'r') as f:
	for i, line in enumerate(f):
		r_on.append(float(line))

sn_min = 7
sn_max = 12
sn_step = 0.1

sn = np.arange(sn_min, sn_max+sn_step, sn_step)

plt.semilogy(sn, simple, label="Simple filtering")
plt.semilogy(sn, r_off, label="Grouping")
plt.semilogy(sn, r_on, label="Grouping with R^2 filtering")

plt.legend()
plt.show()
