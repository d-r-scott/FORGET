#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

r_on_fn = "r_on_FAR_rsq.out"

r_on = []

with open(r_on_fn, 'r') as f:
	for i, line in enumerate(f):
		r_on.append(float(line))

rsq_min = 0
rsq_max = 1
rsq_step = 0.01

rsq = np.arange(rsq_min, rsq_max+rsq_step, rsq_step)

plt.plot(rsq, r_on, label="Grouping with R^2 filtering")

plt.legend()
plt.show()
