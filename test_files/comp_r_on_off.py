#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

r_on = pd.read_csv('r_on.out', sep=' ')
r_off = pd.read_csv('r_off.out', sep=' ')

del r_on['to']
del r_off['to']

r_on['frac'] = r_on['after']/r_on['before']
summary_on = r_on.describe()
print(summary_on)

r_off['frac'] = r_off['after']/r_off['before']
summary_off = r_off.describe()
print(summary_off)


plt.figure()
r_on['frac'].plot.hist(alpha=0.5)
r_off['frac'].plot.hist(alpha=0.5)
plt.show()
