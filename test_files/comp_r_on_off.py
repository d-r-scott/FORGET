#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

r_on = pd.read_csv('r_on.out', sep=' ')
r_off = pd.read_csv('r_off.out', sep=' ')

r_on.columns = ['Reduced', 'number', 'of', 'candidates', 'in', 'file', 'from', 'before', 'to', 'after']
r_off.columns = ['Reduced', 'number', 'of', 'candidates', 'in', 'file', 'from', 'before', 'to', 'after']

r_on = r_on[['before', 'after']]
r_off = r_off[['before', 'after']]

r_on['frac'] = r_on['after']/r_on['before']
summary_on = r_on.describe()
print('R^2 calculation on:')
print(summary_on)
print('\n')

r_off['frac'] = r_off['after']/r_off['before']
summary_off = r_off.describe()
print('R^2 calculation off:')
print(summary_off)
print('\n')

bins = np.linspace(0,1,40)

fig = plt.figure()
gs = gridspec.GridSpec(3,1)

ax = fig.add_subplot(gs[0,0])
plt.title('Fraction of candidates post-grouping')
r_on['frac'].plot.hist(alpha=0.5, label='R^2 calculation on', bins=bins)
r_off['frac'].plot.hist(alpha=0.5, label='R^2 calculation off', bins=bins)
ax.set_xlim([0,1])
plt.legend()

ax = fig.add_subplot(gs[1, 0])
r_on[['frac']].boxplot(vert=False)
ax.set_yticks([])
ax.set_ylabel('R^2 filtering on')
ax.set_xlim([0,1])

ax = fig.add_subplot(gs[2, 0])
r_off[['frac']].boxplot(vert=False)
ax.set_yticks([])
ax.set_ylabel('R^2 filtering off')
ax.set_xlim([0,1])

plt.show()
