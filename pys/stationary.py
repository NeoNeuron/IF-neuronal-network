#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import subprocess
datafile = sys.argv[1]
subprocess.call(['./bin/transpose.out', datafile, './data/tmp/dataT.csv'])
#subprocess.call(['./bin/means.out', './data/tmp/dataT.csv'])
#subprocess.call(['./bin/stds.out', './data/tmp/dataT.csv'])
subprocess.call(['./bin/autocov.out', './data/tmp/dataT.csv', '4000', '3200', '100'])
# import data
means = np.loadtxt('./data/stationary/means.csv', delimiter = ',')
stds = np.loadtxt('./data/stationary/stds.csv', delimiter = ',')
autocov = np.loadtxt('./data/stationary/autocov.csv', delimiter = ',')

means = means / means.max()
stds = stds / stds.max()
autocov = autocov / autocov.max()

means = means[1:4000]
stds = stds[1:4000]

dt = 0.03125
plt.plot(np.arange(len(means))*dt,means,label = 'Mean')
plt.plot(np.arange(len(stds))*dt, stds, label = 'Std')
plt.plot(np.arange(len(autocov))*dt, autocov, label = 'AutoCov')
plt.legend()
plt.xlabel('Time /Delay Time(ms)')
plt.ylabel('Covariance / Mean / Std (Normalized)')
figname = datafile.split('.')[-2]
figname = figname.split('/')[-1] + '.png'
plt.savefig(figname)
