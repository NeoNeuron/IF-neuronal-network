#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import subprocess
datafile = sys.argv[1]
subprocess.call(['./bin/transpose.out', datafile, './data/tmp/dataT.csv'])
subprocess.call(['./bin/means.out', './data/tmp/dataT.csv'])
subprocess.call(['./bin/stds.out', './data/tmp/dataT.csv'])
subprocess.call(['./bin/autocov.out', './data/tmp/dataT.csv'])
# import data
means = np.loadtxt('./data/stationary/means.csv', delimiter = ',')
stds = np.loadtxt('./data/stationary/stds.csv', delimiter = ',')
autocov = np.loadtxt('./data/stationary/autocov.csv', delimiter = ',')

dt = 0.03125
plt.plot(np.arange(len(means))/dt,means,label = 'Mean')
plt.plot(np.arange(len(stds))/dt, stds, label = 'Std')
plt.plot(np.arange(len(autocov))/dt, autocov, label = 'AutoCov')
plt.legend()
figname = datafile.split('.')[0] + '.png'
plt.savefig(figname)
