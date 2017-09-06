#!/usr/bin/python
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import sys

path = sys.argv[1]
tmax = int(sys.argv[2])
opath = './data/tmp/dataT.csv'
subprocess.call(['./transpose.out', path, opath])
subprocess.call(['./bin/means.out', opath])
subprocess.call(['./bin/stds.out', opath])
subprocess.call(['./bin/vector1d.out', opath, './data/tmp/vector1d.csv', '0', '1'])
subprocess.call(['./bin/ac.out', './data/tmp/vector1d.csv', str(tmax)])

auto = np.loadtxt('./data/stationary/ac.csv')
means = np.loadtxt('./data/stationary/means.csv')
stds = np.loadtxt('./data/stationary/stds.csv')
means = means/means.max()
stds = stds/stds.max()
auto = auto/auto.max()
plt.plot(np.arange(tmax + 1)/32, auto, label = "Autocovarience")
plt.plot(np.arange(tmax + 1)/32, means[0:tmax+1], label = "Mean")
plt.plot(np.arange(tmax + 1)/32, stds[0:tmax+1], label = "Std")
plt.legend()
plt.xlabel('Time delay/Time(ms)')
plt.ylabel('Covariance/Mean/Std(normalized)')
# plt.show()
