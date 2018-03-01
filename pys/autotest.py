#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import os, sys

path = sys.argv[1]
tmax = int(sys.argv[2])
opath = sys.argv[3].split('.')
opath = './data/tmp/'+ opath[0] + 'T.csv'
subprocess.call(['./bin/transpose.out', path, opath])
#subprocess.call(['./bin/means.out', opath])
#subprocess.call(['./bin/stds.out', opath])
#subprocess.call(['./bin/vector1d.out', opath, './data/tmp/vector1d.csv', '0', '1', '1000,20000'])
auto = np.zeros((800, tmax + 1));
for i in range(800):
    subprocess.call(['./bin/vector1d.out', opath, './data/tmp/vector1d.csv', str(i), '1', '1000,20000'])
# subprocess.call(['./bin/transpose.out', './data/tmp/vector1d.csv', opath])
    subprocess.call(['./bin/ac.out', './data/tmp/vector1d.csv', str(tmax)])
# subprocess.call(['./bin/ac.out', opath, str(tmax)])

    new_auto = np.loadtxt('./data/stationary/ac.csv')
    auto [i] = new_auto;
#means = np.loadtxt('./data/stationary/means.csv')
#stds = np.loadtxt('./data/stationary/stds.csv')
#means = means/means.max()
#stds = stds/stds.max()
np.savetxt('./data/stationary/acs.csv', auto, fmt = '%10f', delimiter = ',')
auto = auto.sum(0)
auto = auto/auto.max()
plt.plot(np.arange(tmax + 1)/32, auto, label = "Autocovarience")
#plt.plot(np.arange(tmax + 1)/32, means[0:tmax+1], label = "Mean")
#plt.plot(np.arange(tmax + 1)/32, stds[0:tmax+1], label = "Std")
plt.legend()
#plt.xlabel('Time delay/Time(ms)')
plt.xlabel('Time delay(ms')
#plt.ylabel('Covariance/Mean/Std(normalized)')
plt.ylabel('Covariance(normalized)')
plt.savefig('./data/figure/' + sys.argv[3])
os.system('rm ' + opath)
