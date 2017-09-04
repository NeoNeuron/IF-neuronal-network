#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import sys

pathI = sys.argv[1]
subprocess.call(['./bin/means.out', pathI])
subprocess.call(['./bin/stds.out', pathI])
subprocess.call(['./bin/vector1d.out', pathI, './data/tmp/singleI.csv', '0', '1'])
subprocess.call(['./bin/ac.out', './data/tmp/singleI.csv', '20000'])

auto = np.loadtxt('./data/stationary/ac.csv')
means = np.loadtxt('./data/stationary/means.csv')
stds = np.loadtxt('./data/stationary/stds.csv')
means = means/means.max()
stds = stds/stds.max()
auto = auto/auto.max()
plt.plot(np.arange(20001)/32, auto, label = "Autocovarience")
plt.plot(np.arange(20001)/32, means[0:20001], label = "Mean")
plt.plot(np.arange(20001)/32, stds[0:20001], label = "Std")
plt.legend()
plt.xlabel('Time delay/Time(ms)')
plt.ylabel('Covariance/Mean/Std(normalized)')
plt.show()
