#!/bin/python
# this script aims to draw rasterogram of neural network
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys

filepath = sys.argv[1]
f = open(filepath)
tmax = 10000
plt.figure(figsize = (12,6), dpi = 80)
counter = 1
for line in f:
    spike_str = line.replace('\n', '').split(',')
    spike_str.remove('')
    spikes = [float(i) for i in spike_str if float(i) < tmax]
    counter_list = np.ones(len(spikes)) * counter
    plt.scatter(spikes, counter_list, s = 1)
    counter += 1
f.close();
plt.xlim(0, tmax)
plt.ylim(0, counter + 1)
plt.xlabel('Time(ms)')
plt.ylabel('Indices')
plt.savefig('raster')
plt.close()
