#!/bin/python
# this script aims to draw rasterogram of neural network
import numpy as np
import matplotlib.pyplot as plt
import sys

filepath = sys.argv[1]
f = open(filepath)
plt.figure()
counter = 1
for line in f:
    spike_str = line.replace('\n', '').split(',')
    spike_str.remove('')
    spikes = [float(i) for i in spike_str]
    counter_list = np.ones(len(spikes)) * counter
    plt.scatter(spikes, counter_list, s = 1, c = 'b')
    counter += 1
plt.xlabel('time(ms)')
plt.ylabel('indices')
f.close();
plt.show()
