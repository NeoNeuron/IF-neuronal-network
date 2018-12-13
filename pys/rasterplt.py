#!/bin/python
# this script aims to draw rasterogram of neural network
import numpy as np
import matplotlib.pyplot as plt
import sys
import configparser as cp

typepath = sys.argv[1]
types = np.genfromtxt(typepath, delimiter = ',', dtype = 'b')
filepath = sys.argv[2]
f = open(filepath)
#============================
# import config file
#============================
config = cp.ConfigParser()
config.read('doc/config_net.ini')
neuron_num = int(config.get('Network Parameters', 'NeuronNumber'))
tmax = float(config.get('Time', 'MaximumTime'))
dt = 5
t = np.arange(0, tmax, dt)
counts_exc = np.zeros(len(t))
counts_inh = np.zeros(len(t))
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (8,6), dpi = 60, sharex = True)
counter = 1
for line in f:
    spike_str = line.replace('\n', '').split(',')
    spike_str.remove('')
    spikes = [float(i) for i in spike_str if float(i) < tmax]
    counter_list = np.ones(len(spikes)) * counter
    if types[counter - 1]:
        ax1.scatter(spikes, counter_list, s = 1, c = 'b')
        for ele in spikes:
            counts_exc[int(ele/dt)] += 1
    else:
        ax1.scatter(spikes, counter_list, s = 1, c = 'r')
        for ele in spikes:
            counts_inh[int(ele/dt)] += 1
    counter += 1
f.close();
ax1.set_xlim(0, tmax)
ax1.set_ylim(0, counter + 1)
ax1.set_ylabel('Indices')
ax1.grid(linestyle='--')
ax2.plot(t + dt/2, counts_exc, linestyle = '--', linewidth = 1, label='exc neuron')
ax2.plot(t + dt/2, counts_inh, linestyle = '--', linewidth = 1, label='inh neuron')
ax2.plot(t + dt/2, counts_exc + counts_inh, linewidth = 1, label='tot neuron')
ax2.set_ylabel('Mean Firing Counts')
ax2.legend()
ax2.set_xlabel('Time(ms)')
fig.savefig('raster.eps')
plt.close()
