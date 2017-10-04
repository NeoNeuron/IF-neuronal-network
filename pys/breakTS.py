#!/usr/bin/python
import numpy as np
import subprocess
import sys

dt = 0.125 # time step of analysis
auto_time_scale = 40 # unit millisecond
spike_ind = sys.argv[1]
lfp_ind = sys.argv[2]

# load current and spike train
subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike.csv', spike_ind, '500,20000', str(dt)])
spike = np.loadtxt('./data/tmp/singleSpike.csv')
subprocess.call(['./bin/lfp.out', './data/tmp/I.csv', './data/tmp/singleI.csv', lfp_ind, '500,20000', str(dt)])
lfp = np.loadtxt('./data/tmp/singleI.csv')

length = int(auto_time_scale / dt) # length of array within the auto_time_scale;
N = int(len(spike) / length)
spike = spike[0:length*N]
spike = np.reshape(spike, (N,length)).transpose()
np.savetxt('./data/tmp/rasterRST.csv', spike, delimiter = ',')
# subprocess.call(['./bin/transpose.out', './data/tmp/rasterRS.csv', './data/tmp/rasterRST.csv'])

lfp = lfp[0:length*N]
lfp = np.reshape(lfp, (N,length)).transpose()
np.savetxt('./data/tmp/IRST.csv', lfp, fmt = '%.10f', delimiter = ',')
# subprocess.call(['./bin/transpose.out', './data/tmp/IRS.csv', './data/tmp/IRST.csv'])
