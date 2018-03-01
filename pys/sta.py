#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
dt = 0.5
preflag = int(sys.argv[3])
# Preprocessing spike trains and LFPs;
if preflag == 0:
    spike_file = 'raster.csv'
    lfp_file = 'I.bin'
    spike_ind = sys.argv[1]
    lfp_ind = sys.argv[2]
    if lfp_ind == 'all':
        lfp_list = ''
        for i in range(100):
            lfp_list += str(i)
            if i != 99:
                lfp_list += ','
                lfp_ind = lfp_list
    t_range = [500, 600000]
    str_range = str(t_range[0]) + ',' + str(t_range[1])
    subprocess.call(['./bin/spike.out', './data/tmp/' + spike_file, './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
    subprocess.call(['./bin/lfp.out', './data/tmp/' + lfp_file, './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
# Calculating Spike Triggerd Average:
ntd = 30
ptd = 30
drange = str(ntd) + ',' + str(ptd)
subprocess.call(['./bin/sta.out', './data/tmp/singleSpike.bin', './data/tmp/singleI.bin', drange])
data = pd.read_csv('data/sta/sta.csv')
# Plotting:
plotway = 'linear'
if plotway == 'linear':
    plt.plot(data['timelag']*dt, data['sta'], label = 'mi')
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((0,1))
elif plotway == 'log':
    plt.semilogy(data['timelag']*dt, data['sta'], label = 'mi')

plt.xlabel('Delay Time(ms)')
plt.ylabel('Spike Triggered Average')
plt.legend()
plt.savefig('sta')
plt.close()
