# -*- coding: utf-8 -*-
#!/usr/bin/python3
import numpy as np
import subprocess
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

parser = argparse.ArgumentParser(description = "Script to draw the spike triggled average in LFP case, execute after pre_mi_bd.py.")
parser.add_argument('dt', type = float, nargs = 1, help = 'Timing step of time series')
parser.add_argument('trange', metavar = 'range_of_time_window', type = str, nargs = 1, help = 'range of time window for STA')
args = parser.parse_args()
subprocess.call(['./bin/sta.out', './data/tmp/singleSpike.bin', './data/tmp/singleI.bin', args.trange[0]])

dt = args.dt[0]
# import data
dat = np.genfromtxt('./data/sta/sta.csv', delimiter = ',')
# plot data
plt.figure()
plt.plot(dt * dat[:,0], dat[:,1])
ax = plt.gca()
## plot settings:
#xmajorLocator   = MultipleLocator(2) 
#xmajorFormatter = FormatStrFormatter('%d') # set formats of x axis
#xminorLocator   = MultipleLocator(1) 
#
#ymajorLocator   = MultipleLocator(0.1e-2) 
#ymajorFormatter = FormatStrFormatter('%1.1f') # set formats of y axis
#yminorLocator   = MultipleLocator(0.02e-2) 
#
ax = plt.gca()
#ax.xaxis.set_major_formatter(xmajorFormatter)
#ax.xaxis.set_major_locator(xmajorLocator)
#ax.xaxis.set_minor_locator(xminorLocator)
#
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
#ax.yaxis.set_major_locator(ymajorLocator)
#ax.yaxis.set_minor_locator(yminorLocator)

plt.xlabel('Time Delay (ms)')
plt.ylabel('Spike-triggerd Average LFP')
plt.grid()
plt.savefig('sta')
#plt.show()
plt.close()
