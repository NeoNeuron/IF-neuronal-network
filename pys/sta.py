#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import subprocess
import sys
# Calculating Spike Triggerd Average:
drange = sys.argv[1]
#subprocess.call(['./bin/sta.out', './data/tmp/singleSpike.bin', './data/tmp/singleI.bin', drange])
#data = pd.read_csv('data/sta/sta.csv')
data = pd.read_csv('sta.csv')
# Plotting:
dt = 0.5
plotway = 'linear'
if plotway == 'linear':
    plt.plot(data['timelag']*dt, data['sta'], label = 'mi')
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((0,1))
elif plotway == 'log':
    plt.semilogy(data['timelag']*dt, data['sta'], label = 'mi')
# plot settings:
xmajorLocator   = MultipleLocator(25) 
xmajorFormatter = FormatStrFormatter('%d') # set formats of x axis
xminorLocator   = MultipleLocator(5) 

ymajorLocator   = MultipleLocator(0.01e-2) 
#ymajorFormatter = FormatStrFormatter('%1.1f') # set formats of y axis
yminorLocator   = MultipleLocator(0.002e-2) 

plt.xlabel('Time Delay (ms)')
plt.ylabel('Mean LFP')
ax = plt.gca()
ax.xaxis.set_major_formatter(xmajorFormatter)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)

ax.yaxis.get_major_formatter().set_powerlimits((0,1))
ax.yaxis.set_major_locator(ymajorLocator)
ax.yaxis.set_minor_locator(yminorLocator)

plt.grid()
plt.legend()
plt.savefig('sta')
plt.close()
