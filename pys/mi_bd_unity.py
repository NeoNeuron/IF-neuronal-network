import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
dt = 0.5
# Calculating Mutual information:
ntd = 15
ptd = 15
drange = str(ntd) + ',' + str(ptd)
path = sys.argv[1]
spike_ind = sys.argv[2]
lfp_ind = sys.argv[3]
binsize = sys.argv[4]
trange = '500,10000000'
program = './bin/mi_bd_unity.out'
subprocess.call([program, path + 'raster.csv', path + 'I.bin', spike_ind, lfp_ind, trange, str(dt), drange, binsize])
data = pd.read_csv('data/mi/mi_bd_unity.csv')
# Plotting:
plotway = 'linear'
if plotway == 'linear':
  plt.plot(data['timelag']*dt, data['mi'], label = 'mi')
  plt.plot(data['timelag']*dt, data['mi_shuffle'], label = 'mi_shuffle')
  ax = plt.gca()
  ax.yaxis.get_major_formatter().set_powerlimits((0,1))
elif plotway == 'log':
  plt.semilogy(data['timelag']*dt, data['mi'], label = 'mi')
  plt.semilogy(data['timelag']*dt, data['mi_shuffle'], label = 'mi_shuffle')

plt.xlabel('Delay Time(ms)')
plt.ylabel('mutual info')
plt.legend()
plt.savefig('mi')
plt.close()
