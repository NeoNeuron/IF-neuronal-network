import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
dt = 0.5
# Calculating Mutual information:
ntd = 2
ptd = 2
drange = str(ntd) + ',' + str(ptd)
# Calculation parameters; binsize for mi_bd.out, threshold for mi_bd_2bins.out
program = './bin/mi_bd.out'
parameter = sys.argv[1]
subprocess.call([program, './data/tmp/singleSpike.bin', './data/tmp/singleI.bin', drange, parameter])
data = pd.read_csv('data/mi/mi_bd.csv')
subprocess.call([program, './data/tmp/singleSpike_shuffle.bin', './data/tmp/singleI.bin', drange, parameter])
data_shuffle = pd.read_csv('data/mi/mi_bd.csv')
# Plotting:
plotway = 'linear'
if plotway == 'linear':
  plt.plot(data['timelag']*dt, data['mi'], label = 'mi')
  plt.plot(data_shuffle['timelag']*dt, data_shuffle['mi'], label = 'mi_shuffle')
  ax = plt.gca()
  ax.yaxis.get_major_formatter().set_powerlimits((0,1))
elif plotway == 'log':
  plt.semilogy(data['timelag']*dt, data['mi'], label = 'mi')
  plt.semilogy(data_shuffle['timelag']*dt, data_shuffle['mi'], label = 'mi_shuffle')

plt.xlabel('Delay Time(ms)')
plt.ylabel('mutual info')
plt.legend()
plt.savefig('mi')
plt.close()
