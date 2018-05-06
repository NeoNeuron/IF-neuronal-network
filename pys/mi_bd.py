import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
import myio

dt = 0.5
# Preprocessing spike trains and LFPs;
path = sys.argv[1]
spike_ind = sys.argv[2]
lfp_ind = sys.argv[3]
spike_file = 'raster.csv'
lfp_file = 'I.bin'
if lfp_ind == 'all':
  lfp_list = ''
  for i in range(100):
    lfp_list += str(i)
    if i != 99:
      lfp_list += ','
      lfp_ind = lfp_list
t_range = [500, 10000000]
str_range = str(t_range[0]) + ',' + str(t_range[1])
# Generating tmp filename
name_num = uni_num()
tmpname_spike = './data/spike/singleSpike' + name_num + '.bin'
tmpname_spike_shuffle = './data/spike/singleSpikeShuffle' + name_num + '.bin'
tmpname_lfp = './data/lfp/singleI' + name_num + '.bin'

subprocess.call(['./bin/spike.out', path + spike_file, tmpname_spike, spike_ind, str_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', path + spike_file, tmpname_spike_shuffle, spike_ind, str_range, str(dt), 'true'])
subprocess.call(['./bin/lfp.out', path + lfp_file, tmpname_lfp, lfp_ind, str_range, str(dt)])

binsize = sys.argv[4]
drange = sys.argv[5]
# Calculating Mutual information:
# Select target program
program = './bin/mi_bd.out'
# Generating tmp mifile name;
mi_name1 = './data/mi/mi_bd' + uni_num() + '.csv'
print '>> Output mifile --> ' + mi_name1
subprocess.call([program, tmpname_spike, tmpname_lfp, drange, binsize])
data = pd.read_csv(mi_name1)
mi_name2 = './data/mi/mi_bd' + uni_num() + '.csv'
print '>> Output mifile --> ' + mi_name2
subprocess.call([program, tmpname_spike_shuffle, tmpname_lfp, drange, binsize])
data_shuffle = pd.read_csv(mi_name2)
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
# Clear tmp files;
if sys.argv[6] == 'True':
    subprocess(['rm -vf', tmpname_spike, tmpname_spike_shuffle, tmpname_lfp])
