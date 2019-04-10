import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description = "Calculate spike-spike mutual information and plot results.")
parser.add_argument('dir', type = str, help = 'directory of source data and output data')
parser.add_argument('xid', type = str, help = 'index of pre-neuron.')
parser.add_argument('yid', type = str, help = 'index of post-neuron.')
parser.add_argument('delay_range', metavar = 'drange', type = str, help = 'range of time delay, seperated by comma')

args = parser.parse_args()

t_range = '5e2,1e5'
dt = 1
# define tmp files
tmpname_spike1 = '/singleSpike1.csv'
tmpname_spike2 = '/singleSpike2.csv'
tmpname_spike2_shuffle = '/singleSpike2_shuffle.csv'
tmpname_mi = '/mi.csv'
tmpname_mi_shuffle = '/mi_shuffle.csv'

subprocess.call(['./bin/spike.out', args.dir + '/raster.csv', args.dir + tmpname_spike1, args.xid, t_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', args.dir + '/raster.csv', args.dir + tmpname_spike2, args.yid, t_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', args.dir + '/raster.csv', args.dir + tmpname_spike2_shuffle, args.yid, t_range, str(dt), 'true'])
subprocess.call(['./bin/mi_bb.out', args.dir + tmpname_spike1, args.dir + tmpname_spike2, args.dir + tmpname_mi, args.drange])
subprocess.call(['./bin/mi_bb.out', args.dir + tmpname_spike1, args.dir + tmpname_spike2_shuffle,  args.dir + tmpname_mi_shuffle, args.drange])

mi = np.genfromtxt(args.dir + tmpname_mi, delimiter = ',')
mi_shuffle = np.genfromtxt(args.dir + tmpname_mi_shuffle, delimiter = ',')

fig = plt.figure(figsize = (10,6), dpi = 60)
ax = fig.subplots(1,1)
ax.plot(mi[:,0]*dt, mi[:,1], label = 'mi')
ax.plot(mi_shuffle[:,0]*dt, mi_shuffle[:,1], label = 'mi_shuffle')
ax.set_xlabel('Time Lags (ms)')
ax.set_ylabel('Mutual Info')
ax.legend()
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.savefig('mi_bd.eps')
plt.close()
# clean tmp files
os.remove(args.dir + tmpname_spike1)
os.remove(args.dir + tmpname_spike2)
os.remove(args.dir + tmpname_spike2_shuffle)
os.remove(args.dir + tmpname_mi)
os.remove(args.dir + tmpname_mi_shuffle)
