import numpy as np
import matplotlib.pyplot as plt
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description = "Prepare the data files of spike train and local field potential 'LFP' for TDMI calculation.")
parser.add_argument('dir', type = str, help = 'directory of neuronal data')
parser.add_argument('spike_id', type = str, help = 'index of spiking neuron')
parser.add_argument('trange', metavar = 'time range', type = str, help = 'range of time series')
parser.add_argument('bin_size', type = str, help = 'bin size of joint probability distribution')
parser.add_argument('td', metavar = 'time delay', type = str, help = 'negative and positive maximum delay time, seperated by comma')
args = parser.parse_args()

dt = 0.5
# Preprocessing spike trains and LFPs;
spike_file = args.dir + '/raster.csv'
lfp_file = args.dir + '/lfp.bin'
# Generating tmp filename
tmpname_spike = args.dir + '/singleSpike.bin'
tmpname_spike_shuffle = args.dir + '/singleSpikeShuffle.bin'

subprocess.call(['./bin/spike.out', spike_file, tmpname_spike, args.spike_id, args.trange, str(dt), 'false'])
subprocess.call(['./bin/spike.out', spike_file, tmpname_spike_shuffle, args.spike_id, args.trange, str(dt), 'true'])

# Calculating Mutual information:
# Select target program
program = './bin/mi_bd.out'
# Generating tmp mifile name;
mi_name1 = args.dir + '/mi_bd.csv'
print('>> Output mifile --> ' + mi_name1)
subprocess.call([program, tmpname_spike, lfp_file, mi_name1, args.td, args.bin_size])
data = np.genfromtxt(mi_name1, delimiter = ',')
mi_name2 = args.dir + '/mi_bd_shuffle.csv'
print('>> Output mifile --> ' + mi_name2)
subprocess.call([program, tmpname_spike_shuffle, lfp_file, mi_name2, args.td, args.bin_size])
data_shuffle = np.genfromtxt(mi_name2, delimiter = ',')
# Plotting:
plotway = 'linear'
if plotway == 'linear':
  plt.plot(data[:,0]*dt, data[:,1], label = 'mi')
  plt.plot(data_shuffle[:,0]*dt, data_shuffle[:,1], label = 'mi_shuffle')
  ax = plt.gca()
  ax.yaxis.get_major_formatter().set_powerlimits((0,1))
elif plotway == 'log':
  plt.semilogy(data[:,0]*dt, data[:,1], label = 'mi')
  plt.semilogy(data_shuffle[:,0]*dt, data_shuffle[:,1], label = 'mi_shuffle')

plt.xlabel('Delay Time(ms)')
plt.ylabel('mutual info')
plt.grid(linestyle = '--')
plt.legend()
plt.savefig(args.dir + '/mi' + args.spike_id + '.png')
plt.close()
os.remove(tmpname_spike)
os.remove(tmpname_spike_shuffle)
