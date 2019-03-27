import numpy as np
import matplotlib.pyplot as plt
import subprocess
import argparse

parser = argparse.ArgumentParser(description = "Calculate TDMI of spike train and local field potential between two neurons.")
parser.add_argument('dir', type = str, help = 'directory of neuronal data')
parser.add_argument('spike_id', type = str, help = 'index of spiking neuron')
parser.add_argument('lfp_id', type = str, help = 'index of lfp neuron')
parser.add_argument('trange', metavar = 'time range', type = str, help = 'range of time series')
parser.add_argument('bin_size', type = str, help = 'bin size of joint probability distribution')
parser.add_argument('td', metavar = 'time delay', type = str, help = 'negative and positive maximum delay time, seperated by comma')
args = parser.parse_args()

dt = 0.5
# Preprocessing spike trains and LFPs;
spike_file = 'raster.csv'
lfp_file = 'I.bin'
# Generating tmp filename
tmpname_spike = './data/spike/singleSpike.bin'
tmpname_spike_shuffle = './data/spike/singleSpikeShuffle.bin'
tmpname_lfp = './data/lfp/singleI.bin'

subprocess.call(['./bin/spike.out', args.dir + spike_file, tmpname_spike, args.spike_id, args.trange, str(dt), 'false'])
subprocess.call(['./bin/spike.out', args.dir + spike_file, tmpname_spike_shuffle, args.spike_id, args.trange, str(dt), 'true'])
subprocess.call(['./bin/lfp.out', args.dir + lfp_file, tmpname_lfp])

# Calculating Mutual information:
# Select target program
program = './bin/mi_bd.out'
# Generating tmp mifile name;
mi_name1 = './data/mi/mi_bd.csv'
print('>> Output mifile --> ' + mi_name1)
subprocess.call([program, tmpname_spike, tmpname_lfp, mi_name1, args.td, args.bin_size])
data = np.genfromtxt(mi_name1, delimiter = ',')
mi_name2 = './data/mi/mi_bd_shuffle.csv'
print('>> Output mifile --> ' + mi_name2)
subprocess.call([program, tmpname_spike_shuffle, tmpname_lfp, mi_name2, args.td, args.bin_size])
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
plt.savefig('mi.eps')
plt.close()
