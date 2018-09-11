#!/usr/bin/python3
import numpy as np
import subprocess
import argparse

parser = argparse.ArgumentParser(description = "Prepare the data files of spike train and local field potential 'LFP' for TDMI calculation.")
#parser.add_argument('dt', type = float, nargs = 1, default = 0.5, help = 'timing step of discrete time series of spike train and local field potential')
parser.add_argument('dir', type = str, nargs = 1, help = 'directory of raw data')
parser.add_argument('spike_id', type = str, nargs = 1, help = 'index of spike train')
parser.add_argument('lfp_id', type = str, nargs = 1, help = 'index of local field potential')
parser.add_argument('max_t_range', metavar = 'maximum_time_range', type = str, nargs = 1, help = 'maximum of time series')
args = parser.parse_args()
#print args
dt = 0.5
# Preprocessing spike trains and LFPs;
spike_file = 'raster.csv'
lfp_file = 'I.bin'
DIR = args.dir[0]
spike_ind = args.spike_id[0]
lfp_ind = args.lfp_id[0]
if lfp_ind == 'all':
  lfp_list = ''
  for i in range(100):
    lfp_list += str(i)
    if i != 99:
      lfp_list += ','
      lfp_ind = lfp_list
str_range = '500,' + str(args.max_t_range[0])
subprocess.call(['./bin/spike.out', DIR + spike_file, './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', DIR + spike_file, './data/tmp/singleSpike_shuffle.bin', spike_ind, str_range, str(dt), 'true'])
# generate total current of targeted neuron
subprocess.call(['./bin/lfp.out', DIR + lfp_file, './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])

# filter upper bound
freq_up = 300
f = open(DIR + lfp_file, 'rb')
shape = st.unpack('QQ', f.read(16))
lfp = np.array(st.unpack('d'*shape[1], f.read(8*shape[1])))
lfp_fft = np.fft.fft(lfp)
t = np.arange(0,len(lfp))*dt
# set frequency above upper bound to zero
freq = np.fft.fftfreq(len(lfp))*(1000/dt)
lfp_fft[abs(freq)>freq_up] = 0
lfp_new = lfp_fft.ifft(lfp_fft)
# output new data:

