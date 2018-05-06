import numpy as np
import subprocess
import sys
dt = 0.5
# Preprocessing spike trains and LFPs;
spike_file = 'raster.csv'
lfp_file = 'I.bin'
path = sys.argv[1]
spike_ind = sys.argv[2]
lfp_ind = sys.argv[3]
if lfp_ind == 'all':
  lfp_list = ''
  for i in range(100):
    lfp_list += str(i)
    if i != 99:
      lfp_list += ','
      lfp_ind = lfp_list
t_range = [500, 10000000]
str_range = str(t_range[0]) + ',' + str(t_range[1])
subprocess.call(['./bin/spike.out', path + spike_file, './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', path + spike_file, './data/tmp/singleSpike_shuffle.bin', spike_ind, str_range, str(dt), 'true'])
subprocess.call(['./bin/lfp.out', path + lfp_file, './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
