import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
dt = 0.5
# Preprocessing spike trains and LFPs;
spike_file = 'raster.csv'
lfp_file = 'I.bin'
spike_ind = sys.argv[1]
lfp_ind = sys.argv[2]
if lfp_ind == 'all':
  lfp_list = ''
  for i in range(100):
    lfp_list += str(i)
    if i != 99:
      lfp_list += ','
      lfp_ind = lfp_list
t_range = [500, 10000000]
str_range = str(t_range[0]) + ',' + str(t_range[1])
subprocess.call(['./bin/spike.out', './data/dataRepo/triple2/' + spike_file, './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', './data/dataRepo/triple2/' + spike_file, './data/tmp/singleSpike_shuffle.bin', spike_ind, str_range, str(dt), 'true'])
subprocess.call(['./bin/lfp.out', './data/dataRepo/triple2/' + lfp_file, './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
