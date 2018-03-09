import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
# Preprocessing spike trains and LFPs;
dt = 0.5
spike_file = 'raster.csv'
lfp_file = 'I.bin'
spike_ind = '0'
lfp_ind = '1'
ntd = 30
ptd = 30
drange = str(ntd) + ',' + str(ptd)
hist_bin = 10
algs_mode = 'direct'
tmax = np.logspace(3, 7, num = 5)
bls = np.zeros(len(tmax))
for i in range(len(tmax)):
  t_range = [500, 500 + int(tmax[i])]
  str_range = str(t_range[0]) + ',' + str(t_range[1])
  subprocess.call(['./bin/spike.out', './data/tmp/' + spike_file, './data/tmp/singleSpike_shuffle.bin', spike_ind, str_range, str(dt), 'true'])
  subprocess.call(['./bin/lfp.out', './data/tmp/' + lfp_file, './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
  # Calculating Mutual information:
  subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike_shuffle.bin', './data/tmp/singleI.bin', drange, str(dt), str(hist_bin), algs_mode])
  data_shuffle = pd.read_csv('data/mi/mi_bd.csv')
  bls[i] = data_shuffle['mi'].mean()

np.savetxt('bls_tmax_ifnet.txt', bls, fmt = '%.18e')
# Plotting:
plt.loglog(tmax, bls)
plt.xlabel('Length of Time Series(ms)')
plt.ylabel('Baseline of MI')
plt.grid(True, which = 'both', axis = 'both')
plt.savefig('mi_baseline.png')
plt.close()
