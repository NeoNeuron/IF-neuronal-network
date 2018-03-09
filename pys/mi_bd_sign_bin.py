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
t_range = [500, 600000]
str_range = str(t_range[0]) + ',' + str(t_range[1])
# subprocess.call(['./bin/spike.out', './data/tmp/' + spike_file, './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', './data/tmp/' + spike_file, './data/tmp/singleSpike_shuffle.bin', spike_ind, str_range, str(dt), 'true'])
subprocess.call(['./bin/lfp.out', './data/tmp/' + lfp_file, './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
# Calculating Mutual information:
ntd = 30
ptd = 30
drange = str(ntd) + ',' + str(ptd)
hist_bin = np.arange(10,81,2)
algs_mode = 'direct'
bls = np.zeros(len(hist_bin))
for i in range(len(hist_bin)):
  subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike_shuffle.bin', './data/tmp/singleI.bin', drange, str(dt), str(hist_bin[i]), algs_mode])
  data_shuffle = pd.read_csv('data/mi/mi_bd.csv')
  bls[i] = data_shuffle['mi'].mean()

np.savetxt('bls_bins_ifnet.png', bls, fmt = '%.18e')
# Plotting:
plt.semilogy(1.0/hist_bin, bls)
plt.xlabel('Number of bins')
plt.ylabel('Baseline of MI')
plt.grid(True, which = 'both', axis = 'y')
plt.savefig('mi')
plt.close()
