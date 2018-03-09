import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
hist_bin = np.arange(10,81,3)
algs_mode = 'direct'
tmax = np.logspace(3, 7, num = 5)
bls = np.zeros((len(hist_bin), len(tmax)))
for i in range(len(tmax)):
  t_range = [500, 500 + int(tmax[i])]
  str_range = str(t_range[0]) + ',' + str(t_range[1])
  subprocess.call(['./bin/spike.out', './data/tmp/' + spike_file, './data/tmp/singleSpike_shuffle.bin', spike_ind, str_range, str(dt), 'true'])
  subprocess.call(['./bin/lfp.out', './data/tmp/' + lfp_file, './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
  # Calculating Mutual information:
  for j in range(len(hist_bin)):
    subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike_shuffle.bin', './data/tmp/singleI.bin', drange, str(dt), str(hist_bin[j]), algs_mode])
    data_shuffle = pd.read_csv('data/mi/mi_bd.csv')
    bls[j][i] = data_shuffle['mi'].mean()

np.savetxt('bls_map_ifnet.txt', bls, fmt = '%.18e')
# Plotting:
fig = plt.figure()
pcm = plt.imshow(bls, aspect = 'auto', cmap = cm.tab20c, norm = colors.LogNorm(vmin = bls.min(), vmax = bls.max()), extent=[2.5, 7.5, 79.5, 9.5])
plt.xlabel('Length of Time Series(10^x ms)')
plt.ylabel('Number of bins')
fig.colorbar(pcm, extend = 'max')
ax=plt.gca()
ax.set_xticks(np.linspace(3,7,5))
ax.set_yticks(np.linspace(79,10,bls.shape[0]))
plt.savefig('mi_baseline_map.png')
plt.close()
