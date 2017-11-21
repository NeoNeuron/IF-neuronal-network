#!/usr/bin/python
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
s = np.arange(0.005, 0.015, 0.0005)
# prepare path of data;
tmp_dir = './data/tmp/'
saving_dir = './data/dataRepo/'
s_dir = [saving_dir + str(element) + '/' for element in s]

# locate the position of target setting;
indx = 20
mi_max = np.zeros(len(s))
mi_max_ind = np.zeros(len(s))
baselines = np.zeros(len(s))
snr = np.zeros(len(s))
counter = 0
dt = 0.5
spike_ind = 0
lfp_ind = 1
bin_num = 50
for si in s_dir:
    # copy data to working directory;
    subprocess.call(['cp', si + 'I.csv', tmp_dir])
    subprocess.call(['cp', si + 'raster.csv', tmp_dir])
    # prepare spike train and lfp series;
    subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike.csv', str(spike_ind), '500,600000', str(dt), 'false'])
    subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike_shuffle.csv', str(spike_ind), '500,600000', str(dt), 'true'])
    subprocess.call(['./bin/lfp.out', './data/tmp/I.csv', './data/tmp/singleI.csv', str(lfp_ind), '500,600000', str(dt)])
    # excute mutual info calculation program, and import data;
    subprocess.call(['./bin/mi_bd.out', tmp_dir + 'singleSpike.csv', tmp_dir + 'singleI.csv', '5,5', str(dt), str(bin_num)])
    mi = pd.read_csv('./data/mi/mi_bd.csv')
    subprocess.call(['./bin/mi_bd.out', tmp_dir + 'singleSpike_shuffle.csv', tmp_dir + 'singleI.csv', '5,5', str(dt), str(bin_num)])
    mi_shuffle = pd.read_csv('./data/mi/mi_bd.csv')
    mi_max[counter] = mi['mi'].max()
    baselines[counter] = mi_shuffle['mi'].mean()
    snr[counter] = mi_max[counter]/baselines[counter]
    mi_max_ind[counter] = mi['mi'].argmax() - 5
    counter += 1
np.savetxt('baselines.csv', baselines, delimiter = ',', fmt = '%.10f')
np.savetxt('mi_max.csv', mi_max, delimiter = ',', fmt = '%.10f')
np.savetxt('snr.csv', snr, delimiter = ',', fmt = '%.10f')
print mi_max_ind
plt.figure(figsize = (18, 6), dpi = 60)
plt.subplot(1,3,1)
plt.plot(s, mi_max)
plt.xlabel('Strength')
plt.ylabel('Maximum mutual information')
ax = plt.gca()
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.subplot(1,3,2)
plt.plot(s, baselines)
plt.xlabel('Strength')
plt.ylabel('Baseline (mutual info)')
ax = plt.gca()
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.subplot(1,3,3)
plt.plot(s, snr)
plt.xlabel('strength')
plt.ylabel('SNR')
plt.savefig('oct7.png')
