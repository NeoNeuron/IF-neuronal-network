#!/usr/bin/python
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
s = np.arange(0.001, 0.015, 0.0005)
# prepare path of data;
tmp_dir = './data/tmp/'
saving_dir = './data/dataRepo/'
s_dir = [saving_dir + str(element) + '/' for element in s]
# locate the position of target setting;
mi_max = np.zeros(len(s))
mi_max_ind = np.zeros(len(s))
snr = np.zeros(len(s))
counter = 0
dt = 0.5
spike_ind = 0
lfp_ind = 1
mi_mode = 'direct'
bin_num = 20
for si in s_dir:
    # copy data to working directory;
    subprocess.call(['cp', si + 'I.csv', tmp_dir])
    subprocess.call(['cp', si + 'raster.csv', tmp_dir])
    # prepare spike train and lfp series;
    subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike.csv', str(spike_ind), '500,60000', str(dt), 'false'])
    subprocess.call(['./bin/lfp.out', './data/tmp/I.csv', './data/tmp/singleI.csv', str(lfp_ind), '500,60000', str(dt)])
    # excute mutual info calculation program, and import data;
    subprocess.call(['./bin/mi_bd.out', tmp_dir + 'singleSpike.csv', tmp_dir + 'singleI.csv', '5,5', str(dt), str(bin_num), mi_mode])
    mi = pd.read_csv('./data/mi/mi_bd.csv')
    mi_max[counter] = mi['mi'].max()
    mi_max_ind[counter] = mi['mi'].argmax() - 5
    counter += 1
np.savetxt('mi_max.csv', mi_max, delimiter = ',', fmt = '%.10f')
print mi_max_ind
plt.figure(figsize = (10, 8), dpi = 60)
plt.plot(s, mi_max)
plt.xlabel('Strength')
plt.ylabel('Maximum mutual information')
ax = plt.gca()
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.savefig('oct7.png')
