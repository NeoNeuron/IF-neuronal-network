#!/usr/bin/python
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
s = np.arange(0.005, 0.015, 0.0005)
# prepare path of data;
tmp_dir = './data/tmp/'
saving_dir = './data/dataRepo/'
s_dir = [saving_dir + str(element) + '/' for element in s]

# locate the position of target setting;
indx = 20
mis = np.zeros(len(s))
mis_ind = np.zeros(len(s))
counter = 0
for si in s_dir:
    # copy data to working directory;
    subprocess.call(['cp', si + 'I.csv', tmp_dir])
    subprocess.call(['cp', si + 'raster.csv', tmp_dir])
    # prepare spike train and lfp series;
    dt = 0.5
    spike_ind = 0
    lfp_ind = 1
    subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike.csv', str(spike_ind), '500,600000', str(dt), 'false'])
    subprocess.call(['./bin/lfp.out', './data/tmp/I.csv', './data/tmp/singleI.csv', str(lfp_ind), '500,600000', str(dt)])
    # excute mutual info calculation program;
    subprocess.call(['./bin/mi_bd.out', tmp_dir + 'singleSpike.csv', tmp_dir + 'singleI.csv', str(indx), '3,5', '500'])
    # import data of mi_bd.csv
    mi = pd.read_csv('./data/mi/mi_bd.csv')
    mis[counter] = mi['mi'].max()
    mis_ind = mi['mi'].argmax() - 3
    counter += 1

print mis_ind
plt.plot(s, mis)
plt.xlabel('interaction intensity')
plt.ylabel('mutual info')
plt.show()
