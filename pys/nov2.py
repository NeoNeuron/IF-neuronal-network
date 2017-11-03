#!/usr/bin/python
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
# prepare path of data;
tmp_dir = './data/tmp/'
saving_dir = sys.argv[1]

# locate the position of target setting;
totaln = 40
indx = int(sys.argv[2])
ntd = 15
ptd = 15
if (indx - ntd < 0): ntd = indx
if (indx + ptd > 39): ptd = 39 - indx
# copy data to working directory;
subprocess.call(['cp', saving_dir + 'I.csv', tmp_dir])
subprocess.call(['cp', saving_dir + 'raster.csv', tmp_dir])
# prepare spike train and lfp series;
dt = 0.5
spike_ind = 0
lfp_ind = 1
subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike.csv', str(spike_ind), '500,600000', str(dt), 'false'])
subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike_shuffle.csv', str(spike_ind), '500,600000', str(dt), 'true'])
subprocess.call(['./bin/lfp.out', './data/tmp/I.csv', './data/tmp/singleI.csv', str(lfp_ind), '500,600000', str(dt)])
# excute mutual info calculation program;
subprocess.call(['./bin/mi_bd.out', tmp_dir + 'singleSpike.csv', tmp_dir + 'singleI.csv', str(indx), str(ntd) + ',' + str(ptd), '500'])
# import data of mi_bd.csv
mi = pd.read_csv('./data/mi/mi_bd.csv')

# excute mutual info calculation for base level;
subprocess.call(['./bin/mi_bd.out', tmp_dir + 'singleSpike_shuffle.csv', tmp_dir + 'singleI.csv', str(indx), str(ntd) + ',' + str(ptd), '500'])
# import data of mi_bd.csv
mi_shuffle = pd.read_csv('./data/mi/mi_bd.csv')
# mis[counter] = mi['mi'].max()
# mis_ind = mi['mi'].argmax() - 3

plt.plot(mi['timelag']*dt, mi['mi'], label = 'mi')
plt.plot(mi_shuffle['timelag']*dt, mi_shuffle['mi'], label = 'mi_shuffle')
plt.xlabel('Time Delay (ms)')
plt.ylabel('mutual info')
figname = saving_dir.split('/')
figname.remove('')
figname = 'mi_' + figname[-1] + '.png'
plt.savefig(figname)
