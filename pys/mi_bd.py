import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
spike_ind = sys.argv[1]
lfp_ind = sys.argv[2]
dt = 0.5
ntd = 15
ptd = 15
hist_bin = sys.argv[3]
# prepare data container:
t_range = [500, 100000]
algs_mode = sys.argv[4]
drange = str(ntd) + ',' + str(ptd)
str_range = str(t_range[0]) + ',' + str(t_range[1])
subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike_shuffle.bin', spike_ind, str_range, str(dt), 'true'])

subprocess.call(['./bin/lfp.out', './data/tmp/I.bin', './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike.bin', './data/tmp/singleI.bin', drange, str(dt), hist_bin, algs_mode])
data = pd.read_csv('data/mi/mi_bd.csv')
plt.plot(data['timelag']*dt, data['mi'], label = 'mi')
subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike_shuffle.bin', './data/tmp/singleI.bin', drange, str(dt), hist_bin, algs_mode])
data_shuffle = pd.read_csv('data/mi/mi_bd.csv')
plt.plot(data_shuffle['timelag']*dt, data_shuffle['mi'], label = 'mi_shuffle')
plt.xlabel('Delay Time(ms)')
plt.ylabel('mutual info')
plt.legend()
ax = plt.gca()
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.savefig('mi')
plt.close()
