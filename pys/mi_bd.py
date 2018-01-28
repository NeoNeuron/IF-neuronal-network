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
ntd = 35
ptd = 35
hist_bin = 20
# prepare data container:
t_range = [0, 0]
t_range[0] = 500
t_range[1] = 60000
algs_mode = sys.argv[3]
drange = str(ntd) + ',' + str(ptd)
str_range = str(t_range[0]) + ',' + str(t_range[1])
subprocess.call(['./bin/spike.out', './data/tmp/rasterPre.csv', './data/tmp/singleSpike.csv', spike_ind, str_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', './data/tmp/rasterPre.csv', './data/tmp/singleSpike_shuffle.csv', spike_ind, str_range, str(dt), 'true'])

subprocess.call(['./bin/lfp.out', './data/tmp/postI.csv', './data/tmp/singleI.csv', lfp_ind, str_range, str(dt)])
subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike.csv', './data/tmp/singleI.csv', drange, str(dt), str(hist_bin), algs_mode])
data = pd.read_csv('data/mi/mi_bd.csv')
plt.plot(data['timelag']*dt, data['mi'], label = 'mi')
subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike_shuffle.csv', './data/tmp/singleI.csv', drange, str(dt), str(hist_bin), algs_mode])
data_shuffle = pd.read_csv('data/mi/mi_bd.csv')
plt.plot(data_shuffle['timelag']*dt, data_shuffle['mi'], label = 'mi_shuffle')
plt.xlabel('Delay Time(ms)')
plt.ylabel('mutual info')
plt.legend()
ax = plt.gca()
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.savefig('mi')
