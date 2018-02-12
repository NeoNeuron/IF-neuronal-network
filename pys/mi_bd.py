import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
dt = 0.5
preflag = int(sys.argv[5])
# Preprocessing spike trains and LFPs;
if preflag == 0:
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
    t_range = [500, 3000000]
    str_range = str(t_range[0]) + ',' + str(t_range[1])
    subprocess.call(['./bin/spike.out', './data/tmp/' + spike_file, './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
    subprocess.call(['./bin/spike.out', './data/tmp/' + spike_file, './data/tmp/singleSpike_shuffle.bin', spike_ind, str_range, str(dt), 'true'])
    subprocess.call(['./bin/lfp.out', './data/tmp/' + lfp_file, './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
# Calculating Mutual information:
ntd = 30
ptd = 30
drange = str(ntd) + ',' + str(ptd)
hist_bin = sys.argv[3]
algs_mode = sys.argv[4]
subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike.bin', './data/tmp/singleI.bin', drange, str(dt), hist_bin, algs_mode])
data = pd.read_csv('data/mi/mi_bd.csv')
subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike_shuffle.bin', './data/tmp/singleI.bin', drange, str(dt), hist_bin, algs_mode])
data_shuffle = pd.read_csv('data/mi/mi_bd.csv')
# Plotting:
plotway = 'linear'
if plotway == 'linear':
    plt.plot(data['timelag']*dt, data['mi'], label = 'mi')
    plt.plot(data_shuffle['timelag']*dt, data_shuffle['mi'], label = 'mi_shuffle')
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((0,1))
elif plotway == 'log':
    plt.semilogy(data['timelag']*dt, data['mi'], label = 'mi')
    plt.semilogy(data_shuffle['timelag']*dt, data_shuffle['mi'], label = 'mi_shuffle')

plt.xlabel('Delay Time(ms)')
plt.ylabel('mutual info')
plt.legend()
plt.savefig('mi')
plt.close()
