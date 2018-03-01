#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import subprocess
import sys
import mylib

df = float(sys.argv[1]);
duf = float(sys.argv[2]);
frange = sys.argv[3];
frange = [float(i) for i in frange.split(',')]
ufrange = sys.argv[4];
ufrange = [float(i) for i in ufrange.split(',')]
f = frange[0]
uf = ufrange[0]

# prepare pathes;
loading_dir = "./data/tmp/"
saving_dir = './data/list/'

# setting preliminary parameters
time_lb = 1000
time_ub = 20000
time_range = str(time_lb) + ',' + str(time_ub);
negative_time_delay = 100
positive_time_delay = 200
timing_step = 0.25

# initialize output data file:
df_mi_snr = pd.DataFrame({})
df_mi_pos = pd.DataFrame({})
df_mi_ci = pd.DataFrame({})


ll_mi_snr = []
ll_mi_pos = []
ll_mi_ci = []

while uf < ufrange[1]:
    while f < frange[1]:
        u = uf / f;
        # Edit configurations;
        # fini = open('./doc/config_nets.ini', 'r+')
        # fini.seek(182)
        # fini.write(' '*6)
        # fini.seek(182)
        # fini.write('%.5g'%u)
        # fini.seek(277)
        # fini.write(' '*6)
        # fini.seek(277)
        # fini.write('%.5g'%f)
        # fini.close()
        # configuration for
        fini = open('./doc/config_net.ini', 'r+')
        fini.seek(428)
        fini.write(' '*6)
        fini.seek(428)
        fini.write('%.5g'%u)
        fini.seek(508)
        fini.write(' '*6)
        fini.seek(508)
        fini.write('%.5g'%f)
        fini.close()
        # Execute simulation
        subprocess.call(["./bin/net.out", loading_dir])
        # Generate spike train
        subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', '0', time_range, str(timing_step), './data/spike/spike.csv'])
        # Cauculate LFP from single neuron
        subprocess.call(['./bin/lfp.out', './data/tmp/I.csv', '2', time_range, str(timing_step), './data/lfp/lfp.csv'])
        # Calculate mutual information and STA
        subprocess.call(['./bin/mi_bd.out', './data/spike/spike.csv', './data/lfp/lfp.csv', str(negative_time_delay) + ',' + str(positive_time_delay)])
        subprocess.call(['./bin/sta.out', './data/spike/spike.csv', './data/lfp/lfp.csv', str(negative_time_delay) + ',' + str(positive_time_delay)])
		# run analysis
        mi_data = pd.read_csv('./data/mi/mi_bd.csv')
        mi_snr = mylib.snr(mi_data)
        mi_pos = mylib.peakpos(mi_data)
        mi_ci = mylib.evalpeak(mi_data)
        ll_mi_snr.append(mi_snr)
        ll_mi_pos.append(mi_pos)
        ll_mi_ci.append(mi_ci)
        mylib.plot(('%5.g-'%uf) + ('%5.g'%f) + '.png')
        # # Cauculate LFP from multi neuron
        # subprocess.call['./bin/lfp.out', './tmp/postI.txt', '0,1,2,3,4,5,6,7,8,9', str(time_lb) + ',' + str(time_ub), 'lfp.csv']
        # # Calculate mutual information and STA
        # subprocess.call['./bin/mi.out', '1', str(timing_step, str(negative_time_delay) + ',' + str(positive_time_delay)]
        # subprocess.call['./bin/sta.out', './data/raster/raster.csv', './data/lfp/lfp.csv', str(-timing_step * negative_time_delay) + ',' + str(timing_step * positive_time_delay)]
        f += df
    df_mi_snr[str(uf)] = ll_mi_snr
    df_mi_pos[str(uf)] = ll_mi_pos
    df_mi_ci[str(uf)] = ll_mi_ci

    f = frange[0]
    uf += duf
    ll_mi_snr = []
    ll_mi_pos = []
    ll_mi_ci = []

# print data
df_mi_snr.to_csv(saving_dir + "mi_snr.csv", float_format = '%.4f', index  = False)
df_mi_pos.to_csv(saving_dir + "mi_pos.csv", float_format = '%.4f', index  = False)
df_mi_ci.to_csv(saving_dir + "mi_ci.csv", float_format = '%.4f', index  = False)
