#!/usr/bin/python
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
s = np.arange(0.001, 0.015, 0.0005)

# prepare path of data;
tmp_saving = './data/tmp/'
saving_dir = './data/dataRepo/'
# locate the position of target setting;
nrow=11
ncol=30
mis = np.zeros(len(s))
counter = 0
for si in s:
    fini = open('./doc/config_net.ini', 'r+')
    for i in range(nrow - 1):
        line = fini.readline()
    fini.seek(ncol - 1, 1)
    pos_begin = fini.tell()
    fini.write('%.5g'%si)
    pos_end = fini.tell()
    if pos_end < pos_begin + 6:
        fini.write(' '*(pos_begin + 6 - pos_end))
    fini.close()
    # excute neuronal simulation program;
    subprocess.call(['./bin/net.out', tmp_saving])
    # prepare spike train and lfp series;
    # subprocess.call(['./pys/breakTS.py', '0', '1'])
    # excute mutual info calculation program;
    # subprocess.call(['./bin/mi_bd.out', tmp_saving + 'singleSpike.csv', tmp_saving + 'singleI.csv', '20', '0,1', '500'])
    # import data of mi_bd.csv
    # mi = pd.read_csv('./data/mi/mi_bd.csv')
    # mis[counter] = mi['mi'][1]
    # prepare new directory of new data;
    subprocess.call(['mkdir', saving_dir + str(si)])
    subprocess.call(['mv', tmp_saving + 'I.bin', saving_dir + str(si) + '/'])
    subprocess.call(['mv', tmp_saving + 'raster.csv', saving_dir + str(si) + '/'])
    subprocess.call(['cp', './doc/config_net.ini', saving_dir + str(si) + '/'])
# print s
# print mis
