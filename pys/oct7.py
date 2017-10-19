#!/usr/bin/python
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
s=[0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014]
# prepare path of data;
tmp_dir = './data/tmp/'
saving_dir = './data/dataRepo/'
s_dir = [saving_dir + str(element) + '/' for element in s]

# locate the position of target setting;
nrow=11
ncol=30
mis = np.zeros(len(s))
counter = 0
for si in s_dir:
    # prepare new directory of new data;
    subprocess.call(['mv', si, tmp_dir + 'I.csv'])
    subprocess.call(['mv', si, tmp_dir + 'raster.csv'])
    # prepare spike train and lfp series;
    subprocess.call(['./pys/breakTS.py', '0', '1'])
    # excute mutual info calculation program;
    subprocess.call(['./bin/mi_bd.out', tmp_dir + 'singleSpike.csv', tmp_dir + 'singleI.csv', '20', '0,1', '500'])
    # import data of mi_bd.csv
    mi = pd.read_csv('./data/mi/mi_bd.csv')
    mis[counter] = mi['mi'][1]

plt.plot(s, mis)
plt.xlabel('interaction intensity')
plt.ylabel('mutual info')
plt.show()
