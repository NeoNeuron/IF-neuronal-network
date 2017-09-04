#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import mylib
target_ind = 100
# x = np.genfromtxt('./data/x.csv', delimiter = ',', skip_header = target_ind, max_rows = 1)
# np.savetxt('./data/xn.csv', x, fmt = '%.10f')
# subprocess.call(['./bin/surrogate2d.out', './data/y.csv', './data/y_shuffle.csv', '1'])
y = np.genfromtxt('./data/y.csv', delimiter = ',')
y_shuffle = y.copy()
np.random.shuffle(y_shuffle)
np.savetxt('./data/y_shuffle.csv', y_shuffle, fmt = '%.10f', delimiter = ',')
# y_shuffle = np.genfromtxt('./data/y_shuffle.csv', delimiter = ',')
lrange = 20
rrange = 20

subprocess.call(['./bin/mi_dd.out', './data/x.csv', './data/y.csv', str(target_ind), str(lrange) + ',' + str(rrange)])
mis = pd.read_csv('./data/mi/mi_dd.csv')
subprocess.call(['./bin/mi_dd.out', './data/x.csv', './data/y_shuffle.csv', str(target_ind), str(lrange) + ',' + str(rrange)])
mis_shuffle = pd.read_csv('./data/mi/mi_dd.csv')

plt.plot(mis['timelag'], mis['mi'], label = "mi")
plt.plot(mis_shuffle['timelag'], mis_shuffle['mi'], label = "surrogate data")
plt.xlabel('Strength')
plt.ylabel('Mutual Information')
plt.legend()
# plt.savefig('./linear.png')s
plt.show()
