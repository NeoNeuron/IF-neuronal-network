#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import mylib
x = np.zeros(100, 1000)
y = np.zeros(100, 1000)
x[0] = np.random.normal()
y[0] = np.random.normal()
ax = -0.1
ay = -0.1
axy = np.arange(0,1,0.05)
ntd = 100
ptd = 100
mis = []
for k in axy:
    for i in range(9999):
        x[i + 1] = np.random.normal()
        y[i + 1] = k*x[i]**2 + np.random.normal()
    np.savetxt('./data/x.csv', x, fmt = '%.10f')
    np.savetxt('./data/y.csv', y, fmt = '%.10f')
    subprocess.call(['./bin/surrogate.out', './data/y.csv', './data/y_shuffle.csv'])
    subprocess.call(['./bin/mi_dd.out', './data/x.csv', './data/y.csv', str(ntd) + ',' + str(ptd)])
    mi = pd.read_csv('./data/mi/mi_dd.csv')
    subprocess.call(['./bin/mi_dd.out', './data/x.csv', './data/y_shuffle.csv', str(ntd) + ',' + str(ptd)])
    mi_shuffle = pd.read_csv('./data/mi/mi_dd.csv')
    mi['shuffle'] = mi_shuffle['mi']
    mi_val = mylib.snr(mi)
    mis.append(mi_val)
    # plt.plot(mi['timelag'], mi['mi'])
    # plt.plot(mi['timelag'], mi['shuffle'])
    # plt.show()
plt.plot(axy, mis)
plt.savefig('square.png')
plt.show()
