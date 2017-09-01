#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import mylib
x = np.zeros(10000)
y = np.zeros(10000)
x[0] = np.random.normal()
y[0] = np.random.normal()
ax = -0.1
ay = -0.1
axy = 0.5
for i in range(9999):
    x[i + 1] = ax*x[i] + np.random.normal()
    # if i < 3:
    y[i + 1] = ay*y[i] + axy*x[i]*x[i] + np.random.normal()
    # else:
    #     y[i + 1] = ay*y[i] + axy*x[i]**2 + axy*np.exp(x[i - 3]/2) + np.random.normal()

np.savetxt('./data/x.csv', x, fmt = '%.10f')
np.savetxt('./data/y.csv', y, fmt = '%.10f')

ntd = 100
ptd = 100
# Calculate mutual information and STA
subprocess.call(['./bin/cc.out', './data/x.csv', './data/y.csv', str(ntd) + ',' + str(ptd)])
subprocess.call(['./bin/surrogate.out', './data/y.csv', './data/y_shuffle.csv'])
subprocess.call(['./bin/mi_dd.out', './data/x.csv', './data/y.csv', str(ntd) + ',' + str(ptd)])
mi = pd.read_csv('./data/mi/mi_dd.csv')
subprocess.call(['./bin/mi_dd.out', './data/x.csv', './data/y_shuffle.csv', str(ntd) + ',' + str(ptd)])
mi_shuffle = pd.read_csv('./data/mi/mi_dd.csv')
mi['shuffle'] = mi_shuffle['mi']
cc = pd.read_csv('./data/cc/cc.csv')
plt.subplot(1, 2, 1)
plt.plot(mi['timelag'], mi['mi'])
plt.plot(mi['timelag'], mi['shuffle'])
plt.subplot(1, 2, 2)
plt.plot(cc['timelag'], cc['cc'])
plt.show()
