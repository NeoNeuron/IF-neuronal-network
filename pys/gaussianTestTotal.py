#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
num_rand = 60
num_trial = 300000
x = np.zeros((num_rand, num_trial))
y = np.zeros((num_rand, num_trial))
x[0,:] = np.random.normal(size = num_trial)
y[0,:] = np.random.normal(size = num_trial)
ax = 0.1
ay = 0.1
axy = np.arange(0,10,0.5)
target_ind = 30
f = open('./data/mi/mi.csv', 'w')
f.close()
for k in axy:
# k = 0.1
    for i in range(num_rand - 1):
        x[i + 1,:] = ax*x[i,:] + np.random.normal(size = num_trial)
        y[i + 1,:] = ay*y[i,:] + k*x[i,:] + np.random.normal(size = num_trial)
    # np.savetxt('./data/x.csv', x, fmt = '%.10f')
    np.savetxt('./data/xn.csv', x[target_ind,:], fmt = '%.10f')
    # y_shuffle = y.copy()
    # np.random.shuffle(y_shuffle)
    # np.savetxt('./data/y.csv', y, fmt = '%.10f')
    np.savetxt('./data/yn.csv', y[target_ind + 1,:], fmt = '%.10f')
    # subprocess.call(['./bin/mi_dd.out', './data/x.csv', './data/y.csv', str(target_ind), '10,10'])
    subprocess.call(['./bin/mi.out', './data/xn.csv', './data/yn.csv'])
    # mi = pd.read_csv('./data/mi/mi_dd.csv')
        # f = open('mis.csv', 'w')
        # f.close()
        # lrange = 10
        # rrange = 10
        # for i in range(lrange):
        #     np.savetxt('./data/yn.csv', y[target_ind-lrange-i,:], fmt = '%.10f')
        #     subprocess.call(['./bin/mi.out', './data/xn.csv', './data/yn.csv'])
        # for i in range(rrange + 1):
        #     np.savetxt('./data/yn.csv', y[target_ind+i,:], fmt = '%.10f')
        #     subprocess.call(['./bin/mi.out', './data/xn.csv', './data/yn.csv'])
        #
        # mis = np.loadtxt('mis.csv')
        #
        # f = open('mis.csv', 'w')
        # f.close()
        # for i in range(lrange):
        #     np.savetxt('./data/yn.csv', y_shuffle[target_ind-lrange-i,:], fmt = '%.10f')
        #     subprocess.call(['./bin/mi.out', './data/xn.csv', './data/yn.csv'])
        # for i in range(rrange + 1):
        #     np.savetxt('./data/yn.csv', y_shuffle[target_ind+i,:], fmt = '%.10f')
        #     subprocess.call(['./bin/mi.out', './data/xn.csv', './data/yn.csv'])

mis = np.loadtxt('./data/mi/mi.csv')
plt.plot(axy, mis, label = 'Test')
plt.plot(axy, 0.5*np.log(1+axy**2), label = 'Theory')
plt.xlabel('Strength')
plt.ylabel('Mutual Information')
plt.legend()
# plt.savefig('./linear.png')
plt.show()
