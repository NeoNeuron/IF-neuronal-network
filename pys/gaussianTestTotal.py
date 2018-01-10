#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import sys
num_rand = 52
num_trial = 300000
x = np.zeros((num_rand, num_trial))
y = np.zeros((num_rand, num_trial))
x[0,:] = np.random.normal(size = num_trial)
y[0,:] = np.random.normal(size = num_trial)
a = float(sys.argv[1])
b = float(sys.argv[2])
axy = np.arange(0,10,0.5)
target_ind = 50
bins = 150
f = open('./data/mi/mi.csv', 'w')
f.close()
for k in axy:
# k = 0.1
    for i in range(num_rand - 1):
        x[i + 1,:] = a*x[i,:] + np.random.normal(size = num_trial)
        y[i + 1,:] = b*y[i,:] + k*x[i,:] + np.random.normal(size = num_trial)
    np.savetxt('./data/xn.csv', x[target_ind,:], fmt = '%.10f')
    np.savetxt('./data/yn.csv', y[target_ind + 1,:], fmt = '%.10f')
    subprocess.call(['./bin/mi.out', './data/xn.csv', './data/yn.csv', str(bins), str(bins)])

mis = np.loadtxt('./data/mi/mi.csv')
plt.plot(axy, mis, label = 'Test')
if a == b:
    rho2 = axy**2/((1 - a**2)**2 + (1 + a**2)*axy**2)
else:
    rho2 = axy**2*(1 - b**2)/((1 - a**2)*(1 - a*b)**2 + axy**2*(1 - (a*b)**2))

plt.plot(axy, -0.5*np.log(1-rho2), label = 'Theory')
# plt.plot(axy, -0.5*np.log((1 - a**2)/(1 - a**2 + (1 - b**2)*axy**2)), label = 'Theory')
#plt.plot(axy, 0.5*np.log(1 + (1 - b**2)*axy**2), label = 'Theory')
plt.xlabel('Strength')
plt.ylabel('Mutual Information')
plt.legend(loc = 2)
plt.savefig('linear')
plt.close()
# plt.show()
