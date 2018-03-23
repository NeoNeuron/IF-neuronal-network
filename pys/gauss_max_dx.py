#!/usr/bin/python
# This script aims to study the dependence of maximum mutual information of the parematric model of Gaussian random variables;
# The results can be ploted as functions of mutual interacting strength, for both theoretical value and numerical results.
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import sys
num_rand = 2
#num_trial = 10e5*np.logspace(0, 10, num = 11, base = 2)
num_trial = 10e6
a = 0.5
b = 0.5
axy = 1
target_ind = 0
bins = np.arange(10,80,2)

f = open('./data/mi/mi.csv', 'w')
f.write('x_binsize,y_binsize,mi\n')
f.close()
for k in bins:
  x = np.zeros((num_rand, num_trial))
  y = np.zeros((num_rand, num_trial))
  x[0,:] = np.random.normal(size = num_trial)
  y[0,:] = np.random.normal(size = num_trial)
  for i in range(num_rand - 1):
    x[i + 1,:] = a*x[i,:] + np.random.normal(size = num_trial)
    y[i + 1,:] = b*y[i,:] + axy*x[i,:] + np.random.normal(size = num_trial)
  np.savetxt('./data/tmp/xn.csv', x[target_ind,:], fmt = '%.10f')
  np.savetxt('./data/tmp/yn.csv', y[target_ind + 1,:], fmt = '%.10f')
  subprocess.call(['./bin/mi.out', './data/tmp/xn.csv', './data/tmp/yn.csv', str(k), str(k)])

mis = pd.read_csv('./data/mi/mi.csv')
# rho2 = axy**2/((1 - a**2)**2 + (1 + a**2)*axy**2)
rho2 = axy**2/(1+a**2+axy**2)
theory = -0.5*np.log(1-rho2)

plt.plot(mis['x_binsize'], mis['mi'] - theory, label = 'mi~XBinsize')
plt.plot(mis['y_binsize'], mis['mi'] - theory, label = 'mi~YBinsize')
plt.xlabel('Binsize of probability density function')
plt.ylabel('Mutual Information')
plt.legend()
plt.savefig('linear')
plt.close()
