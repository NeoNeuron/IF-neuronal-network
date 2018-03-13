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
num_rand = 52
num_trial = np.logspace(3, 7, num = 5)
a = 0.5
b = 0.5
axy = 1
target_ind = 50
bins = 10
f = open('./data/mi/mi.csv', 'w')
f.close()
for k in num_trial:
  x = np.zeros((num_rand, int(k)))
  y = np.zeros((num_rand, int(k)))
  x[0,:] = np.random.normal(size = int(k))
  y[0,:] = np.random.normal(size = int(k))
  for i in range(num_rand - 1):
    x[i + 1,:] = a*x[i,:] + np.random.normal(size = int(k))
    y[i + 1,:] = b*y[i,:] + k*x[i,:] + np.random.normal(size = int(k))
  np.savetxt('./data/tmp/xn.csv', x[target_ind,:], fmt = '%.10f')
  np.savetxt('./data/tmp/yn.csv', y[target_ind + 1,:], fmt = '%.10f')
  subprocess.call(['./bin/mi.out', './data/tmp/xn.csv', './data/tmp/yn.csv', str(bins), str(bins)])

mis = np.loadtxt('./data/mi/mi.csv')
plt.semilogx(num_trial, mis, label = 'Test')
rho2 = axy**2/((1 - a**2)**2 + (1 + a**2)*axy**2)

plt.semilogx(num_trial, -0.5*np.log(1-rho2)*np.ones(len(num_trial)), label = 'Theory')
plt.xlabel('Length of data')
plt.ylabel('Mutual Information')
plt.legend(loc = 2)
plt.savefig('linear')
plt.close()
