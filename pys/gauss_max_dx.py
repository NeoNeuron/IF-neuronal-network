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
num_trial = 10e7
a = 0.5
b = 0.5
axy = 1
target_ind = 0
bins = np.arange(0.1,1.01,0.02)

x = np.random.normal(size = num_trial)
y1 = np.random.normal(size = num_trial)
y2 = b*y1 + axy*x + np.random.normal(size = num_trial)
np.savetxt('./data/tmp/xn.csv', x, fmt = '%.18f')
np.savetxt('./data/tmp/yn.csv', y2, fmt = '%.18f')

f = open('./data/mi/mi.csv', 'w')
f.close()
for k in bins:
  subprocess.call(['./bin/mi.out', './data/tmp/xn.csv', './data/tmp/yn.csv', str(k), str(k)])

mis = np.loadtxt('./data/mi/mi.csv')
# rho2 = axy**2/((1 - a**2)**2 + (1 + a**2)*axy**2)
rho2 = axy**2/(1+a**2+axy**2)
theory = -0.5*np.log(1-rho2)

plt.plot(bins, mis - theory, label = 'mi~binsize')
plt.xlabel('Binsize of probability density function')
plt.ylabel('Mutual Information')
plt.grid()
plt.savefig('linear')
plt.close()
