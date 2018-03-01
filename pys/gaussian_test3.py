#!/usr/bin/python
# This script aims to find the relationship between baseline level of mutual information and the number of trials when calculating the joint PDF of data.
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import sys
target_ind = 50
lrange = 40
rrange = 0
# num_bin = np.arange(10,51,1)
num_bin = 10
num_trial = np.logspace(3, 8, num = 6)
# print num_trial
a = 0.5
b = 0.5
axy = 1
bls = np.zeros(len(num_trial))

for i in range(len(num_trial)):
  subprocess.call(['python', 'pys/gaussianGenerator.py', str(a), str(b), str(axy), str(int(num_trial[i]))])
  subprocess.call(['./bin/mi_dd.out', './data/tmp/x.csv', './data/tmp/y.csv', str(target_ind), str(lrange) + ',' + str(rrange), str(num_bin), str(num_bin)])
  mi = pd.read_csv('./data/mi/mi_dd.csv')
  bls[i] = mi['mi'][:-10].mean()

# if a == b:
#   rho2 = axy**2/((1 - a**2)**2 + (1 + a**2)*axy**2)
# else:
#   rho2 = axy**2*(1 - b**2)/((1 - a**2)*(1 - a*b)**2 + axy**2*(1 - (a*b)**2))

# plt.semilogy(mis['timelag'], np.ones(len(mis['timelag']))*(-0.5)*np.log(1-rho2), 'r', label = 'Theory')

plt.loglog(num_trial, bls)
plt.xlabel('Number of trials')
plt.ylabel('Baseline of MI')
plt.grid(True, which = 'both', axis = 'both')
plt.savefig('mi_gauss_baseline')
plt.close()
#plt.show()
