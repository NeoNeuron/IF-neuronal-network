#!/usr/bin/python
# This script aims to study the relative paremeters that can effects the significance level of mutual information, including the number of bins of probability distribution functions and the size of data pool.
import numpy as np
import pandas as pd
import subprocess
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
target_ind = 50
lrange = 40
rrange = 0
a = 0.5
b = 0.5
axy = 1
# define test_flag, True for #bin test, False for sizeof data test;
test_flag = True
# first fix the size of data pool:
if test_flag:
  num_bin = np.arange(10,51,1)
  num_trial = 100000
  subprocess.call(['python', 'pys/gaussianGenerator.py', str(a), str(b), str(axy), str(num_trial)])
  bls = np.zeros(len(num_bin))
  for i in range(len(num_bin)):
    subprocess.call(['./bin/mi_dd.out', './data/tmp/x.csv', './data/tmp/y.csv', str(target_ind), str(lrange) + ',' + str(rrange), str(num_bin[i]), str(num_bin[i])])
    mi = pd.read_csv('./data/mi/mi_dd.csv')
    bls[i] = mi['mi'][:-10].mean()
  np.savetxt('bls_bins.txt', bls, fmt = '%.18e')
# then fix the number of bins;
else:
  num_bin = 10
  num_trial = np.logspace(3, 8, num = 6)
  bls = np.zeros(len(num_trial))
  for i in range(len(num_trial)):
    subprocess.call(['python', 'pys/gaussianGenerator.py', str(a), str(b), str(axy), str(int(num_trial[i]))])
    subprocess.call(['./bin/mi_dd.out', './data/tmp/x.csv', './data/tmp/y.csv', str(target_ind), str(lrange) + ',' + str(rrange), str(num_bin), str(num_bin)])
    mi = pd.read_csv('./data/mi/mi_dd.csv')
    bls[i] = mi['mi'][:-10].mean()
  np.savetxt('bls_length.txt', bls, fmt = '%.18e')

if test_flag:
  plt.semilogy(1.0/num_bin, bls)
  plt.xlabel('Number of Bins')
  plt.grid(True, which = 'both', axis = 'y')
else:
  plt.loglog(num_trial, bls)
  plt.xlabel('Number of trials')
  plt.grid(True, which = 'both', axis = 'both')

plt.ylabel('Baseline of MI')
plt.savefig('mi_gauss_baseline')
plt.close()
