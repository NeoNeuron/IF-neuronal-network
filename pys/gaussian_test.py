#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import sys
target_ind = 50
lrange = 20
rrange = 20
num_bin = 10
num_trial = 100000
a = float(sys.argv[1])
b = float(sys.argv[2])
axy = float(sys.argv[3])
generator_flag = sys.argv[4]
if generator_flag == 'true':
    subprocess.call(['python', 'pys/gaussianGenerator.py', str(a), str(b), str(axy), str(num_trial)])

subprocess.call(['./bin/mi_dd.out', './data/tmp/x.csv', './data/tmp/y.csv', str(target_ind), str(lrange) + ',' + str(rrange), str(num_bin), str(num_bin)])
mis = pd.read_csv('./data/mi/mi_dd.csv')

if a == b:
    rho2 = axy**2/((1 - a**2)**2 + (1 + a**2)*axy**2)
else:
    rho2 = axy**2*(1 - b**2)/((1 - a**2)*(1 - a*b)**2 + axy**2*(1 - (a*b)**2))

plt.semilogy(mis['timelag'], np.ones(len(mis['timelag']))*(-0.5)*np.log(1-rho2), 'r', label = 'Theory')

plt.semilogy(mis['timelag'], mis['mi'], label = "mi")
#plt.plot(mis['timelag'], mis['mi_shuffle'], label = "surrogate data")
plt.xlabel('Time delay')
plt.ylabel('Mutual Information')
plt.legend()
plt.savefig('mi_gauss')
plt.close()
#plt.show()
