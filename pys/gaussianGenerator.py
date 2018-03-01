#!/usr/bin/python
import numpy as np
import sys
num_rand = 52
num_trial = int(sys.argv[4])
x = np.zeros((num_rand, num_trial))
y = np.zeros((num_rand, num_trial))
x[0,:] = np.random.normal(size = num_trial)
y[0,:] = np.random.normal(size = num_trial)
ax = float(sys.argv[1])
ay = float(sys.argv[2])
k = float(sys.argv[3])
for i in range(num_rand - 1):
    x[i + 1,:] = ax*x[i,:] + np.random.normal(size = num_trial)
    y[i + 1,:] = ay*y[i,:] + k*x[i,:] + np.random.normal(size = num_trial)
np.savetxt('./data/tmp/x.csv', x, fmt = '%.10f', delimiter = ',')
np.savetxt('./data/tmp/y.csv', y, fmt = '%.10f', delimiter = ',')
