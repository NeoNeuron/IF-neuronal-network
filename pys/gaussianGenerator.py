#!/usr/bin/python
import numpy as np
num_rand = 200
num_trial = 10000
x = np.zeros((num_rand, num_trial))
y = np.zeros((num_rand, num_trial))
x[0,:] = np.random.normal(size = num_trial)
y[0,:] = np.random.normal(size = num_trial)
ax = -0.1
ay = -0.1
# axy = np.arange(0,1,0.05)
# for k in axy:
k = 0.5
for i in range(num_rand - 1):
    x[i + 1,:] = ax*x[i,:] + np.random.normal(size = num_trial)
    y[i + 1,:] = ay*y[i,:] + k*x[i,:] + np.random.normal(size = num_trial)
np.savetxt('./data/tmp/x.csv', x, fmt = '%.10f', delimiter = ',')
np.savetxt('./data/tmp/y.csv', y, fmt = '%.10f', delimiter = ',')
