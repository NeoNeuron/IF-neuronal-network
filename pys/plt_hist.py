#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import subprocess

data = np.loadtxt('./data/tmp/singleI.csv')
length = len(data)
autolen = 40
nrow = length / autolen
plt.figure(figsize = (15,12), dpi = 80)
plt.subplot(2,2,1)
plt.hist(data, 100, normed = True, histype = 'stepfilled')
plt.title('Histogram of all ensembles')

data = data.reshape((nrow, autolen))
plt.subplot(2,2,2)
plt.hist(data[:, 0], 100, normed = True, histype = 'stepfilled')
plt.title('Historgram of first column')
plt.subplot(2,2,3)
plt.hist(data[:, autolen/2], 100, normed = True, histype = 'stepfilled')
plt.title('Historgram of middle column')
plt.subplot(2,2,4)
plt.hist(data[:, -1], 100, normed = True, histype = 'stepfilled')
plt.title('Historgram of last column')
plt.savefig('I_hist')
plt.close()
