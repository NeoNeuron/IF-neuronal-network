#!/usr/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
ntd = 10
ptd = 10
summi = np.zeros(ntd + ptd + 1)
for i in range(20):
    ind = i + 10
    subprocess.call(['./bin/mi_bd.out', 'data/tmp/singleSpike.csv', 'data/tmp/singleI.csv', str(ind), str(ntd) + ',' + str(ptd), '500'])
    data = pd.read_csv('data/mi/mi_bd.csv')
    singlemi = data['mi'] / data['mi'].max()
    summi += singlemi
summi_shuffle = np.zeros(ntd + ptd + 1)
for i in range(20):
    ind = i + 10
    subprocess.call(['./bin/mi_bd.out', 'data/tmp/singleSpike.csv', 'data/tmp/singleI_shuffle.csv', str(ind), str(ntd) + ',' + str(ptd), '500'])
    data = pd.read_csv('data/mi/mi_bd.csv')
    singlemi_shuffle = data['mi'] / data['mi'].max()
    summi_shuffle += singlemi_shuffle
plt.plot(data['timelag']*0.5, summi, label = 'mi')
plt.plot(data['timelag']*0.5, summi_shuffle, label = 'mi_shuffle')
plt.xlabel('timelag(ms)')
plt.ylabel('mutual info')
plt.legend()
plt.show()
