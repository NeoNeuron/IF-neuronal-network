#!/bin/python
import pandas as pd
import matplotlib.pyplot as plt
import sys

sta = pd.read_csv('data/sta.csv')
mi = pd.read_csv('data/mi/mi_sp.csv')
lcc = pd.read_csv('data/lcc/lcc.csv')
plt.figure(figsize = (15,5), dpi = 75)
plt.subplot(1,3,1)
plt.plot(sta['timelag'], sta['sta'])
plt.xlabel('time-delay(ms)')
plt.ylabel('Potential')
plt.grid(True)
plt.subplot(1,3,2)
plt.plot(lcc['timelag'], lcc['lcc'])
plt.xlabel('time-delay(ms)')
plt.ylabel('Pearson\' Correlation')
plt.grid(True)
plt.subplot(1,3,3)
plt.plot(mi['timelag'], mi['ordered'])
plt.plot(mi['timelag'], mi['random'])
plt.xlabel('time-delay(ms)')
plt.ylabel('MI(bits)')
plt.grid(True)
# plt.savefig('./reports/Jul20/' + sys.argv[1])
plt.show()
