#!/bin/python
import pandas as pd
import matplotlib.pyplot as plt

sta = pd.read_csv('data/sta.csv')
mi = pd.read_csv('data/mi/mi_sl.csv')
plt.figure(figsize = (10,5), dpi = 75)
plt.subplot(1,2,1)
plt.plot(sta['timelag'], sta['sta'])
plt.grid(True)
plt.subplot(1,2,2)
plt.plot(mi['timelag'], mi['ordered'])
plt.grid(True)
plt.show()
