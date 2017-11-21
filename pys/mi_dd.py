import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
bins = sys.argv[1]
subprocess.call(['./bin/mi_dd.out', './data/tmp/x.csv', './data/tmp/y.csv', '100', '30,30', bins, bins])
data = pd.read_csv('data/mi/mi_dd.csv')
plt.plot(data['timelag'], data['mi'], label = 'mi')
plt.plot(data['timelag'], data['mi_shuffle'], label = 'mi_shuffle')
plt.xlabel('Time lags')
plt.ylabel('mutual info')
plt.legend()
plt.savefig('mi.png')
