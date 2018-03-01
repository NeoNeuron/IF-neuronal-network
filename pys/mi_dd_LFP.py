import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
# sys.argv[1] = index of LFP neuron;
# sys.argv[2] = range of time delay;
# sys.argv[3] = #bins in histogram;
subprocess.call(['./bin/lfp.out', './data/tmp/I.csv', './data/tmp/I1.csv', sys.argv[1], '500,300000', '0.5'])
subprocess.call(['./bin/mi_dd_LFP.out', './data/tmp/I1.csv', './data/tmp/I1.csv', sys.argv[2], sys.argv[3]])
data = pd.read_csv('data/mi/mi_dd_LFP.csv')
plt.plot(data['timelag']*0.5, data['mi'], label = 'mi')
plt.xlabel('Time lags(ms)')
plt.ylabel('mutual info')
plt.legend()
plt.savefig('mi.png')
