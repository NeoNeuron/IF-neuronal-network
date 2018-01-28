import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
t_range = '500,600000'
dt = 1
spike_ind1 = sys.argv[2]
spike_ind2 = sys.argv[3]
subprocess.call(['./bin/spike.out', sys.argv[1], './data/tmp/singleSpike1.csv', spike_ind1, t_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', sys.argv[1], './data/tmp/singleSpike2.csv', spike_ind2, t_range, str(dt), 'false'])
subprocess.call(['./bin/mi_bb.out', './data/tmp/singleSpike1.csv', './data/tmp/singleSpike2.csv', '20,20'])
data = pd.read_csv('data/mi/mi_bb.csv')
plt.plot(data['timelag'], data['mi'], label = 'mi')
plt.xlabel('Time lags')
plt.ylabel('mutual info')
plt.legend()
ax = plt.gca()
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.savefig('mi.png')
plt.close()
