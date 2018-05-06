import numpy as np
import pandas as pd
import subprocess
import sys
import time
import io
dt = 0.5
ntd = 0
ptd = 0
binsize = 0.03
path = sys.argv[1]
# prepare data container:

t_range = [500, 1e7]
drange = str(ntd) + ',' + str(ptd)
str_range = str(t_range[0]) + ',' + str(t_range[1])
mi_max = np.zeros((100,100))
start = time.time()
name_num = uni_num()
tmpname_spike = 'singleSpike' + name_num + '.bin'
tmpname_lfp = 'singleI' + name_num + '.bin'
for lfp_ind in range(100):
  subprocess.call(['./bin/lfp.out', path + 'I.bin', './data/tmp/' + tmpname_lfp, str(lfp_ind), str_range, str(dt)])
  for spike_ind in range(100):
		if spike_ind != lfp_ind:
			subprocess.call(['./bin/spike.out', path + 'raster.csv', './data/tmp/' + tmpname_spike, str(spike_ind), str_range, str(dt), 'false'])
			subprocess.call(['./bin/mi_bd.out', './data/tmp/' + tmpname_spike, './data/tmp/' + tmpname_lfp, drange, str(binsize)])
			data = pd.read_csv('data/mi/mi_bd.csv')
			mi_max[spike_ind][lfp_ind] = data['mi'].max()
np.savetxt('./data/mi/mi_max' + name_num + '.csv', mi_max, delimiter = ',', fmt = '%.18f')
finish = time.time()
print('totally cost ', finish - start)
