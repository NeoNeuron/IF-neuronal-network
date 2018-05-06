import numpy as np
import pandas as pd
import subprocess
import sys
spike_ind = sys.argv[1]
lfp_ind = sys.argv[2]
dt = 0.5
ntd = 0
ptd = 0
drange = str(ntd) + ',' + str(ptd)
tmax = 1e8
binsize = np.arange(0.001, 0.0301, 0.001)
# prepare data container:
mi_max = np.ndarray(len(binsize))

t_range = [500, tmax]
str_range = str(t_range[0]) + ',' + str(t_range[1])
i = 0
subprocess.call(['./bin/spike.out', './data/dataRepo/1.3kHz/raster.csv', './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
subprocess.call(['./bin/lfp.out', './data/dataRepo/1.3kHz/I.bin', './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
for bin in binsize:
  subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike.bin', './data/tmp/singleI.bin', 'data/mi/mi_bd.csv', drange, str(bin)])
  data = pd.read_csv('data/mi/mi_bd.csv')
  mi_max[i] = data['mi'].max()
  i += 1
np.savetxt('./data/tmp/mi_max_dx.csv', mi_max, delimiter = ',', fmt = '%.18f')
