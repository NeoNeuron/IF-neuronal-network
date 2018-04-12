import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
spike_ind = sys.argv[1]
lfp_ind = sys.argv[2]
dt = 0.5
ntd = 0
ptd = 0
# tmax_list = np.arange(30000, 300001, 30000)
tmax_list = 1e5*np.logspace(12, 13, num = 2, base = 2)
binsize = 0.03
# prepare data container:
mi_max = np.ndarray(len(tmax_list))

t_range = [500, 0]
drange = str(ntd) + ',' + str(ptd)
i = 0
counter = 12
for tmax in tmax_list:
  t_range[1] = int(tmax)
  str_range = str(t_range[0]) + ',' + str(t_range[1])
  subprocess.call(['./bin/spike.out', './data/dataRepo/1.3kHz/raster.csv', './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
  subprocess.call(['./bin/lfp.out', './data/dataRepo/1.3kHz/I.bin', './data/tmp/singleI.bin', lfp_ind, str_range, str(dt)])
  subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike.bin', './data/tmp/singleI.bin', drange, str(dt), str(binsize)])
  subprocess.call(['cp', './data/mi/joint_pdf.csv', './data/mi/joint_pdf_' + str(counter) + '.csv'])
  data = pd.read_csv('data/mi/mi_bd.csv')
  mi_max[i] = data['mi'].max()
  counter += 1
  i += 1
np.savetxt('./data/tmp/mi_max.csv', mi_max, delimiter = ',', fmt = '%.18f')
