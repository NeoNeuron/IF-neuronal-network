import numpy as np
import pandas as pd
import subprocess
import sys
import time
import multiprocessing
import os
import tempfile as tf
#from myio import uni_num 
# Define task for subprocess
def process_main(path, spike_ind, tmpname_lfp, str_range, dt):
    drange = '0,0'
    #binsize = 0.1
    threshold = 0.1
    tmpname_spike = tf.mkstemp(suffix = '.bin', prefix = 'singleSpike', dir = './data/spike/')[1]
    tmpname_mi = tf.mkstemp(suffix = '.csv', prefix = 'mi_bd', dir = './data/mi/')[1]
    subprocess.call(['./bin/spike.out', path + 'raster.csv', tmpname_spike, str(spike_ind), str_range, str(dt), 'false'])
    subprocess.call(['./bin/mi_bd_2bins.out', tmpname_spike, tmpname_lfp, tmpname_mi, drange, str(threshold)])
    data = pd.read_csv(tmpname_mi)
    mi_max = data['mi'].max()
    os.remove(tmpname_spike)
    os.remove(tmpname_mi)
    return mi_max

dt = 0.5
path = sys.argv[1]
process_num = int(sys.argv[2])
# prepare data container:

t_range = [500, 1e7]
str_range = str(t_range[0]) + ',' + str(t_range[1])
mi_max = np.zeros((100,100))
start = time.time()
tmpname_lfp = tf.mkstemp(suffix = '.bin', prefix = 'singleI', dir = './data/lfp/')[1]
for lfp_ind in range(100):
  subprocess.call(['./bin/lfp.out', path + 'I.bin', tmpname_lfp, str(lfp_ind), str_range, str(dt)])
  p = multiprocessing.Pool(process_num)
  result = [p.apply_async(func = process_main, args=(path, spike_ind, tmpname_lfp, str_range, dt)) for spike_ind in range(100)] 
  p.close()
  p.join()
  i = 0
  for res in result:
      mi_max[i][lfp_ind] = (res.get())
      i += 1
tmpname_mis = tf.mkstemp(suffix = '.csv', prefix = 'mi_max', dir = './data/mi/')[1]
np.savetxt(tmpname_mis, mi_max, delimiter = ',', fmt = '%.18f')
print('[-] Output to file -> %s' % tmpname_mis)
finish = time.time()
print('totally cost ', finish - start)
