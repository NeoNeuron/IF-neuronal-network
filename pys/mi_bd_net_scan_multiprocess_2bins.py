import numpy as np
import subprocess
import sys
import time
from multiprocessing import Pool
import os
import tempfile as tf
# Define task for subprocess
def process_main(tmpname_spike, tmpname_lfp, t_range, dt):
    drange = '0,2'
    #binsize = 0.1
    threshold = 0.018
    tmpname_mi = tf.mkstemp(suffix = '.csv', prefix = 'mi_bd', dir = './data/mi/')[1]
    subprocess.call(['./bin/mi_bd_2bins.out', tmpname_spike, tmpname_lfp, tmpname_mi, drange, str(threshold)])
    data = np.genfromtxt(tmpname_mi, delimiter = ',')
    mi_max = data[:,1].max()
    os.remove(tmpname_mi)
    return mi_max

def process_compact(argv):
    return process_main(*argv)

dt = '0.5'
process_num = int(sys.argv[1])
# prepare data container:

t_range = '5e2,1e8'
mi_max = np.zeros((50,50))
start = time.time()
p = Pool(process_num)
for lfp_ind in range(50):
    tmpname_lfp = './data/lfp/lfp_tmp_' + str(lfp_ind) + '.csv'
    args = [('./data/spike/spike_tmp_' + str(spike_ind) + '.csv', tmpname_lfp, t_range, dt) for spike_ind in range(50)]
    result = p.map(process_compact, args)
    mi_max[:,lfp_ind] = result
tmpname_mis = tf.mkstemp(suffix = '.csv', prefix = 'mi_max', dir = './data/mi/')[1]
np.savetxt(tmpname_mis, mi_max, delimiter = ',', fmt = '%.18f')
print('[-] Output to file -> %s' % tmpname_mis)
finish = time.time()
print('totally cost ', finish - start)
