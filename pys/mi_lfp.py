import numpy as np
import pandas as pd
import sys
import subprocess
import random
import tempfile as tf
import multiprocessing
import os

# define subprocessing functions;
def process_main(path, tmpname_spike, tmpname_spike_shuffle, lfp_inds, t_range):
    tmpname_lfp = tf.mkstemp()[1]
    tmpname_mi = tf.mkstemp()[1]
    dt = '0.5'
    d_range = '5,5'
    program = './bin/mi_bd.out'
    parameter = '0.03'
    subprocess.call(['./bin/lfp.out', path + 'I.bin', tmpname_lfp, lfp_inds, t_range, dt])
    subprocess.call([program, tmpname_spike, tmpname_lfp, tmpname_mi, d_range, '0.03'])
    data = pd.read_csv(tmpname_mi)
    mi_max = data['mi'].max()
    subprocess.call([program, tmpname_spike_shuffle, tmpname_lfp, tmpname_mi, d_range, '0.03'])
    data = pd.read_csv(tmpname_mi)
    mi_mean = data['mi'].mean()
    # deleting tmp files
    os.remove(tmpname_lfp)
    os.remove(tmpname_mi)
    return {'mi_max':mi_max, 'mi_mean':mi_mean} 

# inputing parameters:
path = sys.argv[1]
spike_ind = sys.argv[2]
lfp_num = int(sys.argv[3])
process_num = int(sys.argv[4])
trial_num = int(sys.argv[5]) 
# setup general parameters:
dt = 0.5
tmax = 1e7
t_range = '500,' + str(int(tmax))
spike_file = 'raster.csv'
#lfp_file = 'I.bin'
mat_file = 'mat.csv'
lfp_pool = range(1, 100)

## prepare spike train file;
tmpname_spike = tf.mkstemp()[1]
tmpname_spike_shuffle = tf.mkstemp()[1]
subprocess.call(['./bin/spike.out', path + spike_file, tmpname_spike, spike_ind, t_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', path + spike_file, tmpname_spike_shuffle, spike_ind, t_range, str(dt), 'true'])
# import connecting mat file;
mat = np.genfromtxt(path + mat_file, delimiter = ',')
mat = np.delete(mat, -1, 1)
## delayed mutual information setups:
#program = './bin/mi_bd.out'
#param = 0.03
mi_para = np.zeros((trial_num, 2 + lfp_num)) # first column for mi_max, second for mi_mean, third for percentage of direct connection

# start loop:
p = multiprocessing.Pool(process_num)
result = []
for i in range(trial_num):
    parts = random.sample(lfp_pool, lfp_num)
    for j in range(lfp_num):
        mi_para[i][j + 2] = parts[j]
    ## counting connecting percentage:
    #counts = 0
    #for j in parts:
    #    if mat[int(spike_ind)][j] == 1:
    #        counts += 1
    #mi_para[i][2] = counts * 1.0 / lfp_num
    lfp_inds = ','.join(str(ele) for ele in parts)
    result.append(p.apply_async(func = process_main, args=(path, tmpname_spike, tmpname_spike_shuffle, lfp_inds, t_range))) 
p.close()
p.join()
i = 0
# collecting results;
for res in result:
    tmp = res.get()
    mi_para[i][0] = tmp['mi_max']
    mi_para[i][1] = tmp['mi_mean']
    i += 1
np.savetxt('mi_lfp.csv', mi_para, delimiter = ',', fmt = '%.18f')
## removing tmp files;
os.remove(tmpname_spike)
os.remove(tmpname_spike_shuffle)

