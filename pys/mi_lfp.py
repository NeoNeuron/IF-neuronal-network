# -*- coding: utf-8 -*-
import numpy as np
import sys
import subprocess
import random
import tempfile as tf
from multiprocessing import Pool, TimeoutError
import os
import argparse

# define subprocessing functions;
def process_main(path, tmpname_spike, lfp_inds, t_range):
    tmpname_lfp = tf.mkstemp()[1]
    tmpname_mi = tf.mkstemp()[1]
    dt = '0.5'
    d_range = '5,5'
    program = './bin/mi_bd.out'
    parameter = '0.03'
    subprocess.call(['./bin/lfp.out', path + 'I.bin', tmpname_lfp, lfp_inds, t_range, dt])
    subprocess.call([program, tmpname_spike, tmpname_lfp, tmpname_mi, d_range, '0.03'])
    data = np.genfromtxt(tmpname_mi, delimiter = ',')
    mi_max = data[:,1].max()
    # deleting tmp files
    os.remove(tmpname_lfp)
    os.remove(tmpname_mi)
    return mi_max

def process_compact(argv):
    return process_main(*argv)

# inputing parameters:
dt = 0.5
tmax = 1e7
t_range = '500,' + str(int(tmax))
path = sys.argv[1]
spike_ind = '0'
lfp_num = int(sys.argv[2])
process_num = int(sys.argv[3])
trial_num = int(sys.argv[4]) 
# setup general parameters:
spike_file = 'raster.csv'
#lfp_file = 'I.bin'
mat_file = 'mat.csv'
lfp_pool = range(1, 101)

## prepare spike train file;
tmpname_spike = tf.mkstemp()[1]
subprocess.call(['./bin/spike.out', path + spike_file, tmpname_spike, spike_ind, t_range, str(dt), 'false'])
# import connecting mat file;
mat = np.genfromtxt(path + mat_file, delimiter = ',')
mat = np.delete(mat, -1, 1)
## delayed mutual information setups:
#program = './bin/mi_bd.out'
#param = 0.03
mi_para = np.zeros((trial_num, 2)) # first column for mi_max, second for percentage of direct connection in LFP

# start loop:
p = Pool(processes = process_num)
#result = []
args_in = []
for i in range(trial_num):
    parts = random.sample(lfp_pool, lfp_num)
    # counting connecting percentage:
    counts = 0
    for j in parts:
        if mat[0][j] == 1:
            counts += 1
    mi_para[i][0] = counts * 1.0 / lfp_num
    t_range = '500,' + str(int(tmax))
    args_tmp = (path, tmpname_spike, ','.join(str(ele) for ele in parts), t_range)
    args_in.append(args_tmp)
#    result.append(p.apply_async(func = process_main, args=(path, tmpname_spike, lfp_inds), t_range)) 
#p.close()
#p.join()
result = p.map(process_compact, args_in)
#print result
#i = 0
## collecting results;
#for res in result:
#    mi_para[i][1] = res.get()
#    i += 1
mi_para[:,1] = result
np.savetxt('mi_lfp.csv', mi_para, delimiter = ',', fmt = '%.18f')
## removing tmp files;
os.remove(tmpname_spike)
