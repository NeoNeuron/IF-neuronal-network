# -*- coding: utf-8 -*-
import numpy as np
import sys
import subprocess
import tempfile as tf
from multiprocessing import Pool, TimeoutError

# define subprocessing functions;
def process_main(path, tmpname_spike, lfp_inds, t_range):
    tmpname_lfp = tf.mkstemp(dir = './data/tmp/')[1]
    tmpname_mi = tf.mkstemp(dir = './data/tmp/')[1]
    dt = '0.5'
    d_range = '5,5'
    program = './bin/mi_bd_2bins.out'
    parameter = '0.0215'
    subprocess.call(['./bin/lfp.out', path + 'I.bin', tmpname_lfp, lfp_inds, t_range, dt])
    subprocess.call([program, tmpname_spike, tmpname_lfp, tmpname_mi, d_range, parameter])
    dat = np.genfromtxt(tmpname_mi, delimiter = ',')
    # deleting tmp files
    subprocess.call(['rm', '-f', tmpname_lfp])
    subprocess.call(['rm', '-f', tmpname_mi])
    return dat[:,1].max()

def process_compact(argv):
    return process_main(*argv)

# inputing parameters:
dt = 0.5
t_range = '5e2,1e7'
path = sys.argv[1]
spike_ind = '0'
lfp_num = int(sys.argv[2])
process_num = int(sys.argv[3])
trial_num = int(sys.argv[4]) 
# setup general parameters:
spike_file = 'raster.csv'
mat_file = 'mat.csv'
lfp_pool = np.arange(1, 101)

## prepare spike train file;
tmpname_spike = tf.mkstemp(dir = './data/tmp/')[1]
subprocess.call(['./bin/spike.out', path + spike_file, tmpname_spike, spike_ind, t_range, str(dt), 'false'])
# import connecting mat file;
mat = np.genfromtxt(path + mat_file, delimiter = ',')
mat = np.delete(mat, -1, 1)
## delayed mutual information setups:
con_rate = np.zeros(trial_num) # first column for mi_max, second for percentage of direct connection in LFP

# start loop:
p = Pool(processes = process_num)
args_in = []
for i in range(trial_num):
    parts = np.random.choice(lfp_pool, lfp_num, replace = False)
    # counting connecting percentage:
    con_rate[i] = mat[parts,int(spike_ind)].sum() * 1.0 / lfp_num
    # prepare input args;
    args_tmp = (path, tmpname_spike, ','.join(str(ele) for ele in parts), t_range)
    args_in.append(args_tmp)
result = p.map(process_compact, args_in)

#print result
dat = np.transpose(np.vstack((con_rate, np.array(result))))
np.savetxt('mi_lfp_reverse.csv', dat, delimiter = ',', fmt = '%.18f')
## removing tmp files;
subprocess.call(['rm', '-f', tmpname_spike])
