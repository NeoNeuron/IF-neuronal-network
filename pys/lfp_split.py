import subprocess
import sys
import time
from multiprocessing import Pool
import os
# Define task for subprocess
def process_main(path, lfp_ind, t_range, dt):
    filename_lfp = './data/lfp/lfp_tmp_' + lfp_ind + '.csv'
    subprocess.call(['./bin/lfp.out', path + 'I.bin', filename_lfp, lfp_ind, t_range, dt])
    return

def process_compact(argv):
    return process_main(*argv)

path = sys.argv[1]
process_num = int(sys.argv[2])
lfp_max_num = int(sys.argv[3])
# prepare data container:
dt = '0.5'
t_range = '5e2,1e8'
start = time.time()
p = Pool(process_num)
args = [(path, str(lfp_ind), t_range, dt) for lfp_ind in range(lfp_max_num)]
p.map(process_compact, args)
finish = time.time()
print('totally cost ', finish - start)
