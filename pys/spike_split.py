import subprocess
import sys
import time
from multiprocessing import Pool
import os
# Define task for subprocess
def process_main(path, spike_ind, t_range, dt):
    filename_spike = './data/spike/spike_tmp_' + spike_ind + '.csv'
    subprocess.call(['./bin/spike.out', path + 'raster.csv', filename_spike, spike_ind, t_range, dt, 'false'])
    return

def process_compact(argv):
    return process_main(*argv)

path = sys.argv[1]
process_num = int(sys.argv[2])
spike_max_num = int(sys.argv[3])
# prepare data container:
dt = '0.5'
t_range = '5e2,1e8'
start = time.time()
p = Pool(process_num)
args = [(path, str(spike_ind), t_range, dt) for spike_ind in range(spike_max_num)]
p.map(process_compact, args)
finish = time.time()
print('totally cost ', finish - start)
