# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import multiprocessing
import argparse
import configparser as cp
import time
import os

# Define task for subprocess
def process_main(path, spike_id, lfp_name, delay_range, time_range, dt):
    binsize = '0.005'
    if spike_id < 10:
        spike_name = path + '/spike0' + str(spike_id) + '.bin'
        spike_name_shuffle = path + '/spike0' + str(spike_id) + 'shuffle.bin'
        mi_name = path + '/mi0' + str(spike_id) + '.csv'
    else:
        spike_name = path + '/spike' + str(spike_id) + '.bin'
        spike_name_shuffle = path + '/spike' + str(spike_id) + 'shuffle.bin'
        mi_name = path + '/mi' + str(spike_id) + '.csv'
    subprocess.call(['./bin/spike.out', path + '/raster.csv', spike_name, str(spike_id), time_range, str(dt), 'false'])
    subprocess.call(['./bin/spike.out', path + '/raster.csv', spike_name_shuffle, str(spike_id), time_range, str(dt), 'true'])
    subprocess.call(['./bin/mi_bd.out', spike_name, lfp_name, mi_name, delay_range, binsize])
    data = np.genfromtxt(mi_name, delimiter = ',')
    subprocess.call(['./bin/mi_bd.out', spike_name_shuffle, lfp_name, mi_name, delay_range, binsize])
    data_shuffle = np.genfromtxt(mi_name, delimiter = ',')
    mi_bl = data_shuffle[:,1].mean()
    #mi_max_delay = data[data[:,1].argmax(), 0]
    data[:,1] /= mi_bl
    os.remove(spike_name)
    os.remove(spike_name_shuffle)
    os.remove(mi_name)
    return data[:,1]

parser = argparse.ArgumentParser(description = "Prepare the data files of spike train and local field potential 'LFP' for TDMI calculation.")
parser.add_argument('dir', type = str, help = 'directory of source data and output data')
parser.add_argument('id_max', type = int, help = 'maximum range of the spike indices')
parser.add_argument('td', metavar = 'time delay', type = str, help = 'negative and positive maximum delay time, seperated by comma')
parser.add_argument('pn', metavar = 'process number', type = int, help = 'number of processes')
args = parser.parse_args()

# Preprocessing spike trains and LFPs;
config = cp.ConfigParser()
config.read('doc/config_lfp.ini')
lfp_type = config.get('LFP', 'LFPType')
spike_file = 'raster.csv'
lfp_name = 'lfp_' + lfp_type + '.bin'
time_range = config.get('Time', 'TimeRangeMin') + ',' + config.get('Time', 'TimeRangeMax')
dt = float(config.get('Time', 'SamplingTimingStep'))
# Calculating Mutual information:
# Calculation parameters; binsize for mi_bd.out, threshold for mi_bd_2bins.out
program = './bin/mi_bd.out'
parameter = '0.005'
ncols = np.array([int(i) for i in args.td.split(',')])
ncols = ncols.sum() + 1
mi_scan_res = np.empty((args.id_max, ncols))

# start main process
start = time.time()
p = multiprocessing.Pool(args.pn)
result = [p.apply_async(func = process_main, args=(args.dir, spike_id, args.dir + lfp_name, args.td, time_range, dt)) for spike_id in range(args.id_max)] 
p.close()
p.join()
i = 0
for res in result:
  mi_scan_res[i] = (res.get())
  i += 1
np.savetxt(args.dir + "/mi_lfp_" + lfp_type + "_scan.csv", mi_scan_res, delimiter = ',')
print('[-] Output to file -> ' + args.dir + "/mi_lfp_" + lfp_type + "_scan.csv")
finish = time.time()
print('totally cost %3.3f s.' % (finish - start))
