# -*- coding: utf-8 -*-
import numpy as np
import subprocess
import argparse

parser = argparse.ArgumentParser(description = "Prepare the data files of spike train and local field potential 'LFP' for TDMI calculation.")
#parser.add_argument('dt', type = float, nargs = 1, default = 0.5, help = 'timing step of discrete time series of spike train and local field potential')
parser.add_argument('dir', type = str, help = 'directory of raw data')
parser.add_argument('spike_id', type = str, help = 'index of spike train')
parser.add_argument('max_t_range', metavar = 'maximum_time_range', type = str, help = 'maximum of time series')
args = parser.parse_args()
#print args
dt = 0.5
# Preprocessing spike trains and LFPs;
spike_file = 'raster.csv'
lfp_file = 'I.bin'
DIR = args.dir
spike_ind = args.spike_id
str_range = '500,' + args.max_t_range
subprocess.call(['./bin/spike.out', DIR + spike_file, './data/tmp/singleSpike.bin', spike_ind, str_range, str(dt), 'false'])
subprocess.call(['./bin/spike.out', DIR + spike_file, './data/tmp/singleSpike_shuffle.bin', spike_ind, str_range, str(dt), 'true'])
subprocess.call(['./bin/lfp.out', DIR + lfp_file, './data/tmp/singleI.bin'])
