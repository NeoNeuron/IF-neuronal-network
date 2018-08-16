# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import argparse
import sys
parser = argparse.ArgumentParser(description = "Calculate TDMI between given spike train and local field potential 'LFP'")
#parser.add_argument('dt', type = float, nargs = 1, default = 0.5, help = 'timing step of discrete time series of spike train and local field potential')
parser.add_argument('prog', type = str, nargs = 1, help = 'Choice of TDMI calculation program')
parser.add_argument('para', metavar = 'parameter', type = str, nargs = 1, help = 'parameter of TDMI programs')
parser.add_argument('ntd', metavar = 'NTD', type = str, nargs = 1, help = 'negative maximum delay time')
parser.add_argument('ptd', metavar = 'PTD', type = str, nargs = 1, help = 'positive maximum delay time')
parser.add_argument('plt', metavar = 'plot_handle', type = str, nargs = 1, help = 'handle of plotting results')
args = parser.parse_args()
#print args
dt = 0.5
# Calculating Mutual information:
drange = args.ntd[0] + ',' + args.ptd[0]
# Calculation parameters; binsize for mi_bd.out, threshold for mi_bd_2bins.out
program = args.prog[0]
parameter = args.para[0]
subprocess.call([program, './data/tmp/singleSpike.bin', './data/tmp/singleI.bin', './data/mi/mi_bd.csv', drange, parameter])
data = np.genfromtxt('data/mi/mi_bd.csv', delimiter = ',')
subprocess.call([program, './data/tmp/singleSpike_shuffle.bin', './data/tmp/singleI.bin', './data/mi/mi_bd_shuffle.csv', drange, parameter])
data_shuffle = np.genfromtxt('data/mi/mi_bd_shuffle.csv', delimiter = ',')
# Plotting:
plotway = args.plt[0]
if plotway == 'linear':
    plt.plot(data[:,0]*dt, data[:,1], label = 'mi')
    plt.plot(data_shuffle[:,0]*dt, data_shuffle[:,1], label = 'mi_shuffle')
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((0,1))
elif plotway == 'log':
    plt.semilogy(data[:,0]*dt, data[:,1], label = 'mi')
    plt.semilogy(data_shuffle[:,0]*dt, data_shuffle[:,1], label = 'mi_shuffle')

plt.xlabel('Delay Time(ms)')
plt.ylabel('mutual info')
plt.grid()
plt.legend()
plt.savefig('mi')
plt.close()
