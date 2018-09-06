#!/usr/bin/python3
import numpy as np
import subprocess
import sys
import time
from multiprocessing import Pool
import os
import tempfile as tf
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define task for subprocess
def process_main(path, s):
    subprocess.call(['mkdir', '-p', path])
    subprocess.call(['./bin/net_s.out', path, s])
    subprocess.call(['./bin/mi_bd_unity.out', path, '0', '1', '5e2,1e7', '0.5', '15,15', '1.75e-2'])
    print('[-] Main program {' + path + '} finished')
    dat = np.genfromtxt(path + 'mi_bd.csv', delimiter = ',')
    plt.figure(figsize = (10,8), dpi = 80)
    plt.plot(dat[:,0] * 0.5, dat[:,1], label = 'MI')
    plt.plot(dat[:,0] * 0.5, dat[:,2], label = 'MI-shuffle')
    plt.xlabel('Time delay(ms)')
    plt.ylabel('MI')
    plt.legend(fontsize = 20)
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((0,1))
    plt.savefig(path + 'mi_bd.png')
    plt.close()
    return [dat[:,1].max(), dat[:,2].mean()]

def process_compact(argv):
    return process_main(*argv)

process_num = int(sys.argv[1])

# setup parameter ranges;
s_range = np.reshape(np.arange(1e-4, 1e-2, 5e-4), (20,1))

start = time.time()
p = Pool(process_num)
args = [('./data/dataRepo/mi_vs_strength/%.4f/' % s, '%.4f' % s) for s in s_range]
result = p.map(process_compact, args)

# collecting and save data
result = np.array(result)
np.savetxt('./data/dataRepo/mi_vs_strength/mi_vs_strength.csv', np.hstack(s_range, result), delimiter = ',', fmt = '%.8e')
print('[-] Output to file -> ./data/dataRepo/mi_vs_strength/mi_vs_strength.csv')

# fitting with square curve;
def f2(x, a):
    return a*x**2
popt, pcov = curve_fit(f2, s_range, result[:,0])

# plot final results
plt.figure(figsize = (10,8), dpi = 80)
plt.plot(s_range, result[:,0], '*', label = 'MI_max')
plt.plot(s_range, result[:,1], '*', label = 'MI_shuffle')
plt.plot(s_range, f2(s_range, popt[0]), label = 'fitting curve: y=ax^2')
plt.xlabel('Synaptic Strength', fontsize = 20)
plt.ylabel('Max mutual information(nats)', fontsize = 20)
plt.legend(fontsize = 20)
plt.grid()
ax = plt.gca()
ax.xaxis.get_major_formatter().set_powerlimits((0,1))
ax.yaxis.get_major_formatter().set_powerlimits((0,1))
plt.savefig('./data/dataRepo/mi_vs_strength/mi_vs_strength.png')
print('[-] Figure to file -> ./data/dataRepo/mi_vs_strength/mi_vs_strength.png')
finish = time.time()
print('totally cost ', finish - start) 
