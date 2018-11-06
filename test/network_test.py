#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import struct as st

# compile test program
p = subprocess.call(['g++', '--std=c++11', '-O2', 'src/neuron.cpp', 'src/network.cpp', 'src/math_helper.cpp', 'src/get-config.cpp', 'test/network_test.cpp', '-o', 'test/network_test.out'])

if p == 0:
    # excecute test program
    subprocess.call(['./test/network_test.out'])
    # import test data and pre-process
    f = open('tmp/data_network_test.bin', 'rb')
    shape = st.unpack('Q'*2, f.read(8*2))
    dat = np.empty(shape)
    for i in range(shape[0]):
        dat[i] = np.array(st.unpack('d'*shape[1], f.read(8*shape[1])))
    f.close()
    dat = np.array([abs(x - dat[-1]) for x in dat])
    dat_mean = dat.mean(1)
    # plot figure
    dt = np.ones(dat.shape[0] - 1)
    for i in range(len(dt)):
        dt[i] = 0.25/(2**i)
    fig, ax = plt.subplots(figsize=(6,5), dpi=72)
    ax.loglog(dt, dat_mean[:-1], 'o', markerfacecolor = 'None', label = 'data')
    c_est = dat_mean[-2]/dt[-2]**1
    ax.plot(dt, c_est*dt**1, label = '1st-order')
    c_est = dat_mean[-2]/dt[-2]**4
    ax.plot(dt, c_est*dt**4, label = '4th-order')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_yticks(np.logspace(-14,-6,num=5))
    plt.legend()
    plt.grid()
    plt.xlabel('Timing step (ms)', fontsize = 12)
    plt.ylabel('Relative deviation', fontsize = 12)
    plt.savefig('network_test.png')

    # clean test file;
    subprocess.Popen('rm -f test/network_test.out tmp/data_network_test.bin tmp/spiketrain.csv', shell = True)
