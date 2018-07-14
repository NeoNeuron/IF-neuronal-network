#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import struct as st

# compile test program
p = subprocess.call(['g++', '-O2', './src/neuron.cpp', 'src/network.cpp', 'src/math_helper.cpp', 'src/connectivity_matrix.cpp', 'test/network_test.cpp', '-o', 'test/network_test.out'])

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
        dt[i] = 0.5/(2**i)
    plt.loglog(dt, dat_mean[:-1], '.', label = 'data')
    c_est = dat_mean[0]/dt[0]**1
    plt.loglog(dt, c_est*dt**1, label = '1st-order')
    c_est = dat_mean[0]/dt[0]**4
    plt.loglog(dt, c_est*dt**4, label = '4th-order')
    plt.legend()
    plt.grid()
    plt.xlabel('Timing step (ms)')
    plt.ylabel('Relative deviation')
    plt.show()

    # clean test file;
    subprocess.Popen('rm -f test/network_test.out tmp/data_network_test.bin tmp/spiketrain.csv', shell = True)
