#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import subprocess

# compile test program
p = subprocess.call(['g++', '--std=c++11', '-O2', './src/neuron.cpp', './src/math_helper.cpp', 'test/neuron_test.cpp', '-o', 'test/neuron_test.out'])

if p == 0:
    # excecute test program
    subprocess.call(['./test/neuron_test.out'])
    # import test data and pre-process
    arr = np.genfromtxt('tmp/data_new.txt', delimiter = ',')
    arr = np.delete(arr, -1, 1)
    for i in range(arr.shape[0]):
        arr[i] = abs(arr[i] - arr[-1])
    arr_mean = arr.sum(1) / arr.shape[1]

    # plot figure
    dt = np.ones(arr.shape[0] - 1)
    for i in range(len(dt)):
        dt[i] = 0.5/(2**i)
    plt.loglog(dt, arr_mean[:-1], '.', label = 'data')
    c_est = arr_mean[0]/dt[0]**1
    plt.loglog(dt, c_est*dt**1, label = '1st-order')
    c_est = arr_mean[0]/dt[0]**4
    plt.loglog(dt, c_est*dt**4, label = '4th-order')
    plt.legend()
    plt.grid()
    plt.xlabel('Timing step (ms)')
    plt.ylabel('Relative deviation')
    plt.show()

    # clean test file;
    subprocess.Popen('rm -f test/neuron_test.out tmp/data_new.txt', shell = True)
