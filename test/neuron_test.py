#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import subprocess

# compile test program
p = subprocess.call(['g++', '--std=c++11', '-O2', './src/neuron.cpp', './src/math_helper.cpp', 'test/neuron_test.cpp', '-o', 'test/neuron_test.out'])

if p == 0:
    # clear previous data files
    subprocess.Popen('rm -f tmp/data_neuron_raster.csv tmp/data_new.txt', shell = True)
    # excecute test program
    subprocess.call(['./test/neuron_test.out'])
    # import test data and pre-process
    dat = np.genfromtxt('tmp/data_new.txt', delimiter = ',')
    dat = np.delete(dat, -1, 1)
    # import spike trains;
    dat_spike = np.genfromtxt('tmp/data_neuron_raster.csv', delimiter = ',')
    dat_spike = dat_spike[:,:-1]
    #============================================================ 
    dat = np.array([abs(x - dat[-1]) for x in dat])
    dat_mean = dat.mean(1)
    dat_spike = np.array([abs(x - dat_spike[-1]) for x in dat_spike])
    dat_spike_mean = dat_spike.mean(1)

    #============================================================ 
    # plot figure
    #============================================================ 
    dt = np.ones(dat.shape[0] - 1)
    for i in range(len(dt)):
        dt[i] = 0.5/(2**i)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12,6), dpi = 72)
    # subplot 1 ==================
    ax1.plot(dt, dat_mean[:-1], 'o', markerfacecolor = 'None', label = 'data')
    c_est = dat_mean[-2]/dt[-2]**1
    ax1.plot(dt, c_est*dt**1, label = '1st-order')
    c_est = dat_mean[-2]/dt[-2]**4
    ax1.plot(dt, c_est*dt**4, label = '4th-order')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    #ax1.set_yticks(np.logspace(-14,-6,num=5))
    ax1.legend()
    ax1.grid()
    ax1.set_title('Convergence of potential')
    ax1.set_xlabel('Timing step (ms)', fontsize = 12)
    ax1.set_ylabel('Relative deviation', fontsize = 12)
    # subplot 2 ==================
    ax2.plot(dt, dat_spike_mean[:-1], 'o', markerfacecolor = 'None', label = 'data_spike')
    c_est = dat_spike_mean[-2]/dt[-2]**1
    ax2.plot(dt, c_est*dt**2, label = '2nd-order')
    c_est = dat_spike_mean[-2]/dt[-2]**4
    ax2.plot(dt, c_est*dt**4, label = '4th-order')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    #ax2.set_yticks(np.logspace(-14,-6,num=5))
    ax2.set_title('Convergence of Spike trains')
    plt.legend()
    plt.grid()
    ax2.set_xlabel('Timing step (ms)', fontsize = 12)
    ax2.set_ylabel('Relative deviation', fontsize = 12)
    plt.show()

    # clean test file;
    #subprocess.Popen('rm -f test/neuron_test.out tmp/data_new.txt', shell = True)
