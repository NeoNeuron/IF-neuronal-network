import numpy as np
import matplotlib.pyplot as plt
import argparse
import configparser as cp
import subprocess
import subprocess

p = subprocess.call(['g++', '--std=c++11', '-O2', 'src/neuron.cpp', 'src/math_helper.cpp', 'src/get-config.cpp', 'src/main_gain_curve.cpp', '-o', 'bin/gain_curve.out'])

if p == 0:
    # excecute test program
    subprocess.call(['./bin/gain_curve.out', './tmp/'])
    #============================================================ 
    # load configureations: 
    #============================================================ 
    config = cp.ConfigParser()
    config.read('doc/gain_curve.ini')
    neuron_type = config.get('Neuron', 'NeuronType')
    p_rate = float(config.get('Driving Settings', 'DrivingRate'))
    s_min  = float(config.get('Driving Settings', 'DrivingStrengthMin'))
    s_max  = float(config.get('Driving Settings', 'DrivingStrengthMax'))
    ds     = float(config.get('Driving Settings', 'ds'))
    tmax   = float(config.get('Time', 'MaximumTime'))

    #============================================================ 
    # import test data and pre-process
    #============================================================ 
    dat = np.array([])
    # import spike trains;
    f = open('./tmp/gain_curve.csv')
    for line in f:
        line = line.replace('\n', '').split(',')
        line.remove('')
        dat = np.append(dat, len(line)/tmax*1e3)
    f.close()
    #============================================================ 
    # plot figure
    #============================================================ 
    s = np.arange(s_min, s_max + ds, ds)
    if len(s) > len(dat):
        s = s[0:len(dat)]
    fs = s*p_rate*1e3
    fig, ax = plt.subplots(figsize = (7,6), dpi = 60)
    ax.plot(fs, dat)
    ax.set_xlabel(r'$f\nu$ (Hz)')
    ax.set_ylabel('Firing Rate (Hz)')
    ax.xaxis.get_major_formatter().set_powerlimits((0,1))
    ax.set_title('Gain Curve of ' + neuron_type + ' Neuron (' + str(p_rate) + ' kHz)')
    ax.grid(linestyle = '--')
    fig.savefig('gain_curve.eps')
    plt.close()

    # clean test file;
    #subprocess.Popen('rm -f test/network_test.out tmp/data_network_test.bin tmp/spiketrain.csv', shell = True)
