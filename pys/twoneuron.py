#!/usr/bin/python
import subprocess
import os
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd

# Define loading directory and figure saving directory:
saving_dir = './data/tmp/'

# define trial times;
num_trial = 2

# prepare output files
f = open('./data/tmp/n1i.csv', 'w')
f.close()
f = open('./data/tmp/n2i.csv', 'w')
f.close()
f = open('./data/tmp/n1s.csv', 'w')
f.close()
f = open('./data/tmp/n2s.csv', 'w')
f.close()
for i in range(num_trial):
    subprocess.call(['./bin/net.out', saving_dir])
    subprocess.call(['./bin/split.out', './data/tmp/raster.csv', './data/tmp/n1s.csv', './data/tmp/n2s.csv'])
    subprocess.call(['./bin/transpose.out', './data/tmp/I.csv', './data/tmp/dataT.csv'])
    subprocess.call(['./bin/split.out', './data/tmp/dataT.csv', './data/tmp/n1i.csv', './data/tmp/n2i.csv'])
# Setting loops for local field potentials;
# target_neuron_indice_list = random.sample(all_neuron, 1)
