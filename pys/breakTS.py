#!/usr/bin/python
import numpy as np
import subprocess
import sys

dt = 0.5
spike_ind = sys.argv[1]
lfp_ind = sys.argv[2]

subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike.csv', spike_ind, '500,600000', str(dt)])
subprocess.call(['./bin/lfp.out', './data/tmp/I.csv', './data/tmp/singleI.csv', lfp_ind, '500,600000', str(dt)])
