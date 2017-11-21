import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
spike_ind = sys.argv[1]
lfp_ind = sys.argv[2]
dt = 0.5
ntd = 15
ptd = 15
tmax_list = np.arange(30000, 600001, 30000)
hist_bin_list = np.arange(5, 105, 5)
# prepare data container:
mi_max = np.ndarray((len(hist_bin_list), len(tmax_list)))
snr = np.ndarray((len(hist_bin_list), len(tmax_list)))
peak_dev = np.ndarray((len(hist_bin_list), len(tmax_list)))

t_range = [0, 0]
t_range[0] = 500
drange = str(ntd) + ',' + str(ptd)
i = 0
for tmax in tmax_list:
    t_range[1] = tmax
    str_range = str(t_range[0]) + ',' + str(t_range[1])
    subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike.csv', spike_ind, str_range, str(dt), 'false'])
    subprocess.call(['./bin/spike.out', './data/tmp/raster.csv', './data/tmp/singleSpike_shuffle.csv', spike_ind, str_range, str(dt), 'true'])
    subprocess.call(['./bin/lfp.out', './data/tmp/I.csv', './data/tmp/singleI.csv', lfp_ind, str_range, str(dt)])
    j = 0
    for hist_bin in hist_bin_list:
# hist_bin = 50
        subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike.csv', './data/tmp/singleI.csv', drange, str(dt), str(hist_bin)])
        data = pd.read_csv('data/mi/mi_bd.csv')
        mi_max[j][i] = data['mi'].max()
# plt.plot(data['timelag']*dt, data['mi'], label = 'mi')
        subprocess.call(['./bin/mi_bd.out', './data/tmp/singleSpike_shuffle.csv', './data/tmp/singleI.csv', drange, str(dt), str(hist_bin)])
        data = pd.read_csv('data/mi/mi_bd.csv')
        snr[j][i] = mi_max[j][i] / (data['mi'].mean()+ data['mi'].std() * 3)
        peak_dev[j][i] = (mi_max[j][i] - data['mi'].mean()) / data['mi'].std()
# plt.plot(data['timelag']*dt, data['mi'], label = 'mi_shuffle')
        j += 1
    i += 1
# np.savetxt('baselines.csv', baselines, delimiter = ',', fmt = '%.10f')
np.savetxt('mi_max.csv', mi_max, delimiter = ',', fmt = '%.10f')
np.savetxt('snr.csv', snr, delimiter = ',', fmt = '%.10f')
np.savetxt('peak_dev.csv', peak_dev, delimiter = ',', fmt = '%.10f')
plt.imshow(snr)
plt.xlabel('Maximum time(ms)')
plt.ylabel('Bin number of histograms')
plt.colorbar()
ax = plt.gca()
ax.xaxis.get_major_formatter().set_powerlimits((0,1))
plt.savefig('snr')
