import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn

loading_dir = './results/Apr22/t3mix-55/'
all_data = pd.read_csv(loading_dir + 'asc.csv', sep = '\t')
# data strucure:{'index', 'connections', 'connecting portion', 'mean firing rate', 'signal noise ratio', 'peak time', 'time constant'}
init_g = False
init_r = False
plt.figure(1)

for i in range(len(all_data)):
	if all_data.loc[i,'time constant'] > 0:
		if init_g == False:
			plt.scatter(all_data.loc[i,'connections'],all_data.loc[i,'signal noise ratio'], c='g', alpha = 0.2, label = 'effective event')
			init_g = True
		else:
			plt.scatter(all_data.loc[i,'connections'],all_data.loc[i,'signal noise ratio'], c='g', alpha = 0.2)
	else:
		if init_r == False:
			plt.scatter(all_data.loc[i,'connections'],all_data.loc[i,'signal noise ratio'], c='r',alpha = 0.8, label = 'ineffective event')
			init_r = True
		else:
			plt.scatter(all_data.loc[i,'connections'],all_data.loc[i,'signal noise ratio'], c='r',alpha = 0.8)

plt.xticks(range(0,12,1))
plt.xlabel('Number of Neurons', fontsize = 15)
plt.ylabel('Signal-noise Ratio', fontsize = 15)
plt.legend(loc = 2)
plt.savefig(loading_dir + 'snr-num.png')

init_g = False
init_r = False
plt.figure(2)
for i in range(len(all_data)):
	if all_data.loc[i,'time constant'] > 0:
		if init_g == False:
			plt.scatter(all_data.loc[i, 'mean firing rate'], all_data.loc[i, 'signal noise ratio'], c='g', alpha = 0.4, label = 'effective event')
			init_g = True
		else:
			plt.scatter(all_data.loc[i, 'mean firing rate'], all_data.loc[i, 'signal noise ratio'], c='g', alpha = 0.4)
	else:
		if init_r == False:
			plt.scatter(all_data.loc[i, 'mean firing rate'], all_data.loc[i, 'signal noise ratio'], c='r',alpha = 0.8, label = 'ineffective event')
			init_r = True
		else:
			plt.scatter(all_data.loc[i, 'mean firing rate'], all_data.loc[i, 'signal noise ratio'], c='r',alpha = 0.8)
plt.xlabel('Mean firing rate (Hz)', fontsize = 15)
plt.ylabel('Signal-noise Ratio', fontsize = 15)
plt.legend()
plt.savefig(loading_dir + 'snr-rate.png')

plt.figure(3)
plt.scatter(all_data['mean firing rate'], all_data['signal noise ratio'], c = all_data['connecting portion'], cmap = cm.jet, alpha = 0.5)

plt.colorbar()
plt.xlabel('Mean firing rate (Hz)', fontsize = 15)
plt.ylabel('Signal-noise Ratio', fontsize = 15)
plt.savefig(loading_dir + 'snr-rate-portion.png')

plt.figure(4)
plt.scatter(all_data['mean firing rate'], all_data['signal noise ratio'], c = all_data['connections'], cmap = cm.jet, alpha = 0.6)

plt.colorbar()
plt.xlabel('Mean firing rate (Hz)', fontsize = 15)
plt.ylabel('Signal-noise Ratio', fontsize = 15)
plt.savefig(loading_dir + 'snr-rate-num.png')

plt.show()
