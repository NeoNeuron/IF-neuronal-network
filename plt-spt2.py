import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import seaborn

num = 10
loading_dir = './tdmi/figure-eps/Apr24/t6-91/'
all_data = pd.DataFrame({'connecting portion':[], 'mean firing rate':[], 'signal noise ratio':[], 'peak time':[], 'time constant':[], 'neuron number':[]})

for i in range(num):
	data = pd.read_csv(loading_dir + 'tdmi-' + str(i+1) + '/ana1.csv', sep = '\t')
	slength = len(data['connecting portion'])
	data['neuron number'] = pd.Series(np.ones(slength).astype(int) * (i + 1), index = data.index)
	all_data = all_data.append(data, ignore_index = True)

plt.figure(1)
plt.scatter(all_data['mean firing rate'], all_data['signal noise ratio'], c = all_data['connecting portion'], cmap = cm.jet, alpha = 0.6)

plt.colorbar()
plt.xlabel('Mean firing rate (Hz)', fontsize = 16)
plt.ylabel('signal-noise ratio', fontsize = 16)
plt.savefig(loading_dir + 'mrate-snr1.png')

plt.figure(2)
plt.scatter(all_data['mean firing rate'], all_data['signal noise ratio'], c = all_data['neuron number'], cmap = cm.jet, alpha = 0.6)

plt.colorbar()
plt.xlabel('Mean firing rate (Hz)', fontsize = 16)
plt.ylabel('signal-noise ratio', fontsize = 16)
plt.savefig(loading_dir + 'mrate-snr2.png')

plt.show() 