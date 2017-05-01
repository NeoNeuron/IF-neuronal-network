import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn

num = 10
loading_dir = './tdmi/figure-eps/Apr24/t6-91/'
all_data = pd.DataFrame({'mean firing rate':[], 'signal noise ratio':[], 'peak time':[], 'time constant':[]})

for i in range(num):
	data = pd.read_csv(loading_dir + 'tdmi-' + str(i+1) + '/ana1.csv', sep = '\t')
	# print data[0][:]
	all_data = all_data.append(data, ignore_index = True)

init_g = False
init_r = False

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

plt.xlabel('Mean firing rate (Hz)', fontsize = 16)
plt.ylabel('signal-noise ratio', fontsize = 16)
plt.legend(loc = 2)
plt.savefig(loading_dir + 'mrate-snr.png')
plt.show() 