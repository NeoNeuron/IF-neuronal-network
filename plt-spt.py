import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn

num = 10
loading_dir = './tdmi/figure-eps/Apr24/t6-91/'
all_data = []

for i in range(num):
	data = np.genfromtxt(loading_dir + 'tdmi-' + str(i+1) + '/ana1.csv', skip_header = 1, usecols = [2, 4])
	# print data[0][:]
	all_data.append(data)
init_g = False
init_r = False

for i in range(num - 1):
	for j in range(len(all_data[i][:])):
		if all_data[i][j][1] > 0:
			if init_g == False:
				plt.scatter(i+1,all_data[i][j][0],c='g', alpha = 0.2, label = 'effective event')
				init_g = True
			else:
				plt.scatter(i+1,all_data[i][j][0],c='g', alpha = 0.2)
		else:
			if init_r == False:
				plt.scatter(i+1,all_data[i][j][0],c='r',alpha = 0.8, label = 'ineffective event')
				init_r = True
			else:
				plt.scatter(i+1,all_data[i][j][0],c='r',alpha = 0.8)

if all_data[-1][1] > 0:
	plt.scatter(num,all_data[-1][0],c='g', alpha = 0.8)
else:
	plt.scatter(num,all_data[-1][0],c='r',alpha = 0.8) 

all_data_box = [line[:,0].copy() for line in all_data[:-1]]
plt.boxplot(all_data_box)
plt.xlim([0.5,0.5 + num])
plt.xlabel('neuronal number', fontsize = 16)
plt.ylabel('signal-noise ratio', fontsize = 16)
plt.legend(loc = 2)
plt.savefig(loading_dir + 'snr-boxplot.png')
plt.show() 