# this program aims to analyze multiple tdmi signal from multiple trials;
import numpy as np
from scipy.optimize	import curve_fit 

def tdmi_parameters():
	# load signals
	signal_order = np.loadtxt("./file-dat/tdmi_ordered.dat")
	signal_rand =  np.loadtxt("./file-dat/tdmi_rand.dat")
	# calculate noise level
	noise_mean = np.mean(signal_rand[:,1])
	noise_std = np.std(signal_rand[:,1])
	# allocate maximum mutual information signal
	signal_max_ind = np.argmax(signal_order[:,1])
	signal_max = signal_order[signal_max_ind,1]
	signal_max_time = signal_order[signal_max_ind,0]
	# calculate the signal-noise ratio in TDMI data;
	sn_ratio = signal_max / noise_mean
	# find the duriation of sufficient signal
	#if type(signal_max_ind) == int:
	ind = signal_max_ind
	while signal_order[ind, 1] >= noise_mean + noise_std:
		ind -= 1
	signal_front_ind = ind
	ind = signal_max_ind
	while signal_order[ind, 1] >= noise_mean + noise_std:
		ind += 1
	signal_back_ind = ind
	# elif type(signal_max_ind) == list:
	# 	ind = signal_max_ind[0]
	# 	while signal_order[ind, 1] >= noise_mean + noise_std:
	# 		ind -= 1
	# 	signal_front_ind = ind
	# 	ind = signal_max_ind[-1]
	# 	while signal_order[ind, 1] >= noise_mean + noise_std:
	# 		ind += 1
	# 	signal_back_ind = ind
	# calculate time constant for the decay;

	def exp_template(x, a, b, c):
		return a*np.exp(-b * x) + c

	popt, pcov = curve_fit(exp_template, signal_order[signal_max_ind:signal_back_ind, 0], signal_order[signal_max_ind:signal_back_ind, 1], p0 = (1, 1, 0))

	decay_tau = 1 / popt[1]

	# print noise_mean
	# print '\n'
	# print noise_std
	# print '\n'
	# print signal_max_time
	# print '\n'
	# print sn_ratio
	# print '\n'
	# print decay_tau
	# print '\n'
	return sn_ratio, signal_max_time, decay_tau

a, b, c = tdmi_parameters()

print a
print b
print c
