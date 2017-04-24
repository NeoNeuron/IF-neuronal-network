#!/usr/bin/python
import numpy as np 
import matplotlib.pyplot as plt 
import time
from scipy.optimize	import curve_fit
import random


def load_raster(loading_dir, filename, index): 
	"""
	convert raster data for certain neuron;
	Return an 1d array of spiking time points
	"""
	f = open(loading_dir + filename)
	counter = 0
	for line in f:
		if counter == index:
			raster = line
			break
		counter += 1
	raster = raster.replace('\n', '').split('\t')
	for s in raster:
		if s == '':
			raster.remove(s)
	raster = [float(s) for s in raster]
	raster = np.delete(raster, 0)
	return raster

def raster_to_bool(raster, dt):
	tmax = raster[-1]
	T = int(np.floor(tmax / dt)) + 1
	raster_bool = np.zeros(T)
	for num in raster:
		raster_bool[int(np.floor(num / dt))] = 1
	return raster_bool

def lfp(loading_dir, filename, index_list):
	post_current = open(loading_dir + filename)
	current = []
	for line in post_current:
		line = line.replace('\n', '').split('\t')
		line = [float(line[num]) for num in index_list]
		current.append(sum(line) / 3)
	post_current.close()
	return current

def load_matrix(loading_dir, filename):
	conMat = open(loading_dir + filename)
	matrix = []
	for line in conMat:
		line = line.replace('\n', '').split('\t')
		for s in line:
			if s == '':
				line.remove(s)
		line = [float(strs) for strs in line]
		matrix.append(line)
	conMat.close()
	matrix = np.array(matrix)
	return matrix

def isi_stat(raster, bins = 50, hist_plot = False):
	raster_temp = np.delete(np.insert(raster, 0, 0), -1)
	isi = raster - raster_temp
	[n, bins] = np.histogram(isi, bins = bins, normed = True)
	if hist_plot == True:
		plt.hist(isi, bins = bins, normed = 1)
		plt.show()
	return n, bins

def mean_rate(loading_dir, filename, index, tmax):
	""" 
	loading_dir: directory that neuronal data locate;
	index: index of pre-network neuron;
	tmax: maximum time considered in mean_rate() operation, unit in second; 
	"""
	# load raster data;
	raster = load_raster(loading_dir = loading_dir, filename = filename, index = index)
	raster = raster[raster < (tmax * 1000)]
	mean_rate = len(raster)*1.0 / tmax
	return mean_rate

def tdmi_parameters():
	"""
	Returns:
	sn_ratio: signal-noise ratio
	peak_time: time point of maximum mutual information
	time_constant: time constant for signal decay;
	"""
	# load signals
	signal_order = np.loadtxt("./tdmi/file-dat/tdmi_ordered.dat")
	signal_rand =  np.loadtxt("./tdmi/file-dat/tdmi_rand.dat")
	# calculate noise level
	noise_mean = signal_rand[:,1].mean()
	noise_std = signal_rand[:,1].std()
	# allocate maximum mutual information signal
	signal_max_ind = np.argmax(signal_order[:,1])
	signal_max = signal_order[signal_max_ind,1]
	peak_time = signal_order[signal_max_ind,0]

	# calculate the signal-noise ratio in TDMI data;
	sn_ratio = signal_max / noise_mean
	# find the duriation of sufficient signal
	#if type(signal_max_ind) == int:
	ind = signal_max_ind
	while signal_order[ind, 1] >= noise_mean + noise_std:
		ind -= 1
		if ind == 0:
			break
	signal_front_ind = ind
	ind = signal_max_ind
	while signal_order[ind, 1] >= noise_mean + noise_std:
		ind += 1
		if ind == len(signal_order): 
			break
	signal_back_ind = ind
	# Judge whether fit or not
	if  signal_back_ind - signal_max_ind > 3:
		fexp = lambda x, a, b, c: a*np.exp(-b * x) + c
		try:
			popt, pcov = curve_fit(fexp, signal_order[signal_max_ind:signal_back_ind, 0], signal_order[signal_max_ind:signal_back_ind, 1], p0 = (1, 1, 0))
		except RuntimeError:
			popt = [0, 0, 0]
		# decayed time constant
		if popt[1] == 0:
			time_constant = 0
		else:
			time_constant = 1 / popt[1]
	else:
		peak_time = -1
		time_constant = -1
	return sn_ratio, peak_time, time_constant

def MakeTitle(saving_filename):
	# split filename into subunit and reassumble them into acceptable title string;
	# sub_unit[0] = 'tdmi'
	# sub_unit[1] = index of neuron
	# sub_unit[2] = connecting order
	# sub_unit[3] = time range
	# sub_unit[4] = str labels of classification
	# sub_unit[5] = expected occupancy for historgram in mutual information calculation
	# sub_unit[6] = timing step for MI
	# sub_unit[7] = maximum negative time delay
	# sub_unit[8] = maximum positive time delay
	sub_unit = saving_filename.split('-')
	title = sub_unit[0].upper()
	title += ' #' + sub_unit[1]
	if sub_unit[2] == '1':
		title += ' $1^{st}$ '
	elif sub_unit[2] == '2':
		title += ' $2^{nd}$ '
	sub_sub_unit = sub_unit[3].split('_')
	title += sub_sub_unit[0] + '~' + sub_sub_unit[1] + ' ms '
	title += sub_unit[4] + '\n'
	title += 'expected occupancy = ' + sub_unit[5]
	sub_sub_unit = sub_unit[6].split('_')
	title += ' dt = ' + sub_sub_unit[0] + '.' + sub_sub_unit[1] + ' ms'
	title += ' NTD = ' + sub_unit[7]
	title += ' PTD = ' + sub_unit[8]
	return title

def DivideNeuronalFunction(neuron_types, neuron_list):
	counter = 0
	for i in neuron_list:
		if neuron_types[i] == 1:
			counter += 1 
	neuron_numbers = [counter, len(neuron_list) - counter]
	return neuron_numbers

# Create text in TDMI plot, including type of target neuron and the number of neuron it connected as well as their type;
def CreateText(loading_dir, neuron_index, order, classification, num):
	# loading files;
	pre_net_types = np.genfromtxt(loading_dir + 'preNeuron.txt', dtype = int, usecols = 0)
	post_net_types = np.genfromtxt(loading_dir + 'postNeuron.txt', dtype = int, usecols = 0)
	con_mat = np.genfromtxt(loading_dir + 'conMat.txt', dtype = int)
	# create text variable
	text = '#' + str(neuron_index) + ' neuron '
	if pre_net_types[neuron_index] == 1:
		text += 'is excitatory\n'
	else:
		text += 'is inhibitory\n'
	# consider the order of connection
	neuron_list_1 = [i for i in range(np.size(con_mat[neuron_index, :])) if con_mat[neuron_index, i] == 1]
	if order == 1:
		neuron_list_all = neuron_list_1
	else:
		neuron_list_2 = []
		for ind in neuron_list_1:
			neuron_list_2 += [j for j in range(np.size(con_mat[neuron_index, :])) if con_mat[ind, j] == 1]
		del ind
		neuron_list_2 = np.unique(neuron_list_2)
		neuron_list_all = np.setdiff1d(neuron_list_2, neuron_list_1)
	# consider the classification of connection
	if classification == 'exc':
		neuron_list = [ind for ind in neuron_list_all if post_net_types[ind] == 1]
	elif classification == 'inh':
		neuron_list = [ind for ind in neuron_list_all if post_net_types[ind] == 0]
	else:
		neuron_list = neuron_list_all
	# select num neurons from neuron pool above
	if np.size(neuron_list) > num and num > 0:
		neuron_list = random.sample(neuron_list, num)
	else:
		pass	

	text += 'LFP is generated by ' + str(np.size(neuron_list)) + ' neurons\n'
	# classify neuronal type of all this neurons
	neuron_numbers = DivideNeuronalFunction(post_net_types, neuron_list = neuron_list)
	text += 'including ' + str(neuron_numbers[0]) + ' excitatroy and ' + str(neuron_numbers[1]) + ' inhibitory neurons'
	return text

def PlotTdmi(saving_filename, figure_text):
	data_ordered = np.loadtxt("./tdmi/file-dat/tdmi_ordered.dat")
	data_rand = np.loadtxt("./tdmi/file-dat/tdmi_rand.dat")
	# basic plot
	fig = plt.figure(0, figsize=(10,8), dpi=60)
	plt.plot(data_ordered[:,0], data_ordered[:,1], label = "tdmi-original")
	plt.plot(data_rand[:,0], data_rand[:,1], label = "tdmi-swapped")
	# setting axis range;
	x_max = np.max(data_ordered[:,0])
	x_min = np.min(data_ordered[:,0])
	if np.max(data_ordered[:,1]) > np.max(data_rand[:,1]):
		y_max = np.max(data_ordered[:,1])
	else:
		y_max = np.max(data_rand[:,1])
	if np.min(data_ordered[:,1] < np.min(data_rand[:,1])):
		y_min = np.min(data_ordered[:,1])
	else:
		y_min = np.min(data_rand[:,1])
	abs_diff = y_max - y_min;
	y_min -= abs_diff * 0.1
	y_max += abs_diff * 0.1
	plt.axis([x_min, x_max, y_min, y_max])
	# setting labels and title
	plt.xlabel("Time-delay(ms)")
	plt.ylabel("Mutual Information(bits)")
	title = MakeTitle(saving_filename = saving_filename)
	# title = saving_filename
	plt.title(title)
	plt.legend()
	plt.grid(True)
	# add text
	x_text_pos = (x_max - x_min) * 0.05 + x_min
	y_text_pos = (y_max - y_min) * 0.9 + y_min
	plt.text(x_text_pos, y_text_pos, figure_text, weight = 'light')
	saving_dir = "./tdmi/figure-eps/"
	plt.savefig(saving_dir + saving_filename + '.eps')
	plt.close(0)
	del data_rand
	del data_ordered

# start = time.clock()
# load current data:
# loading_dir = '/media/kyle/Drive/ResearchData/Apr18/test1/'
# aaa = load_raster(loading_dir = loading_dir, index = 10)
# [b, c] = isi_stat(aaa, 10, True)


# finish = time.clock()

# print finish - start
