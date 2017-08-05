#!/usr/bin/python
#coding:utf-8
import numpy as np
import matplotlib.pyplot as plt
import time
# from scipy.optimize	import curve_fit
import random
import os
import pandas as pd

def load_spike_train(path, index):
	"""
	Load spike train of certain neuron, labeled by index;
	path: path of data file storing raster data
	index: neuronal index of those whoes spike train are in raster data files
	Return:
		spikes: 1 d array of spiking time points;
	"""
	spikes = np.genfromtxt(path, delimiter = ',', skip_header = index, max_rows = 1)
	spikes =  np.array([i for i in spikes if str(i) != 'nan'])
	return spikes

def raster_to_bool(raster, dt):
	"""
	Convert spike train to binary time series;
	raster: original spike train, which is a float array;
	dt: width of time period of each element in the binary element;
	Return:
		raster_bool: 1d bool array;
	"""
	tmax = raster[-1]
	T = int(np.floor(tmax / dt)) + 1
	raster_bool = np.zeros(T)
	for num in raster:
		raster_bool[int(np.floor(num / dt))] = 1
	return raster_bool

def lfp(path, index_list):
	"""
	Load local field potential from data files of currents
	path: path of current files
	index_list: list of targeting neuronal index;
	Return:
		current: 1d float array, local field potential;
	"""
	currents = np.genfromtxt(path, delimiter = ',', usecols = index_list)
	current = np.mean(currents, axis = 1)
	return current

def load_matrix(path):
	"""
	Return:
	matrix: an ndarray-typed data;
	"""
	conMat = open(path)
	matrix = []
	for line in conMat:
		line = line.replace('\n', '').split(',')
		for s in line:
			if s == '':
				line.remove(s)
		line = [float(strs) for strs in line]
		matrix.append(line)
	conMat.close()
	matrix = np.array(matrix)
	return matrix

def isi_stat(raster, bins = 50, hist_plot = False):
	"""
	raster: double vectors
	bin: number of bins in histogram
	hist_plot: whether show the histogram, true for plot, false for not;
	return:
		n: frequency of each bin
		bins: frontier edge of bins
	"""
	raster_temp = np.delete(np.insert(raster, 0, 0), -1)
	isi = raster - raster_temp
	n, bins = np.histogram(isi, bins = bins)
	if hist_plot == True:
		plt.hist(isi, bins = bins)
		plt.show()
	return n, bins

def mean_rate(path, index, tmax):
	"""
	loading_dir: directory that neuronal data locate;
	index: index of pre-network neuron;
	tmax: maximum time considered in mean_rate() operation, unit in second;
	"""
	# load raster data;
	raster = load_spike_train(path = path, index = index)
	raster = raster[raster < (tmax * 1000)]
	mean_rate = len(raster)*1.0 / tmax
	return mean_rate

def tdmi_parameters(mi):
	"""
	mi: DataFame containing 'timelag', 'ordered' and 'random'
	Returns:
	snr: signal-noise ratio
	peak_position: postition of mutual information peak, True for postitive and False for negative;
	"""
	time_series = mi['timelag']
	signal_ordered = mi['ordered']
	signal_random = mi['random']
	# calculate noise level
	noise_mean = signal_random.mean()
	noise_std = signal_random.std()
	# allocate maximum mutual information signal
	signal_max_ind = signal_ordered.argmax()
	signal_max = signal_ordered[signal_max_ind]
	peak_time = time_series[signal_max_ind]
	if peak_time > 0:
		peak_position = True
	else:
		peak_position = False

	if noise_mean != 0:
		# calculate the signal-noise ratio in TDMI data;
		snr = signal_max / noise_mean
	else:
		snr = 'inf'
		peak_position = 'nan'
	return snr, peak_position

def MakeTitle(saving_filename):
	"""
	Split filename into subunit and reassumble them into acceptable title string;
	unit[0] = 'tdmi'
	unit[1] = index of neuron
	unit[2] = connecting order
	unit[3] = lower bond of time interval
	unit[4] = upper bond of time interval
	unit[5] = expected occupancy for historgram in mutual information calculation
	unit[6] = timing step for TDMI
	unit[7] = maximum negative time delay
	unit[8] = maximum positive time delay
	"""
	unit = saving_filename.split('-')
	title = unit[0].upper()
	title += ' #' + unit[1]
	if unit[2] == '1':
		title += ' $1^{st}$ '
	elif unit[2] == '2':
		title += ' $2^{nd}$ '
	title += unit[3] + '~' + unit[4] + ' ms\n'
	title += 'Expected Occupancy = ' + unit[5]
	title += ' timing step = ' + unit[6][0] + '.' + unit[6][1:] + ' ms'
	title += ' NTS = ' + unit[7]
	title += ' PTS = ' + unit[8]
	return title

def DivideNeuronalTypes(neuron_types, neuron_list):
	counter = 0
	for i in neuron_list:
		if neuron_types[i] == 1:
			counter += 1
	num_exc = counter
	num_inh = len(neuron_list) - counter
	return num_exc, num_inh

def CreateText(loading_dir, neuron_index, connecting_list):
	"""
	Create text in TDMI plot, including type of target neuron and the number of neuron it connected as well as their type;
	"""
	# loading files;
	prenet_type = np.genfromtxt(loading_dir + 'preNeuron.txt', dtype = int, usecols = 0, delimiter = ',')
	postnet_type = np.genfromtxt(loading_dir + 'postNeuron.txt', dtype = int, usecols = 0, delimiter = ',')
	conMat = load_matrix(loading_dir = loading_dir, filename = 'conMat.txt')
	# create text variable
	text = '#' + str(neuron_index) + ' neuron '
	if prenet_type[neuron_index] == 1:
		text += 'is excitatory\n'
	else:
		text += 'is inhibitory\n'
	num_exc, num_inh = DivideNeuronalTypes(neuron_types = postnet_type, neuron_list = connecting_list)

	text += 'LFP is generated by ' + str(np.size(connecting_list)) + ' neurons\n' + 'including ' + str(num_exc) + ' excitatroy and ' + str(num_inh) + ' inhibitory neurons'
	return text

def PlotTdmi(time_series, signal_order, signal_rand, saving_dir, saving_filename, figure_text = None):
	# basic plot
	fig = plt.figure(0, figsize=(10,8), dpi=60)
	plt.plot(time_series, signal_order, label = "TDMI-original")
	plt.plot(time_series, signal_rand, label = "TDMI-swapped")
	# setting axis range;
	x_max = time_series[len(time_series) - 1]
	x_min = time_series[0]
	if max(signal_order) > max(signal_rand):
		y_max = max(signal_order)
	else:
		y_max = max(signal_rand)
	if min(signal_order) < min(signal_rand):
		y_min = min(signal_order)
	else:
		y_min = min(signal_rand)
	abs_diff = y_max - y_min;
	y_min -= abs_diff * 0.1
	y_max += abs_diff * 0.1
	plt.axis([x_min, x_max, y_min, y_max])
	# setting labels and title
	plt.xlabel("Time-delay(ms)")
	plt.ylabel("Mutual Information(bits)")
	# title = MakeTitle(saving_filename = saving_filename)
	title = saving_filename
	plt.title(title)
	plt.legend()
	plt.grid(True)
	# add text
	x_text_pos = (x_max - x_min) * 0.05 + x_min
	y_text_pos = (y_max - y_min) * 0.9 + y_min
	if figure_text != None:
		plt.text(x_text_pos, y_text_pos, figure_text, weight = 'light')
	plt.savefig(saving_dir + saving_filename + '.eps')
	plt.close(0)

def plot(filename):
	sta = pd.read_csv('data/sta.csv')
	mi = pd.read_csv('data/mi/mi_sl.csv')
	lcc = pd.read_csv('data/lcc/lcc.csv')
	plt.figure(figsize = (15,5), dpi = 75)
	plt.subplot(1,3,1)
	plt.plot(sta['timelag'], sta['sta'])
	plt.xlabel('time-delay(ms)')
	plt.ylabel('LFP')
	plt.grid(True)
	plt.subplot(1,3,2)
	plt.plot(lcc['timelag'], lcc['lcc'])
	plt.xlabel('time-delay(ms)')
	plt.ylabel('Pearson\' Correlation')
	plt.grid(True)
	plt.subplot(1,3,3)
	plt.plot(mi['timelag'], mi['ordered'])
	plt.plot(mi['timelag'], mi['random'])
	plt.xlabel('time-delay(ms)')
	plt.ylabel('MI(bits)')
	plt.grid(True)
	plt.savefig('./data/figure/' + filename)
	# plt.show()
