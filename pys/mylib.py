#!/usr/bin/python
#coding:utf-8
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize	import curve_fit
import random
import os
import pandas as pd


def Compile():
	"""
	compile *.cpp files
	"""
	os.system("g++ ./multi-network/*.cpp ./io/io.cpp -o ./multi-network/two-network-system.out")
	os.system("g++ ./lfp/*.cpp ./io/io.cpp -o ./lfp/calculate-lfp.out")
	os.system("g++ ./tdmi/*.cpp ./io/io.cpp -o ./tdmi/calculate-tdmi.out")


def load_raster(loading_dir, filename, index):
	"""
	convert raster data for certain neuron;
	Return an 1 d array of spiking time points
	"""
	f = open(loading_dir + filename)
	counter = 0
	for line in f:
		if counter == index:
			raster = line
			break
		counter += 1
	f.close()
	raster = raster.replace('\n', '').split(',')
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
		line = line.replace('\n', '').split(',')
		line = [float(line[num]) for num in index_list]
		current.append(sum(line) / 3)
	post_current.close()
	return current

def load_matrix(loading_dir, filename):
	"""
	Return:
	matrix: an ndarray-typed data;
	"""
	conMat = open(loading_dir + filename)
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
	raster_temp = np.delete(np.insert(raster, 0, 0), -1)
	isi = raster - raster_temp
	n, bins = np.histogram(isi, bins = bins)
	if hist_plot == True:
		plt.hist(isi, bins = bins)
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

def tdmi_parameters(mi):
	"""
	Returns:
	sn_ratio: signal-noise ratio
	peak_time: time point of maximum mutual information
	time_constant: time constant for signal decay;
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

	if noise_mean != 0:
		# calculate the signal-noise ratio in TDMI data;
		sn_ratio = signal_max / noise_mean
		# find the duriation of sufficient signal
		# if type(signal_max_ind) == int:
		ind = signal_max_ind
		if ind >= noise_mean + noise_std:
			while signal_ordered[ind] >= noise_mean + noise_std:
				if ind == len(signal_ordered) - 1:
					break
				ind += 1
			signal_back_ind = ind
			if  signal_back_ind - signal_max_ind > 3:
				fexp = lambda x, a, b, c: a*np.exp(-b * x) + c
				try:
					popt, pcov = curve_fit(fexp, time_series.ix[signal_max_ind:signal_back_ind].values, signal_ordered.ix[signal_max_ind:signal_back_ind].values, p0 = (1, 1, 0))
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
		else:
			peak_time = -1
			time_constant = -1
	else:
		sn_ratio = 0
		peak_time = -1
		time_constant = -1
	return sn_ratio, peak_time, time_constant

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

# start = time.clock()
# load current data:
# loading_dir = '/media/kyle/Drive/ResearchData/Apr18/test1/'
# aaa = load_raster(loading_dir = loading_dir, index = 10)
# [b, c] = isi_stat(aaa, 10, True)


# finish = time.clock()

# print finish - start
