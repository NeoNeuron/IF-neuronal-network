#!/usr/bin/python
import numpy as np 
import matplotlib.pyplot as plt 
import time

def load_raster(loading_dir, filename, index): # convert raster data to binary trains;
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



start = time.clock()
# load current data:
loading_dir = '/media/kyle/Drive/ResearchData/Apr18/test1/'
# aaa = load_raster(loading_dir = loading_dir, index = 10)
# [b, c] = isi_stat(aaa, 10, True)


finish = time.clock()

print finish - start
