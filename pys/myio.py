#!/usr/bin/python
import numpy as np
import struct as st
import datetime
import random

def read_1d(filename, dtype):
  f = open(filename, 'rb')
  data_len = st.unpack('Q',f.read(16))[1]
  data1d = np.emtpy(data_len)
  for i in range(data_len):
    data1d[i] = st.unpack(dtype, f.read(8))[0]
  f.close()
  return data1d

def read_2d(filename, dtype):
  f = open(filename, 'rb')
  data_shape = st.unpack('Q',f.read(16))
  data2d = np.emtpy(data_shape)
  for i in range(data_shape[0]):
    for j in range(data_shape[1]):
      data2d[i][j] = st.unpack(dtype, f.read(8))[0]
  f.close()
  return data1d

def uni_num():
  nowTime=datetime.datetime.now().strftime("%Y%m%d%H%M%S") 
  randomNum=random.randint(0,1000) 
  if randomNum < 10:
    randomNum = 2*str(0) + str(randomNum)
  elif randomNum >= 10 and randomNum < 100:
    randomNum = str(0) + str(randomNum)
  uniqueNum = str(nowTime) + str(randomNum)
  return uniqueNum
