#!/usr/bin/python
import numpy as np
import pandas as pd
import mylib
import matplotlib.pyplot as plt

# By graph;
# mi = pd.read_csv('./data/mi/mi_sl.csv')
# peak_value = mi['ordered'].max()
# print peak_value
# # pmi_ordered = hist(mi['ordered'])
# # pmi_random = hist(mi['random'])
# num = len(mi['ordered'])
# plt.hist(mi['ordered'], bins = 50, label = 'random', weights = 1.0/num*np.ones(num))
# plt.axvline(x=peak_value, ymin=0, ymax=1, linewidth=2, color = 'r')
# # plt.hist(mi['random'], bins = 50, label = 'ordered', normed = True)
# # plt.plot(mi['ordered'])
# plt.show()

# By numbers;
mi = pd.read_csv('./data/mi/mi_sl.csv')
peak_value = mi['ordered'].max()
[noise_mean, noise_std] = mylib.baseline(mi)
deviate_level = (peak_value - noise_mean) / noise_std
print deviate_level
