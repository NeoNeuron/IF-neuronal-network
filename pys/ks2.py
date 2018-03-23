#!/usr/bin/python
# Peakcock test, two-dimensional KS test for two emperical data;
import numpy as np
import matplotlib.pyplot as plt

def cdf2d(H, xedges, yedges, orient = (False, False)):
    # Prerequirest: H is normed;
    cdf = np.empty(H.shape)
    xstep = xedges[1] - xedges[0]
    ystep = yedges[1] - yedges[0]
    A = xstep * ystep
    # orientation selection:
    if orient[0]:
        H = np.flipud(H)
    if orient[1]:
        H = np.fliplr(H)
    # Start Calculation;
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            if j == 0:
                cdf[i][0] = H[i][0] * A
            else:
                cdf[i][j] = cdf[i][j - 1] + H[i][j] * A
        if i != 0:
            cdf[i] += cdf[i - 1]
    return cdf

def Peacock2d(X, Y, bin_num_1, bin_num_2):
    # Find the range of histograms;
    max1 = X[0].max() if X[0].max() >= Y[0].max() else Y[0].max()
    min1 = X[0].min() if X[0].min() <= Y[0].min() else Y[0].min()
    max2 = X[1].max() if X[1].max() >= Y[1].max() else Y[1].max()
    min2 = X[1].min() if X[1].min() <= Y[1].min() else Y[1].min()
    # Calculate histograms;
    H1, edges1, edges2 = np.histogram2d(X[0], X[1], bins = [bin_num_1, bin_num_2], range = [[min1, max1], [min2, max2]], normed = True)
    H2, edges1, edges2 = np.histogram2d(Y[0], Y[1], bins = [bin_num_1, bin_num_2], range = [[min1, max1], [min2, max2]], normed = True)
    C1ff = cdf2d(H1, edges1, edges2, orient = (False, False))
    C1ft = cdf2d(H1, edges1, edges2, orient = (False, True))
    C1tf = cdf2d(H1, edges1, edges2, orient = (True, False))
    C1tt = cdf2d(H1, edges1, edges2, orient = (True, True))

    C2ff = cdf2d(H2, edges1, edges2, orient = (False, False))
    C2ft = cdf2d(H2, edges1, edges2, orient = (False, True))
    C2tf = cdf2d(H2, edges1, edges2, orient = (True, False))
    C2tt = cdf2d(H2, edges1, edges2, orient = (True, True))

    diff_ff = abs(C1ff - C2ff).max()
    diff_tf = abs(C1tf - C2tf).max()
    diff_ft = abs(C1ft - C2ft).max()
    diff_tt = abs(C1tt - C2tt).max()

    diff_max = max(diff_ff, diff_ft, diff_tf, diff_tt)
    x_len = X.shape[1]
    y_len = Y.shape[1]
    p = np.exp(-2*x_len*y_len/(x_len + y_len)*diff_max**2)
    return diff_max, p

num_rand = 2
num_trial = 1e5*2**10
x = np.random.normal(size = num_trial)
y1 = np.random.normal(size = num_trial)
ay = 0.5
k = 1
y2 = ay*y1 + k*x + np.random.normal(size = num_trial)

# Settings for histograms;
bin_num_1 = 10
bin_num_2 = 10

# import original data; Suppose D are 2 by N arrays;
D = np.concatenate(([x], [y2]))
D_val = np.zeros(11)
p_val = np.zeros(11)
for i in range(11):
  length = int(1e5*2**i )
  D_tmp = D[:,:length]
  D_val[i], p_val[i] = Peacock2d(D_tmp, D, bin_num_1, bin_num_2)

np.savetxt('./data/tmp/D_val.csv', D_val, fmt = '%.18f')
np.savetxt('./data/tmp/p_val.csv', p_val, fmt = '%.18f')
print('D value:\n')
print D_val
print('p-value:\n')
print p_val
