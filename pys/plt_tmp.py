import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('data/mi/mi_bd.csv')
plt.plot(data['timelag'], data['mi'])
plt.xlabel('timelag')
plt.ylabel('mutual info')
plt.show()
