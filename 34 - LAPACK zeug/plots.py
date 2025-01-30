import numpy as np
from matplotlib import pyplot as plt

dim, time = np.loadtxt("time_data.csv", delimiter=',', unpack=True)
plt.scatter(dim,time)
plt.show()