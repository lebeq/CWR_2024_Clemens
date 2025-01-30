import numpy as np
from matplotlib import pyplot as plt

t, infected = np.loadtxt('infect_data.csv', float, delimiter=',', unpack=True)
plt.scatter(t,infected)
plt.show()