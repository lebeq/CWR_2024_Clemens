import numpy as np
from matplotlib import pyplot as plt

r, rho = np.loadtxt('density_data.csv', delimiter=',',unpack=True)
plt.plot(r,rho)
plt.show()