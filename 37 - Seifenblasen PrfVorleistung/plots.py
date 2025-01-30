import numpy as np
from matplotlib import pyplot as plt

d, vol, var = np.loadtxt("volume_data.csv", delimiter=',',unpack=True)
d1, mean = np.loadtxt("mean_data.csv", delimiter=',', unpack=True)
d2, ana_sol = np.loadtxt("analytic_data.csv", delimiter=',', unpack=True)

plt.errorbar(d,vol,yerr=var, alpha = 0.5, fmt = 'o')
plt.plot(d2, ana_sol, 'g-')
plt.show()

plt.scatter(d1,mean,alpha=0.5)
plt.show()
