import numpy as np
from matplotlib import pyplot as plt

#Energie plot
x, E = np.loadtxt("energy_data.csv", float, delimiter=',', unpack=True)

plt.plot(x,E)
plt.show()

#Wellenfkt für n = 0
x, psi_0 = np.loadtxt("psi_0_data.csv", float, delimiter=',', unpack = True)

plt.plot(x,psi_0)
plt.show()

x, psi_1 = np.loadtxt("psi_1_data.csv", float, delimiter=',', unpack = True)

plt.plot(x,psi_1)
plt.show()