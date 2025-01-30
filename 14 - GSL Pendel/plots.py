import numpy as np
from matplotlib import pyplot as plt

data_pos = np.loadtxt("pos_file.csv", float, delimiter = ',')
data_energy = np.loadtxt("energy_file.csv", float, delimiter=',')

time = data_pos[:,0]
positions = data_pos[:,1:]
energy = data_energy[:,1]

#plot positions
plt.plot(time,positions, color='k')
plt.xlabel("t [s]")
plt.ylabel("position")
plt.show()

#plot energy
plt.plot(time,energy)
plt.xlabel("t [s]")
plt.ylabel("E_total")
plt.show()
