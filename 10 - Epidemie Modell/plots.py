import numpy as np
from matplotlib import pyplot as plt

#auslesen von data
t, S, E, I, R = np.loadtxt("SEIRdata.csv", float, delimiter=",", unpack=True)

plt.plot(t,S, 'b-', label = 'Susceptible')
plt.plot(t,E, 'y--', label = 'Exposed')
plt.plot(t,I, 'r--', label = 'Infected')
plt.plot(t,R, 'g--', label = 'Removed')
plt.xlabel('Time')
plt.legend()
plt.show()