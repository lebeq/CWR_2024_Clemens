import numpy as np
from matplotlib import pyplot as plt

x1, psi_1 = np.loadtxt("psi_1_data.csv",delimiter=',',unpack=True)
x2, psi_2 = np.loadtxt("psi_2_data.csv",delimiter=',',unpack=True)
x3, psi_3 = np.loadtxt("psi_3_data.csv",delimiter=',',unpack=True)

plt.plot(x1,psi_1,'r',label='n = 1')
plt.plot(x2,psi_2,'g',label='n = 2')
plt.plot(x3,psi_3,'k',label='n = 3')
plt.xlabel('pos')
plt.ylabel('Wellenfunktion')
plt.legend()
plt.savefig('./wellenfktionen.png')
plt.show()
