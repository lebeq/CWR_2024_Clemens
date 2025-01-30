import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import *

data_pos = np.loadtxt("pos_file.csv", float, delimiter = ',')
data_energy = np.loadtxt("energy_file.csv", float, delimiter=',')
residuums = np.loadtxt("residuum_file.csv", float, delimiter=',')

#time = data_pos[:,0]
#positions = data_pos[:,1]
#energy = data_energy[:,1]
dt = residuums[:,0]
res = residuums[:,1]

def linfunc(x,m,b):
    return m*x+b

popt, pcov = curve_fit(linfunc, np.log(dt), np.log(res))

plt.plot(np.log(dt), linfunc(np.log(dt),*popt))
plt.suptitle(f'sterinung ist {popt[0]}')
plt.show()

'''
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
'''
#plot rewsiduums
plt.loglog(dt,res)
plt.xlabel('dt')
plt.ylabel('x_num - x_exact')
plt.show()
