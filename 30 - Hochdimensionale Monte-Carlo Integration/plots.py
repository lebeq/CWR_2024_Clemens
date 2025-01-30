import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import *

def linfunc(x,m,b):
    return m*x+b

'''
dim, data, var = np.loadtxt("mc_data.csv", delimiter=',', unpack=True)
ana_dim, ana_data = np.loadtxt("ana_data.csv",delimiter=',',unpack=True)
plt.errorbar(dim,data, yerr=var, alpha=0.5, fmt='o')
plt.plot(ana_dim,ana_data, 'r')
plt.show()
'''
N, data2, var2, abs_diff = np.loadtxt("hypersphere_data.csv", delimiter = ',', unpack=True)
ana_N, ana_dat = np.loadtxt("analytical_hypersphere.csv",delimiter=',', unpack=True)
plt.errorbar(N,data2,yerr=var2,alpha=0.5,fmt='o')
plt.plot(ana_N,ana_dat,'r')
plt.xlabel('N')
plt.ylabel('Volume')
plt.show()

#line fit zu abweichung
popt, pcov = curve_fit(linfunc,np.log(N), np.log(abs_diff))

plt.scatter(N,abs_diff)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('log(N)')
plt.ylabel('log(sqrt(var))')
plt.show()

plt.plot(np.log(N), linfunc(np.log(N),*popt))
plt.xlabel('log(N)')
plt.ylabel('linfunc')
plt.suptitle(f'steigung ist {popt[0]}')
plt.grid()
plt.show()
