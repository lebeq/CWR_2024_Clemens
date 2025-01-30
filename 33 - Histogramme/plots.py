import numpy as np
from matplotlib import pyplot as plt

x,D = np.loadtxt("hist_data.csv", delimiter=',', unpack=True)

'''#Intervallbreiten
dx = y-x
#Gesamtzahl von Samples
N = np.sum(H)
#Relative Dichte
D = H/dx/N'''

#plotten
plt.bar(x,D, align='edge', width=(x[1]-x[0]), edgecolor='k')
plt.xlabel('Zufallszahl erzeugt mit Polarmethode')
plt.ylabel('Frequenz')
plt.savefig('.\Gaussian.jpg')
plt.show()
