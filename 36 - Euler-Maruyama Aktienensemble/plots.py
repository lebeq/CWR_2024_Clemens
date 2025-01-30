import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('stock_data.csv', float, delimiter=',',unpack=False)
time = data[:,0]
#die einzelnen stocks auslesen und plotten
stonky = []
for i in range(1,201):
    #stonky.append(data[:,i])
    #print(i,np.shape(stonky[i-1]))
    plt.plot(time,data[:,i], alpha = 0.5)
plt.xlabel('Zeit [Monate]')
plt.ylabel('Stonks')
plt.show()