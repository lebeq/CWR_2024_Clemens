import numpy as np
from matplotlib import pyplot as plt


data = np.loadtxt("heat_data.csv", delimiter=',')
plt.imshow(data, aspect = 'auto')
plt.show()

plt.plot(data[::10,:].T, alpha=0.5)
plt.show()