import numpy as np
from matplotlib import pyplot as plt

F = np.loadtxt('potential_data.csv', float, delimiter=',')
plt.matshow(F,10)
plt.contour(F, 14, colors='black')
plt.show()