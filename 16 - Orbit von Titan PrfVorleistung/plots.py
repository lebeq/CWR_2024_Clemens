import numpy as np
from matplotlib import pyplot as plt

dateien = ['titan_data.csv', 'exc_titan_data.csv', 'adapt_aphelion.csv']
plot_png = ['./Teil3.png', './Teil4.png', './Teil5.png']

for i in dateien:
    t, x, y, r, v = np.loadtxt(i, float, delimiter=',',unpack=True)
    
    figure, axis = plt.subplots(1,1, figsize=(6,6))
    #figure(figsize = (6,6))
    axis.plot(x,y, label='trajektorie')
    axis.set_xlabel('x(t)')
    axis.set_ylabel('y(t)')
    axis.set_title('Trajektorie')
    trstr = i.replace('.csv', '') + 'Trajektorie.png'
    figure.tight_layout()
    plt.savefig(trstr)
    plt.show()
    
    fig, ax = plt.subplots(1,2)

    ax[0].plot(t,r, label='Abs. Abstand')
    ax[0].set_xlabel('t')
    ax[0].set_ylabel('r(t)')
    ax[0].set_title('Abstand Titan - Saturn')

    ax[1].plot(t,v, label='Geschw. Titan')
    ax[1].set_xlabel('t')
    ax[1].set_ylabel('v(t)')
    ax[1].set_title('Geschwindigkeit von Titan')

    fig.tight_layout()
    str = i.replace('csv', 'png')
    plt.savefig(str)
    plt.show()