import numpy as np
from matplotlib import pyplot as plt

steppers = ['euler_data.csv', 'RK2_data.csv', 'RK4_data.csv', 'VV_data.csv']
energies = ['euler_energy.csv', 'RK2_energy.csv', 'RK4_energy.csv', 'VV_energy.csv']

for i in steppers:
    m1_x, m1_y, m2_x, m2_y, m3_x, m3_y = np.loadtxt(i, float, delimiter=',', unpack = True)
    plt.plot(m1_x, m1_y, label='Masse 1')
    plt.plot(m2_x,m2_y, 'r', label='Masse 2')
    plt.plot(m3_x,m3_y, 'g', label='Masse 3')
    plt.suptitle(i.replace('.csv', ''))
    filename = i.replace('data.csv', 'Trajektorie.png')
    plt.savefig(f'./{filename}')
    plt.legend(loc = 1, fontsize = 'small')
    plt.show()


fig, ax = plt.subplots(2,2, figsize=(10,8))
for i in range(4):
    t, energy = np.loadtxt(energies[i], float, delimiter=',', unpack = True)
    max_energy = np.amax(energy)
    min_energy = np.amin(energy)
    ax[i//2,i%2].plot(t, energy, label='E_total')
    ax[i//2,i%2].axhline(y=max_energy, color = 'g', linestyle= '--', label='maximum')
    ax[i//2,i%2].set_ylim(min_energy,max_energy)
    ax[i//2,i%2].set_title(energies[i].replace('.csv', ''))
fig.tight_layout()
handles, labels = plt.gca().get_legend_handles_labels()
fig.legend(handles, labels, loc=1, fontsize = 'small')
fig.supxlabel('t [s]')
fig.supylabel('E_tot [J]')
plt.autoscale(False)
plt.savefig('./Energies.png')
plt.show()

figure, axis = plt.subplots(2,2)
for i in range(4):
    m1_x, m1_y, m2_x, m2_y, m3_x, m3_y = np.loadtxt(steppers[i], float, delimiter=',', unpack = True)
    axis[i//2,i%2].plot(m1_x,m1_y, label='Masse 1')
    axis[i//2,i%2].plot(m2_x,m2_y, 'r', label='Masse 2')
    axis[i//2,i%2].plot(m3_x,m3_y, 'g', label='Masse 3')
    axis[i//2,i%2].set_xlabel('x(t) [m]')
    axis[i//2,i%2].set_ylabel('y(t) [m]')
    axis[i//2,i%2].set_title(steppers[i].replace('.csv', ''))
figure.tight_layout()
handles, labels = plt.gca().get_legend_handles_labels()
figure.legend(handles, labels, loc = 1, fontsize = 'small')
plt.savefig('./Zusammensetztung_Trajektortien.png')
plt.show()