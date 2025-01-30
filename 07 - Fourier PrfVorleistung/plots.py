import numpy as np
from matplotlib import pyplot as plt

'''
Hier ist der code der die Daten aus Aufgabe 7 plottet.
'''

############ Aufgabenteil 4 ################

figure, axis = plt.subplots(2,2)

#Daten auslesen
data10 = np.loadtxt("fourier_int_10.csv", delimiter=",") #M=10
data30 = np.loadtxt("fourier_int_30.csv", delimiter=",") #M=30
data100 = np.loadtxt("fourier_int_100.csv", delimiter=",") #M=100
data500 = np.loadtxt("fourier_int_500.csv", delimiter=",") #M=500
anadata = np.loadtxt("ana_sol.csv", delimiter = ",") #analytische LSG
imagdata = np.loadtxt('imag_fourier_int_100.csv', delimiter=',') #imaginärteil, M=100

# x, y, z = np.loadtxt("daten.csv", float, delimiter=",", unpack=True)

#analytische LSG x,y data
xdata_ana = anadata[:,0]
ydata_ana = anadata[:,1]

#Plot von M = 10
xdata10 = data10[:,0]
ydata10 = data10[:,1]
axis[0,0].plot(xdata10, ydata10, 'b-', label = 'numerische Lsg')
axis[0,0].plot(xdata_ana,ydata_ana, 'r-', label = 'analytische Lsg')
axis[0,0].set_title('M = 10')
#Plot von M = 30
xdata30 = data30[:,0]
ydata30 = data30[:,1]
axis[0,1].plot(xdata30, ydata30, 'b-', label = 'numerische Lsg')
axis[0,1].plot(xdata_ana,ydata_ana, 'r-', label = 'analytische Lsg')
axis[0,1].set_title('M = 30')
#Plot von M = 100
xdata100 = data100[:,0]
ydata100 = data100[:,1]
axis[1,0].plot(xdata100, ydata100, 'b-', label = 'numerische Lsg')
axis[1,0].plot(xdata_ana,ydata_ana, 'r-', label = 'analytische Lsg')
axis[1,0].set_title('M = 100')
#Plot von M = 500
xdata500 = data500[:,0]
ydata500 = data500[:,1]
axis[1,1].plot(xdata500, ydata500, 'b-', label = 'numerische Lsg')
axis[1,1].plot(xdata_ana,ydata_ana, 'r-', label = 'analytische Lsg')
axis[1,1].set_title('M = 500')

#Axenbeschriftung global
figure.supxlabel('k')
figure.supylabel('Realteil Integral')

#Layout finetuning
plt.subplots_adjust(left = 0.13, bottom = 0.1, hspace = 0.4)
handles, labels = plt.gca().get_legend_handles_labels() #wahrscheinlich nicht nötig, weil alle gleich, holt fir labels und farben von den subplots
figure.legend(handles, labels, loc = 1, fontsize = 'small') #fügt legende hinzu, loc = 1 steht für oben rechts

#Plots speichern und anzeigen
plt.savefig('./Realteil_vergleiche.png')
plt.show()


################ Aufgabenteil 5 ##############

imdat = np.loadtxt("imag_fourier_int_100.csv", delimiter=',')
xdata = imdat[:,0]
ydata = imdat[:,1]
plt.plot(xdata,ydata, label = 'Imaginärteil')
plt.xlabel('k')
plt.ylabel('Imaginärteil Integral')
plt.title('M=100')
plt.legend(loc = 1, fontsize = 'small')


plt.savefig('./Imaginärteil_plot.png')
plt.show()