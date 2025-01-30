import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import *

#auslesen von Daten
data_const = np.loadtxt("p5_abs_err_const.csv", float, delimiter = ',')
data_var = np.loadtxt("p5_abs_err_var.csv", float, delimiter=',')
dt_const = data_const[:,0]
abs_err_const = data_const[:,1]
dt_var = data_var[:,0]
abs_err_var = data_var[:,1]

def linfunc(x,m,b):
    '''
    Für curve_fit
    f(x)=m*x+b
    m ist die Steigung aka das Konvergenzverhalten
    '''
    return m*x+b

#Plots von a=var und a=const
figure, axis = plt.subplots(1,2, figsize=(9,7))
axis[0].loglog(dt_var,abs_err_var, label = 'a = var')
axis[0].set_title('Höhenabhängige Beschleunigung')
axis[1].loglog(dt_const,abs_err_const, color = 'r', label = 'a = const')
axis[1].set_title('Konstante Beschleunigung')
figure.supxlabel('log(dt)')
figure.supylabel('log(Betrag von Fehler)')
figure.legend(loc=1)

#plots speichern und anzeigen
plt.savefig('./Abs_Error.png')
plt.show()

#fitting von linearer Funktion zum Data
popt_const, pcov_const = curve_fit(linfunc, np.log(dt_const), np.log(abs_err_const))
popt_var, pcov_var = curve_fit(linfunc, np.log(dt_const), np.log(abs_err_const))

#plotten von gefitteter Geraden
figure2, axis2 = plt.subplots(1,2, figsize=(9,7))
axis2[0].plot(np.log(dt_var), linfunc(np.log(dt_var), *popt_const), '--g')
axis2[0].set_title(f'Steigung für g const ist {popt_const[0]}')

axis2[1].plot(np.log(dt_var), linfunc(np.log(dt_var), *popt_var), '--g')
axis2[1].set_title(f'Steigung g var ist {popt_var[0]}')
figure2.supxlabel('log(x)')
figure2.supylabel('m*log(x)+b')

#plots speichern und anzeigen
plt.savefig('./fitted_data.png')
plt.show()

#Ausgabe von Steigung im Terminal 
print(popt_const[0], popt_var[0])
