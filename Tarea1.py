#######################################################
'''
Metodos Numericos para la Ciencia e Ingenieria FI3104
Tarea01
Maximiliano Dirk Vega Aguilera
18.451.231-9
'''
#######################################################

import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as ac
import scipy as sci

#######################################################

a=np.loadtxt('sun_AM0.dat')  # Se abre el archivo

n = len(a)
w = np.zeros(n)     # w = Wavelength (nm)
e = np.zeros(n)     # e = Power (W*m-2*nm-1)
w = a[:,0]          # obtiene w de a
e = a[:,1]          # obtiene e de a

w1 = 10. * w          # w1 = Wavelength (Angstrom)
e1 = e * ((10 ** (7)) * (10. ** (-4.)) * (10 ** (-1)))  # e1 = Power (erg*s-1*cm-2*A-1)

plt.plot(np.log(w1),np.log(e1))
plt.plot(np.log10(w1),np.log10(e1))
plt.title('$\lambda$ vs Potencia en escala Log-Log')
plt.xlabel('Longitud de Onda ($\lambda$) $[\AA]$')
plt.ylabel('Potencia  $[erg \ s^{-1} cm^{-2} \AA^{-1}]$')
plt.show()

#Nota: El flujo es la integral de la Potencia en todo lambda
#######################################################
#Parte de calcular la luminosidad




#######################################################
#Parte funcion de Plank





#######################################################
