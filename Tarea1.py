'''
Metodos Numericos para la Ciencia e Ingenieria FI3104
Tarea01
Maximiliano Dirk Vega Aguilera
18.451.231-9
'''
#######################################################

import numpy as np
import matplotlib.pyplot as plt
import astropy as ast
import scipy as sci

#######################################################

a=np.loadtxt('sun_AM0.dat')  # Se abre el archivo

n = len(a)
w = np.zeros(n)     # w = Wavelength (nm)
e = np.zeros(n)     # e = Power (W*m-2*nm-1)
w = a[:,0]
e = a[:,1]

w1 = 10. * w          # w1 = Wavelength (Angstrom)
e1 = e * ((10 ** (7)) * (10. ** (-4.)) * (10 ** (-1)))  # e1 = Power (erg*s-1*cm-2*A-1)


plt.plot(np.log(w1),np.log(e1))
plt.xlabel('Longitud de Onda $[\AA]$')
plt.ylabel('Potencia $[erg \ s^{-1} cm^{-2} \AA^{-1}]$')
plt.show()
