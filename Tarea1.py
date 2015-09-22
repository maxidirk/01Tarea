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


print w
print e
