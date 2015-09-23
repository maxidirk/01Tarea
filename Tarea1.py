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
'''
plt.xscale('log')
plt.yscale('log')
#plt.plot(np.log(w1),np.log(e1))  #Ln-Ln
#plt.plot(np.log10(w1),np.log10(e1))  #Log-Log
plt.plot(w1,e1)
plt.title('Flujo por $\lambda$ vs $\lambda$ en escala Log-Log')
plt.xlabel('Longitud de Onda ($\lambda$) $[\AA]$')
plt.ylabel('Flujo por $\lambda$ $[erg \ s^{-1} cm^{-2} \AA^{-1}]$')
plt.show()
'''

#Nota: El flujo es la integral de la Potencia en todo lambda (?)
#######################################################
#Parte de calcular la luminosidad
'''
Integral usando metodo del trapecio
int f(x) dx = dx/2 * (f(xo) + f(xo + dx))
en lo que sigue w1 es x y e1 es f(x)
'''
Ft = 0   #Flujo total a obtener

for i in range(n-1):
    Ft += ( (w1[i+1]-w1[i]) / 2.) * (e1[i] + e1[i+1])

print Ft #Flujo total en cgs

R = 1.5 * (10**13) # Distancia aproximada al Sol en cm
L = 4. * np.pi * (R**2) * Ft  # LUminosidad

print L







#######################################################
#Parte funcion de Plank





#######################################################
