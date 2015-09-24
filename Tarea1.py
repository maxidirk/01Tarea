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
from astropy import units as au
from scipy import integrate as sciI

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
plt.plot(w1,e1)
plt.title('Flujo por $\lambda$ vs $\lambda$ en escala Log-Log')
plt.xlabel('Longitud de Onda ($\lambda$) $[\AA]$')
plt.ylabel('Flujo por $\lambda$ $[erg \ s^{-1} cm^{-2} \AA^{-1}]$')
plt.show()
'''

#######################################################

'''
Integral usando metodo del trapecio
int f(x) dx = dx/2 * (f(xo) + f(xo + dx))
en lo que sigue w1 es x y e1 es f(x)
'''

Ft = 0   #Flujo total a obtener

for i in range(n-1):  #Se calcula la integral
    Ft += ( (w1[i+1]-w1[i]) / 2.) * (e1[i] + e1[i+1])

Ft = Ft * au.erg / (au.s*(au.cm**2))   #Se agregan las unidades
print 'Flujo Total Recibido= ', Ft #Flujo recibido total en cgs

'''
Usando L=2*pi*d^2*Ft
'''
d = 1.5 * (10**13) * au.cm   # Distancia aproximada al Sol en cm
L = 4. * np.pi * (d**2) * Ft   # Luminosidad total en cgs

print 'Luminosidad Total = ', L


#######################################################

'''
Flujo emitido a partir de la ecuacion de Plank
P=(2pih/c^2)*(KBT/h)^4 * int x^3/(e^x -1)
'''

T = 5778.*au.K #Temperatura del Sol [K]
Cp = ((2.*np.pi*ac.h.cgs)/(ac.c.cgs**2)) * (((ac.k_B.cgs*T)/ac.h.cgs)**4) #Termino cte de P
Ip = 0 #Termino integral de P

'''
con el cambio de variable x=tan(x) => dx=sec(x)^2
usando el metodo del trapecio
se puede intergrar la funcion:
tan(x)^3/( exp( tan(x) ) -1 ) * 1/cos(x)^2
entre 0 y pi/2
'''

def f_Ip (x):  #Se define la funcion a integrar
    return (np.tan(x)**3/(np.exp(np.tan(x))-1)) *(1/(np.cos(x)**2))

tol = 0.001   #Cuanto quiero acercarme al limite
paso= 0.001   #Tamanho de los intervalos a sumar
x = np.arange(tol , np.pi/2. -tol  , paso) #arreglo con los argumentos a usar

for i in range(len(x)-1):   #Se calcula la integral
    Ip += ((x[i+1]-x[i]))/2.*(f_Ip(x[i+1])+f_Ip(x[i]))

P = Cp * Ip  #Flujo Total Emitido

print 'Flujo Total Emitido = ', P #Flujo total emitido

'''
Usando L= 4*pi*R^2*P
'''
R = np.sqrt(L/(4.*np.pi*P)) #Radio del Sol obtenido a partir de l y P calulados

print 'Radio del Sol = ', R #Radio del Sol

#######################################################

'''
A continuacion se calculan las integrales calculadas anteriormente usando
los metodos de integracion predefinidos en scipy
'''

f_Ipx2 = lambda x: (x**3/(np.exp(x)-1)) #no se usa la version x=tan(x) porque tira un error

Ft_sci_trapz = sciI.trapz(e1, x=w1) * au.erg / (au.s*(au.cm**2)) #Primera integral

fiP_sci_quad = Cp * sciI.quad(f_Ipx2, 0, 100)[0]  #Segunda integral * termino constante


print 'Flujo Recibido Calculado por scipy = ', Ft_sci_trapz

print 'Flujo Emitido Calculado por scipy = ', fiP_sci_quad

#Diferencia entre lo calculado por mi y el calculado por scipy
dFr = Ft - Ft_sci_trapz   #Diferencia de flujo recibido
dFe = P - fiP_sci_quad   #Diferencia de flujo emitido

print 'Diferencia de Flujo Recibido = ', dFr
print 'Diferencia de Flujo Emitido = ', dFe

#######################################################

'''
Velocidad de ejecucion obtenida con timeit:

Flujo Recibido Mio =

Flujo Recibido scipy = 1000 loops, best of 3: 324 micros per loop

Flujo Emitido Mio =

Flujo Emitido scipy = 100 loops, best of 3: 4.07 ms per loop

'''

#######################################################
