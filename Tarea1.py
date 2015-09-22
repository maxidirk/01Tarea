#Tarea01
#Maximiliano Dirk Vega Aguilera 18.451.231-9

#######################################################

import numpy as np
#import matplotlib.pyplot as plt
from pylab import *
import astropy as ast
import scipy as sci

#######################################################

a=np.loadtxt('sun_AM0.dat')



plot(log(a[0][:]),log(a[1][:]))
show()

	

