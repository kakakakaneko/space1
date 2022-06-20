import numpy as np
from matplotlib import pyplot as plt
from astropy import constants as const

#####################-------constants-----------###################

pi = np.pi
r_sun = const.R_sun.value
r_earth = const.R_earth.value
au = const.au.value
sigma = const.sigma_sb.value
L_sun = const.L_sun.value
h = const.h.value
k = const.k_B.value
c = const.c.value
L_tra = 5.24e-4*L_sun

####################----------------------------###################

T_eq = []
dist = [0.01111,0.01521,0.02144,0.02817,0.0371,0.0451,0.063]
pla = ["b","c","d","e","f","g","h"]

for i in range(len(dist)):
    E_in = (L_tra/(4*pi*(dist[i]*au)**2))*pi
    T = (E_in/(sigma*4*pi))**(1/4)
    T_eq.append(T)
    
plt.title("TRAPPIST-1")
plt.xlabel("distance (au)")
plt.ylabel("T_eq(K)")
plt.scatter(dist, T_eq)
for i in range(len(dist)):
    plt.scatter([dist[i]],[T_eq[i]],label=pla[i])
plt.axhline(y=273.15, linestyle='dashed', color='Red')
plt.axhline(y=373.15, linestyle='dashed', color='Red')
plt.legend()
plt.show()