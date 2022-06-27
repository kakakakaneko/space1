import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from astropy import constants as const
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit

#####################-------constants-----------###################

pi = np.pi
r_sun = const.R_sun.value
r_earth = const.R_earth.value
au = const.au.value
sigma = const.sigma_sb.value
L_sun = const.L_sun.value
M_sun = const.M_sun.value
h = const.h.value
k = const.k_B.value
c = const.c.value
ryd = const.Ryd.value

####################----------------------------###################

data = pd.read_csv('51peg.dat', header=None,delim_whitespace=True)
t = list(data[0])
v = list(data[1])
error = list(data[2])
vg = 0

for vs in v:
    vg += vs
    
vg = vg/(len(v))

print("Vg = ", vg)

for i in range(len(v)):
    v[i] = v[i]-vg

frequency, velocity = LombScargle(t, v).autopower(nyquist_factor=20)
plt.plot(frequency,velocity)
plt.xlabel("Frequency")
plt.ylabel("Velocity")
plt.show()
estimated_timescale = 1/frequency[np.argmax(velocity)]
print("Time_scale = ", estimated_timescale)



phase = (t/estimated_timescale)%1
plt.errorbar(phase,v,yerr=error,marker="o",capthick=1, capsize=10,lw = 1, linestyle="None", label="data",zorder=1)
plt.xlabel("phase")
plt.ylabel("velocity")
#plt.show()

def func(t,a,b):
    return a*np.sin(2*pi*t+b)

popt, pcov = curve_fit(func,phase,v,p0=[0.70,0])
x_fit = np.linspace(0,1,100)
y_fit = func(x_fit,*popt)
plt.plot(x_fit,y_fit,label="fit",zorder=2,linewidth=3)
plt.legend()
plt.show()
print("A = ",popt[0])

###--------------------------------------------------

# task4: max(min) is the A(振幅)
# task5:solve by using kepler's law
# 656.1122764193188 nm
# 486.00909364393993 nm
# 433.93669075351784 nm
# 410.0701727620744 nm

###--------------------------------------------------