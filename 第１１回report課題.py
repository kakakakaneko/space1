import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from astropy import constants as const
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit
from scipy import optimize as opt

#####################-------constants-----------###################

pi = np.pi
G = const.G.value
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
M_jup = const.M_jup.value

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
plt.xlabel("Frequency(1/day)")
plt.ylabel("Power")
plt.show()
estimated_timescale = 1/frequency[np.argmax(velocity)]
print("Time_scale = ", estimated_timescale)



phase = (t/estimated_timescale)%1
plt.errorbar(phase,v,yerr=error,marker="o",capthick=1, capsize=10,lw = 1, linestyle="None", label="data",zorder=1)
plt.xlabel("Phase")
plt.ylabel("Velocity(km/s)")
#plt.show()

def func(t,a,b):
    return a*np.sin(2*pi*t+b)

popt, pcov = curve_fit(func,phase,v,p0=[0.70,0])
x_fit = np.linspace(0,1,100)
y_fit = func(x_fit,*popt)
vr = popt[0]*1000
plt.plot(x_fit,y_fit,label="fit",zorder=2,linewidth=3)
plt.legend()
plt.show()
print("A = ",popt[0]*1000,"m/s")

def function(e):
    return (1+e)**2-2*pi*M_sun*G*(e**3)/(estimated_timescale*24*60*60*(vr**3))

eps = opt.newton(function,0.6)
print("eps = ",eps)
print("D = ",M_sun*G*(eps**2)/((vr**2)*(1+eps))/au,"au")
print("M_acc/M_jup = ",eps*M_sun/M_jup)
print("M_jup/M_sun = ",M_jup/M_sun)

###--------------------------------------------------

# task4: max(min) is the A(振幅)
# task5:solve by using kepler's law
# 656.1122764193188 nm
# 486.00909364393993 nm
# 433.93669075351784 nm
# 410.0701727620744 nm

###--------------------------------------------------