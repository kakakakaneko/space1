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
pc = const.pc.value

####################----------------------------###################


data = pd.read_csv('SgrAstardata.dat', header=None)
t = list(data[0])
x = list(data[1])
x_error = list(data[2])
y = list(data[3])
y_error = list(data[4])



r = []
theta = []

for i in range(len(t)):
    r.append(np.sqrt(x[i]**2+y[i]**2))
    if (y[i]<=0)&(i<=15):
        theta.append(np.arccos(x[i]/r[i]))
    elif (y[i]>=0)&(i<=15):
        theta.append(2*pi-np.arccos(x[i]/r[i]))
    else:
        theta.append(2*pi+np.arccos(x[i]/r[i]))
    


plt.axes().set_aspect("equal")
plt.grid()
plt.errorbar(x,y,xerr=x_error,yerr=y_error,marker="o",capthick=1, capsize=10,lw = 1, label="Raw data",linestyle="None",zorder=1)
plt.xlabel("x(arcsec)")
plt.ylabel("y(arcsec)")


def func(t,a,b,eps):
    return a/(1+eps*np.cos(t+b))

popt, pcov = curve_fit(func,theta,r)
t_fit = np.linspace(0,3*pi,100)
r_fit = func(t_fit,*popt)
x_fit = []
y_fit = []
for i in range(len(t_fit)):
    x_fit.append(r_fit[i]*np.cos(-t_fit[i]))
    y_fit.append(r_fit[i]*np.sin(-t_fit[i]))

        
semimajor = (func(2*pi-popt[1],*popt)+func(pi-popt[1],*popt))/2
print("epsilon = ",popt[2])
print("semi-major axis(ac) = ",semimajor,"arcsec")
plt.plot(x_fit,y_fit,label="Curve fit",zorder=2,linewidth=3)
plt.legend()
plt.show()

a_real = 8.5e3*pc*np.tan(semimajor/3600/180*pi)/au
print("semi-major axis(au) = ",a_real,"au")

index  = [0,0]
index[0] = np.argmax(r)
index[1] = np.argmin(r)

plt.plot(t_fit,r_fit,label="Curve fit")
plt.scatter(theta,r,label="Raw data")
plt.xlabel("φ")
plt.ylabel("r(arcsec)")
plt.legend()
plt.show()

time_scale = abs(t[index[0]]-t[index[1]+1])*2
print("time_scale = ",time_scale,"year")

M_sgrA = ((2*pi)**2)*((a_real*au)**3)/G/(time_scale*365*24*3600)**2/M_sun
print("M_sgrA = ",M_sgrA,"M_sun")

R_g = 2*G*M_sgrA*M_sun/(c**2)
print("R_g = ",R_g/au,"au")
r_near = 8.5e3*pc*np.tan(func(2*pi-popt[1],*popt)/3600/180*pi)
print("r_near/R_g = ",r_near/R_g)
###--------------------------------------------------

# task4: max(min) is the A(振幅)
# task5:solve by using kepler's law
# 656.1122764193188 nm
# 486.00909364393993 nm
# 433.93669075351784 nm
# 410.0701727620744 nm

###--------------------------------------------------