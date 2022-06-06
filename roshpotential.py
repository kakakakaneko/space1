from scipy import optimize as opt
import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

m1 = 1
m2 = 0.001
mu1 = m1/(m1+m2)
mu2 = m2/(m1+m2)
e = 1e-3

def rosh(X,Y):
    return -mu1/np.sqrt((X+mu2)**2+Y**2)-mu2/np.sqrt((X-mu1)**2+Y**2)-(X**2+Y**2)/2

x = np.linspace(-1.25,1.25,250)
y = np.linspace(-1.25,1.25,250)
X, Y = np.meshgrid(x,y)
Z = rosh(X,Y)


fig, ax = plt.subplots()
cont_ls = ax.contour(X, Y, Z, np.linspace(-2,-1,100), cmap='PuOr')
ax.scatter([0.9313099885409696],[0])
ax.scatter([1.069892950509403],[0])
ax.scatter([-1.0004162503620317],[0])
ax.scatter([1/2-mu2],[np.sqrt(3)/2])
ax.scatter([1/2-mu2],[-np.sqrt(3)/2])

ax.set_aspect('equal','box')
plt.colorbar(cont_ls)
plt.savefig("contour_PuOr_linspace.png", bbox_inches = 'tight', pad_inches = 0.1)