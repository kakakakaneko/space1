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
e = -1e-5

def L1(x):
    return mu1*((x-mu1)**2)-mu2*((x+mu2)**2)-x*((x-mu1)**2)*((x+mu2)**2)
    
def L2(x):
    return mu1*((x-mu1)**2)+mu2*((x+mu2)**2)-x*((x-mu1)**2)*((x+mu2)**2)

def L3(x):
    return mu1*((x-mu1)**2)+mu2*((x+mu2)**2)+x*((x-mu1)**2)*((x+mu2)**2)

xinit = -5
L1solution = opt.newton(L1,xinit)
L2solution = opt.newton(L2,xinit)
L3solution = opt.newton(L3,xinit)
L4_x = 1/2 - mu2
L4_y = np.sqrt(3)/2
L5_x = 1/2 - mu2
L5_y = - np.sqrt(3)/2

print("x of L1 = ", L1solution)
print("x of L2 = ", L2solution)
print("x of L3 = ", L3solution)

def function(t,u):
    x, y, vx, vy = u
    r1 = ((x+mu2)**2 + y**2)**(1/2)
    r2 = ((x-mu1)**2 + y**2)**(1/2)
    d2xdt2 = 2*vy + x - (x+mu2)/(r1**3)*mu1-(x-mu1)/(r2**3)*mu2
    d2ydt2 = -2*vx + y*(1-mu1/(r1**3)-mu2/(r2**3))
    return [vx, vy, d2xdt2,d2ydt2]

tini = 0
tfin = 1000
domain = (tini, tfin)
t = np.linspace(tini, tfin, num=1000)
uini = [L3solution, 0, e, 0]
usol = solve_ivp(fun=function, t_span=domain,\
                 y0=uini, t_eval=t, rtol=1e-4, atol=1e-6)


x = usol.y[0]
y = usol.y[1]

  
plt.plot(x,y,c="blue")
plt.plot(mu1,0,marker = "o",color="green",label="jupiter")
plt.plot(-mu2,0,marker ="o",color="orange",label="sun")
plt.plot(L1solution,0,marker ="o",label="L1")
plt.plot(L2solution,0,marker ="o",label="L2")
plt.plot(L3solution,0,marker ="o",label="L3")
plt.plot(L4_x,L4_y,marker ="o",label="L4")
plt.plot(L5_x,L5_y,marker ="o",label="L5")
plt.legend()

######
# x of L1 =  0.9313099885409696
# x of L2 =  1.069892950509403
# x of L3 =  -1.0004162503620317
######