import piston
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo

steps = 1000
Npart = 2000

T = 200 
gravityg = 1
w = np.sqrt(3*gravityg**2/(Npart*T))

deltat = 5*2*np.pi/w/steps
X0 = 1.02

mu = 0.1 # this parameter is related to the ratio of period time to equilibration time. high values drive the system out of equilibrium. It is also the total mass of the gas. 

mpart = mu/Npart

t,x = piston.simulate(steps, deltat, Npart, mpart, T, X0, gravityg)
plt.plot(t*w/(2*np.pi),x)

def f(t, c, a, phi):
    return c + a*np.cos(w*t+phi)

popt, _ = spo.curve_fit(f, t, x,p0=[x.mean(),X0-x.mean(), 0])

plt.plot(t*w/(2*np.pi),f(t,*popt))
plt.xlabel(r'$t \omega_s/(2\pi)$')
plt.ylabel(r'$L/L_0$')
plt.show()

X0 = 1.5
T = 100
deltat*=4
for mu in [0.1, 1, 10]:
    t,x = piston.simulate(steps, deltat, Npart, mu/Npart, T, X0, gravityg)
    plt.plot(t*w/(2*np.pi),x, label='$\\mu={}$'.format(mu))

plt.xlabel(r'$t \omega_s/(2\pi)$')
plt.ylabel(r'$L/L_0$')
plt.legend()
plt.show()
