import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import scipy.interpolate as spinter
import scipy.integrate as spi

# 1D house with a heating wall at x=0 and a leaky wall at x=L
# the boltzmann transport equation for this scenario
# (1D collisionless gas) can be solved analytically up
# to an integral equation for the boundary condition.

L = 1
σ = 0.9

Tenv = 1
nenv = 1
m = 0.1



def fenv(v):
    return np.exp(-m*v**2/2/Tenv)/(2*np.pi*Tenv/m)**0.5



def Theat(t):
    rate = 0.1
    traise = 2
    
    return Tenv+rate*np.minimum((t-2),t*0+traise)*(t>2)
    
def integralabsvelheat(maxv,t=0):
    return 1/(2*np.pi*m/Theat(t))**0.5 * (1-np.exp(-m/Theat(t)/2*maxv**2))

meanabsvelenv = 2/(2*np.pi*m/Tenv)**0.5
meanabsvelheat = 2*integralabsvelheat(1e8)

ni0 = nenv*meanabsvelenv/meanabsvelheat

ts = np.linspace(1e-7,8,200)
nis = np.ones_like(ts)*ni0

def fheat(v,t):
    return np.exp(-m*v**2/2/Theat(t))/(2*np.pi*Theat(t)/m)**0.5

def residual(nis,ts,detailed=False):
    lhs = nis*integralabsvelheat(1e8,ts) - (1-σ)*nenv*meanabsvelenv/2 - σ*ni0*integralabsvelheat(2*L/ts)
    rhs = σ * (2*L)**2 * np.array([np.trapz(1/(t-ts[:it])**3*fheat(2*L/(t-ts[:it]),ts[:it])*ni0,ts[:it]) for it, t in enumerate(ts)])

    residual = np.sum((rhs-lhs)**2)
    if not detailed:
        return residual
    else:
        return rhs-lhs

res = spo.minimize(residual,nis,method='CG',args=(ts,))

nis = res.x
print(res.message)
print(f'Residual: {residual(nis,ts)}')

plt.plot(ts,nis*integralabsvelheat(1e8,ts), label='$n_i$')
plt.plot(ts,Theat(ts), label='$T^\mathrm{heat}$')
#plt.plot(ts,ni0+0*nis)
#plt.plot(ts,residual(nis*0+ni0,ts,True))
plt.legend()
plt.show()

def integN(nis, ts):
    nit = spinter.interp1d(ts,nis,'quadratic',fill_value=ni0, bounds_error=False)
    def integrand(v,x,t):
        return nit(t-x/v)*fheat(v,t-x/v) + σ*nit(t-(2*L-x)/v)*fheat(v,t-(2*L-x)/v)

    N = 0*ts
    Nerr = 0*ts
    for i, t in enumerate(ts):
        N[i], Nerr[i] = spi.dblquad(integrand, 0, L, lambda x: 0, lambda x: 10*(Tenv/m), epsabs=1e-3, epsrel=1e-3, args=(t,)) 
    return N+(1-σ)*nenv/2, Nerr

N, Nerr = integN(nis, ts)

plt.errorbar(ts, N, Nerr)
plt.show()
