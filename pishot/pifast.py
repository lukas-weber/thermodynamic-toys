import numpy as np
import matplotlib.pyplot as plt

def measure_pi(Nsamples,Ntries=500):
    x = np.random.rand(Ntries,Nsamples,2)
    pi = np.mean(np.sum(x**2, axis=2) < 1, axis=1)*4
    sigmaPi = np.std(pi, ddof=1)
    return sigmaPi

Ns = np.logspace(1, 4, 100, dtype=int)
sigmaPis = np.zeros(len(Ns))

for i, N in enumerate(Ns):
    sigmaPis[i] = measure_pi(N)

plt.plot(Ns, sigmaPis, '.', label='Error')
plt.plot(Ns, 1/Ns**0.5, label='$1/\sqrt{N}$')

# p = np.pi/4
# plt.plot(Ns, 4./Ns * np.sqrt(Ns*p*(1-p)))

plt.xlabel('N')
plt.ylabel('$\sigma_\pi$')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
