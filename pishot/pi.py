import numpy as np
import random
import matplotlib.pyplot as plt

def measure_pi(Nsamples, Ntries=500):
    pis = []
    for _ in range(Ntries):
        hits = 0
        for _ in range(Nsamples):
            x = random.random()
            y = random.random()

            if x**2 + y**2 < 1:
                hits += 1
        pis.append(4*hits/Nsamples)
    sigmaPi = np.std(pis, ddof=1)
    return sigmaPi

Ns = np.logspace(1, 4, 100, dtype=int)
sigmaPis = np.zeros(len(Ns))

for i, N in enumerate(Ns):
    sigmaPis[i] = measure_pi(N)

plt.plot(Ns, sigmaPis, '.', label='Error')
plt.plot(Ns, 1/Ns**0.5, label='$1/\sqrt{N}$')

plt.xlabel('N')
plt.ylabel('$\sigma_\pi$')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
