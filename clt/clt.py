import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps

def check_histogram(name, rand, mean, std, bins):
    Nsamples, Nsum = rand.shape
    sums = np.sum(rand,axis=1)

    mean *= Nsum
    std *= Nsum**0.5

    x = mean+np.linspace(-5*std, 5*std,100)
    plt.title('{}: $N_{{samples}}={}$, $N_{{sum}} = {}$'.format(name,Nsamples,Nsum))
    plt.hist(sums, bins=20, density=True)
    plt.plot(x, sps.norm.pdf(x, loc=mean, scale=std), label='Gaussian')
    plt.xlabel('X')
    plt.ylabel('p(X)')
    plt.legend()
    plt.show()

Nsamples = 100000
Nsum = 300


xcoin = np.random.randint(0,2,(Nsamples,Nsum))*2-1
check_histogram('Coin', rand=xcoin, mean=0, std=1, bins=18)

xuniform = np.random.rand(Nsamples,Nsum)*2-1
check_histogram('Uniform', rand=xuniform, mean=0, std=2/12**0.5, bins=18)

λ = 2
xpoisson = np.random.poisson(λ, (Nsamples,Nsum))
check_histogram('Poisson',rand=xpoisson, mean=λ, std=λ**0.5, bins=18)
