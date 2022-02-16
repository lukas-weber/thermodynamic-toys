import numpy as np
import matplotlib.pyplot as plt
import glob
import re

# mag^2
files = glob.glob("data*.csv")
for i, f in enumerate(files):
	L = int(f[7:-4])
	data = np.genfromtxt(f,skip_header=1)

	plt.errorbar(data[:,0],data[:,1],yerr=data[:,2],label="L=%g"%L)

plt.xlabel('$T/J$')
plt.ylabel('$M^2$')
plt.legend()
plt.show()


# susceptibility*L
files = glob.glob("data*.csv")
for i, f in enumerate(files):
	L = int(f[7:-4])
	data = np.genfromtxt(f,skip_header=1)

	chiL = data[:,1:]/data[:,0,np.newaxis]*L
	plt.errorbar(data[:,0],chiL[:,0],yerr=chiL[:,1],label="L=%g"%L)

plt.xlabel('$T/J$')
plt.ylabel('$\chi L$')
plt.legend()
plt.show()

