import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize as op

plt.rcParams['errorbar.capsize'] = 3

# data = np.loadtxt("data.txt")

# X,Y = data[:,0],np.log10(np.abs(data[:,1]))

# plt.scatter(X,Y)
# plt.xlabel('$\log_{10}N_{points}$')
# plt.ylabel('$\log_{10}$Err')
# plt.show()

data = np.loadtxt("time_data.txt")

X,Y,sY = data[:,0],data[:,1],data[:,2]
Y_OMP,sY_OMP = data[:,3],data[:,4]

plt.errorbar(X,Y,sY,fmt = 'bx--', label='hand-written')
# plt.errorbar(X,Y_OMP,sY_OMP,fmt = 'rx--', label='OMP function')

plt.axvline(4, label='Number of cores',color = 'k', lw=0.2)
'''
fitting = lambda X,A,T,offset,delta : A*np.sin(2*np.pi/T*X + delta) + offset
X,Y,sY = X[3:],Y[3:],sY[3:]
seed = [2e-3,4,0.025,-np.pi/2]
coef,cov = op.curve_fit(fitting,X,Y,sigma = sY,p0 =seed)
X = np.linspace(1,X[-1],num = 500)
# plt.plot(X,fitting(X,*seed),label='Seed',alpha = 0.1)
plt.plot(X,fitting(X,*coef),label=f'{coef[0]:.2e}$\sin(2\pi/${coef[1]:.2e} + {coef[3]:.2e}) + {coef[2]:.2e}')
'''
plt.xlabel('$N_{integrals}$')
plt.ylabel('$Time$')
plt.legend()
plt.show()