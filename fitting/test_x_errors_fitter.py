import numpy as np
import matplotlib.pyplot as plt
import x_errors_fitter as xef


W = 23.6*10**(-6)
rho = 1.383
eps = 0.553 

def birks(p, x):
    A_b = p[0]
    k_b = p[1]
    return (A_b/W)*(x/(1+(k_b/(eps*rho))*x))

p0 = np.array([0.8,0.0486]) # starting params for fit
p_to_fit = np.array([0.84,0.05]) # fake data
x = np.arange(10,20,1.0)
y = birks(p_to_fit,x)
n = x.shape[0]


x_cov = np.identity(n)*x*0.01
y_cov = np.identity(n)*(y*0.02)**2


plt.figure()
plt.errorbar(x,y,yerr=np.sqrt(np.diag(y_cov)),xerr=np.sqrt(np.diag(x_cov)),linestyle='None',linewidth=2)

params, param_sigma, cov, chi2_ndof = xef.fit(birks, x, x_cov, y, y_cov, p0)
print(params, param_sigma, cov, chi2_ndof)
y_fit = birks(params,x)
plt.plot(x,y_fit)

plt.show()
