import numpy as np
import matplotlib.pyplot as plt


# data systematics (the new ones only apply to data)
alpha = 0.920
alpha_sigma = 0.009
beta = 0.195
beta_sigma = 0.004
Ab = 0.83
Ab_sigma = 0.013
kb = 0.053
kb_sigma = 0.003


# non-recombination calibration systematic
alpha_norec = 0.017
beta_norec = 0.026
Ab_norec = 0.010
kb_norec = 0.006


# SCE wrong correction map systematic
alpha_SCE = 0.004
beta_SCE = 0.002
Ab_SCE = 0.007
kb_SCE = 0.002

# E field systematc
alpha_E = 0.0
beta_E = 0.002
Ab_E = 0.0
kb_E = 0.001


# add in quadrature
alpha_error = np.sqrt(alpha_sigma**2 + alpha_norec**2 + alpha_SCE**2 + alpha_E**2)
beta_error = np.sqrt(beta_sigma**2 + beta_norec**2 + beta_SCE**2 + beta_E**2)
Ab_error = np.sqrt(Ab_sigma**2 + Ab_norec**2 + Ab_SCE**2 + Ab_E**2)
kb_error = np.sqrt(kb_sigma**2 + kb_norec**2 + kb_SCE**2 + kb_E**2)


print(f"  best fit parameters alpha: {alpha:.3f} +/-{alpha_error:.3f}, beta: {beta:.3f} +/-{beta_error:.3f}")
print(f"  best fit parameters A Birks: {Ab:.3f} +/-{Ab_error:.3f}, k Birks: {kb:.3f} +/-{kb_error:.3f}")



# figure out the % uncertainty in the energy resolution
W = 23.6*10**(-6)
rho = 1.383
# some detector specific parameters                                                      
eps = 0.553 


def dQdx(dEdx):
    return (1.0/(beta*W))*np.log((beta*dEdx) + alpha)

def dEdx(dEdx_prime,a,b):
    return (1.0/b)*(np.exp(b*W*dQdx(dEdx_prime))-a)


E_bins = np.arange(1.0,20.0,0.1) # MeV/cm


dEdx_upper = dEdx(E_bins,alpha+alpha_error,beta+beta_error)
dEdx_alpha = dEdx(E_bins,alpha+alpha_error,beta)
dEdx_beta = dEdx(E_bins,alpha,beta+beta_error)
dEdx_upa_downb = dEdx(E_bins,alpha+alpha_error,beta-beta_error)
dEdx_downa_upb = dEdx(E_bins,alpha-alpha_error,beta+beta_error)
dEdx_lower = dEdx(E_bins,alpha-alpha_error,beta-beta_error)

plt.figure()
plt.plot(E_bins,dQdx(E_bins))
plt.xlabel("dE/dx")
plt.ylabel("dQ/dx")


#plt.rc('text', usetex=True)
#plt.figure()
#plt.plot(E_bins, (dEdx_upper-E_bins)/E_bins, label=r"$\alpha+\sigma_{\alpha}, \beta+\sigma_{\beta}$")
#plt.plot(E_bins, np.zeros_like(E_bins), label="baseline")
#plt.plot(E_bins, (dEdx_lower-E_bins)/E_bins, label=r"$\alpha-\sigma_{\alpha}, \beta-\sigma_{\beta}$")
#plt.plot(E_bins, (dEdx_upa_downb-E_bins)/E_bins, label="up a, down b")
#plt.plot(E_bins, (dEdx_downa_upb-E_bins)/E_bins, label="down a, up b")
#plt.plot(E_bins, (dEdx_alpha-E_bins)/E_bins, label="up a")
#plt.plot(E_bins, (dEdx_beta-E_bins)/E_bins, label="up b")
#plt.title("dE/dx Resolution from Recombination")
#plt.xlabel("dE/dx [MeV/cm]")
#plt.ylabel("Fractional Uncertainty")
#plt.legend()


# do something better, throw n_expts toy experiments for the different dE/dx values
n_expts = 10000
alphas = np.random.normal(alpha,alpha_error,n_expts)
betas = np.random.normal(beta,beta_error,n_expts)

experiments = np.empty((n_expts,E_bins.shape[0]))


plt.figure()
for i in range(n_expts):
    dEdx_loop = dEdx(E_bins,alphas[i],betas[i])
    thing = (dEdx_loop-E_bins)/E_bins
    experiments[i] = thing
    #plt.plot(E_bins, thing, 'o', markersize=1,color="mediumseagreen",alpha=0.1)


lower_bounds = np.zeros_like(E_bins)
upper_bounds = np.zeros_like(E_bins)
for j in range(E_bins.shape[0]):
    the_slice = experiments[:,j]
    slice_mu = np.mean(the_slice)
    slice_sigma = np.sqrt(np.mean((the_slice-slice_mu)**2))
    lower_bounds[j] = slice_mu - slice_sigma
    upper_bounds[j] = slice_mu + slice_sigma

plt.plot(E_bins,lower_bounds,linewidth=2,color="red",label="lower bound one sigma")
plt.plot(E_bins,upper_bounds,linewidth=2,color="red",label="upper bound one sigma")
plt.ylim((-0.1,0.1))
plt.xlim((0.0,20.0))
plt.title("Recombination dE/dx Uncertainty")
plt.xlabel("dE/dx [MeV/cm]")
plt.ylabel("Fractional Uncertainty")
#plt.legend()
plt.savefig("dEdx_uncertainty.pdf")

plt.show()
