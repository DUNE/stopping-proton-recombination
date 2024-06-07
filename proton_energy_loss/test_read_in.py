import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

data = np.load("energy_loss_vs_rr.npz")
print(data['rr'].shape, data['dEdx'].shape)

rrs = np.arange(data['rr'][0],data['rr'][-1],0.1)
rrs_upper = np.arange(data['rr_upper'][0],data['rr_upper'][-1],0.1)
rrs_lower = np.arange(data['rr_lower'][0],data['rr_lower'][-1],0.1)
fn = interpolate.interp1d(data['rr'],data['dEdx'],kind='cubic')
fn_lower = interpolate.interp1d(data['rr_lower'],data['dEdx'],kind='cubic')
fn_upper = interpolate.interp1d(data['rr_upper'],data['dEdx'],kind='cubic')

argo = lambda x: 17.0*(x**(-0.42))

cols = ['orchid','lightseagreen','mediumpurple','red']


plt.scatter(data['rr'],data['dEdx'],color=cols[0],label="histogram values")
plt.plot(rrs,fn(rrs),color=cols[0],label="interpolated")
plt.xlabel("Residual Range [cm]")
plt.ylabel("dE/dx [MeV/cm]")
plt.legend()

plt.figure()
plt.scatter(data['rr'],data['dEdx'],color=cols[0],label="histogram values")
plt.plot(rrs,fn(rrs),color=cols[0], label="interpolated")
plt.plot(rrs,argo(rrs), color=cols[3], label="ArgoNeuT Approx")
plt.legend()
plt.xlabel("Residual Range [cm]")
plt.ylabel("dE/dx [MeV/cm]")
plt.xlim((0.0,10.0))
plt.ylim((0.0,45.0))


plt.figure()
plt.scatter(data['rr'],data['dEdx'],color=cols[0],label="Histogram Values (bin center)")
plt.plot(rrs,fn(rrs),color=cols[0],label="Interpolated (bin center)")
plt.scatter(data['rr_lower'],data['dEdx'],color=cols[1],label="Histogram Values (bin lower edge)")
plt.plot(rrs_lower,fn_lower(rrs_lower),color=cols[1],label="Interpolated (bin lower edge)")
plt.scatter(data['rr_upper'],data['dEdx'],color=cols[2],label="Histogram Values (bin upper edge)")
plt.plot(rrs_upper,fn_upper(rrs_upper),color=cols[2],label="Interpolated (bin upper edge)")
plt.plot(rrs_lower,argo(rrs_lower),color=cols[3],label="ArgoNeuT Approximation")
plt.legend()
plt.xlim((0.0,10.0))
plt.ylim((0.0,45.0))
plt.xlabel("Residual Range [cm]")
plt.ylabel("dE/dx [MeV/cm]")
plt.title("Comparing Numerical to ArgoNeuT Approximation")

plt.figure()
plt.scatter(data['rr'],data['dEdx'],color=cols[0],label="Histogram Values (bin center)")
plt.plot(rrs,fn(rrs),color=cols[0],label="Interpolated (bin center)")
plt.scatter(data['rr_lower'],data['dEdx'],color=cols[1],label="Histogram Values (bin lower edge)")
plt.plot(rrs_lower,fn_lower(rrs_lower),color=cols[1],label="Interpolated (bin lower edge)")
plt.scatter(data['rr_upper'],data['dEdx'],color=cols[2],label="Histogram Values (bin upper edge)")
plt.plot(rrs_upper,fn_upper(rrs_upper),color=cols[2],label="Interpolated (bin upper edge)")
plt.plot(rrs_lower,argo(rrs_lower),color=cols[3],label="ArgoNeuT Approximation")
plt.legend()
plt.xlim((100.0,120.0))
plt.ylim((2.0,3.0))
plt.xlabel("Residual Range [cm]")
plt.ylabel("dE/dx [MeV/cm]")
plt.title("Comparing Numerical to ArgoNeuT Approximation")

plt.figure()
plt.scatter(data['rr'],data['dEdx'],color=cols[0],label="Histogram Values (bin center)")
plt.plot(rrs,fn(rrs),color=cols[0],label="Interpolated (bin center)")
plt.scatter(data['rr_lower'],data['dEdx'],color=cols[1],label="Histogram Values (bin lower edge)")
plt.plot(rrs_lower,fn_lower(rrs_lower),color=cols[1],label="Interpolated (bin lower edge)")
plt.scatter(data['rr_upper'],data['dEdx'],color=cols[2], label="Histogram Values (bin upper edge)")
plt.plot(rrs_upper,fn_upper(rrs_upper),color=cols[2],label="Interpolated (bin upper edge)")
plt.plot(rrs_lower,argo(rrs_lower),color=cols[3],label="ArgoNeuT Approximation")
plt.legend()
plt.xlim((2.5,50.0))
plt.ylim((3.0,13.0))
plt.xlabel("Residual Range [cm]")
plt.ylabel("dE/dx [MeV/cm]")
plt.title("Comparing Numerical to ArgoNeuT Approximation")


plt.show()
