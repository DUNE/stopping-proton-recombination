import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import interpolate

file = uproot.open("/vols/dune/awaldron/protons/dEdx_from_Mike/DEDXvsRRvsPitch.root")

# to find the most probable value of dEdx
def get_MPV(profile,bins):
    bin_width = bins[1] - bins[0]
    bins = bins[:-1]
    fn = interpolate.interp1d(bins,profile,kind='cubic')
    m = np.argmax(profile)
    mpv = optimize.fmin(lambda x: -1.0*fn(x), 0.5*(bins[m]+bins[m+1])) 
    return mpv



protons = file['proton_DEDXvsRRvsPitch']
values = protons.values()
rr = protons.axis(0).edges()
rr_centers = rr + 0.5*(rr[1]-rr[0])
rr_negshift = rr - 0.5*(rr[1]-rr[0])
rr_negshift = rr_negshift[:-1]
rr_centers = rr_centers[:-1]
rr_lower = rr[:-1]
rr_upper = rr[1:]
print(rr_centers,rr)
dEdx = protons.axis(1).edges() # I think these are the right way around
dEdx_centers = dEdx + 0.5*(dEdx[1]-dEdx[0])
pitch = protons.axis(2).edges()

# find the appropriate pitch
target_pitch = 0.7 # cm, we have cut 0.5 - 0.7 cm.  Two relevant bins, try?
target_i = 0
for i in range(pitch.shape[0]):
    if pitch[i] < target_pitch:
        target_i = i
    else:
        break

dEdx_vs_rr = values[:,:,target_i]
dEdx_vs_rr = np.swapaxes(dEdx_vs_rr,0,1)
print(pitch[target_i], pitch[target_i+1])


arranged = np.arange(0.0,300.0,0.01)
plt.pcolormesh(rr_negshift,dEdx_centers,dEdx_vs_rr)
plt.plot(arranged,17.0*arranged**(-0.42),color="red",label="ArgoNeuT Approximation")
plt.legend()
plt.xlim(0.0,10.0)
plt.ylim(0.0,45.0)
plt.xlabel("Residual Range [cm]")
plt.ylabel("dE/dx [MeV/cm]")
plt.title("Numerical Energy Loss Compared to ArgoNeuT Approximation")
plt.savefig("num_vs_argo.png")

MPVs = np.zeros_like(rr)
means = np.zeros_like(rr)

# now look at the dEdx profiles
for j in range(rr.shape[0]-1):
    profile = dEdx_vs_rr[:,j]
    #plt.figure()
    #plt.plot(dEdx[1:],profile)
    MPV = get_MPV(profile,dEdx)
    MPVs[j] = MPV
    means[j] = np.average(dEdx[:-1],weights=profile)
    

plt.figure()
plt.plot(rr_centers,MPVs[:-1],label="Numerical (MPV)")
plt.plot(rr_centers,means[:-1],color="gold",linestyle='--',label="Numerical (mean)")
plt.plot(arranged,17.0*arranged**(-0.42),color="red",label="ArgoNeuT Approximation")
plt.xlim(0.0,10.0)
plt.ylim(0.0,45.0)
plt.xlabel("Residual Range [cm]")
plt.ylabel("dE/dx [MeV/cm]")
plt.legend()
plt.savefig("num_MPV_vs_Argo.pdf")

print(rr[30],MPVs[30]) # just to compare things.


# now look at the differences
plt.figure()
plt.plot(arranged,np.zeros_like(arranged),color="midnightblue",linestyle="--")
centers_argo = 17.0*rr_centers**(-0.42)
plt.plot(rr_centers,100.0*(MPVs[:-1]-centers_argo)/centers_argo,color="midnightblue")
plt.xlim(0.0,120.0)
plt.ylim(-6.0,6.0)
plt.xlabel("Residual Range [cm]")
plt.ylabel("(Numerical - ArgoNeuT)/ArgoNeuT [%]")
plt.savefig("residuals_full.pdf")

# now look at the differences
plt.figure()
plt.plot(arranged,np.zeros_like(arranged),color="midnightblue",linestyle="--")
centers_argo = 17.0*rr_centers**(-0.42)
plt.plot(rr_centers,100.0*(MPVs[:-1]-centers_argo)/centers_argo,color="midnightblue")
plt.xlim(2.5,50.0)
plt.ylim(-6.0,6.0)
plt.xlabel("Residual Range [cm]")
plt.ylabel("(Numerical - ArgoNeuT)/ArgoNeuT [%]")
plt.savefig("residuals.pdf")


# output the values
#rr_centers, MPVs
np.savez("energy_loss_vs_rr.npz", rr=rr_centers, rr_lower=rr_lower, rr_upper=rr_upper, dEdx=MPVs[:-1])

plt.show()


