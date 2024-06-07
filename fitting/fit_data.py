import numpy as np
import uproot3 as uproot
import math
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import fit_helpers as fh
from matplotlib.colors import LogNorm
import scipy.odr as od
from scipy.ndimage import gaussian_filter1d
from rr_to_dEdx_converter import converter
import x_errors_fitter as xef

USE_MC = False
USE_BIRKS = True

# for saving plots
TAIL = "data"
if USE_MC:
    TAIL = "mc"

# some external parameters of nature        
W = 23.6*10**(-6)
rho = 1.383
# some detector specific parameters                 
eps = 0.553 # min 0.549, max 0.558, mean = 0.553
if USE_MC:
    eps = 0.4867


# whether to use fixed width bins or binned by end point resolution
USE_FIXED_BINS = False
NBINS = 30
#SIGMA_RR = 0.5 # uncertainty in residual range (wire spaceing for now
SIGMA_RR = 0.3 # uncertainty in residual range


# yz direction: 1%, 1.5% for data
# x direction: 0.3%
sigma2_dqdx_uncorr = 0.015**2
sigma2_dqdx_corr = 0.003**2

if USE_MC:
    sigma2_dqdx_uncorr = 0.01**2 
    sigma_dqdx_corr = 0.003**2

sigma2_dqdx_tot = sigma2_dqdx_uncorr + sigma2_dqdx_corr 
#SIGMA_DQDX = np.sqrt(0.003**2 + 0.015**2)

#MIN_BIN_WIDTH = 0.2 # (in dE/dx) only used in variable binning
MIN_BIN_WIDTH = 0.3
#MIN_BIN_WIDTH = 0.02 # (in dE/dx) only used in variable binning

PLOT_PROFILES = False
REMOVE_FAILED_FITS = True # false to keep the fall back values
USE_UNCALIBRATED_DQDX = False


def make_covariance(corr_sig2, sig2_diag, values):
    # the diagonal part is already correct
    n = values.shape[0]
    diag = np.identity(n)*sig2_diag
    # the off-diagonal part needs to be the values**2 times the 
    #  %**2 correlated uncertainty (just a number)
    off_diag = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            print(i,j)
            off_diag[i,j] = values[i]*values[j]*corr_sig2
    cov = diag + off_diag
    print(cov)
    return cov, np.sqrt(sig2_diag)

def make_covariance_dE(sigma):
    n = sigma.shape[0]
    cov = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            cov[i][j] = sigma[i]*sigma[j]
    return cov

def recalculate_dedx_hyp(df):
    df['rr_shift'] = df['rr_shift']*(1.0 - df['rr_fit_is_good'])
    good_mask = (df['rr_chi2'] < 1e11)
    df['rr_shift'] = df['rr_shift']*good_mask
    df['resrange_new'] = df['resrange'] - df['rr_shift']
    df['dedx_hyp'] = (df['resrange_new']**(-0.42))*17.0
    return df


def get_df_from_file(use_mc):
    #mc_file = "../output/mc_proton_recombination.root"
    #mc_file = "/vols/dune/awaldron/protons/root_files/mc_proton_recombination.root" # STANDARD - USE
   # mc_file = "/vols/dune/awaldron/protoDUNE/sungbin_endpoint/mc_proton_recombination_0_1mm_full_stats.root" # STANDARD - USE
    #mc_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/mc_proton_recombination_0_1mm_1cal_const_argoneut_box.root" # STANDARD - USE
    #mc_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/mc_proton_recombination_0_1mm_argoneut_box.root" # STANDARD - USE
    mc_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/mc_proton_recombination_my_cal.root" # STANDARD - USE
    #mc_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/mc_proton_recombination.root" # STANDARD - USE
    #mc_file = "/vols/dune/awaldron/protons/root_files/mc_proton_recombination_deltas.root" # DELTA STUDY
    #mc_file = "/vols/dune/awaldron/protons/root_files/mc_proton_recombination_wrongSCEcorrection.root" # SYSTEMATIC STUDY
    # sce off for checking
    #mc_file = "data/mc_proton_recombination_sce_off.root"
    #data_file = "../output/data_proton_recombination.root"
    #data_file = "/vols/dune/awaldron/protons/root_files/data_proton_recombination.root"
    #data_file = "/vols/dune/awaldron/protoDUNE/sungbin_endpoint/data_proton_recombination_2mm.root"
    #data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_0_1mm_argoneut_box.root"
    #data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_0_1mm_new_box.root"
    #data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_box_argoneut_2.root" # THIS NEEDS REMAKING
    #data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_0_1mm_1cal_const_iter_3.root" # THIS NEEDS REMAKING
    #data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_0_1mm_1cal_const_norecocal.root"
    data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_box_iter_3.root"
    #data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_endpoint_test_box_iter.root" # THIS NEEDS REMAKING
    #data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_0_1mm_1cal_const_argoneut_box.root"
    if USE_BIRKS:
        #data_file = "/vols/dune/awaldron/protons/root_files/data_proton_recombination_birks.root"
        #data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_0_1mm_birks.root"
        data_file = "/dune/app/users/waldron/dev/2023/April/stopping-protons/output/data_proton_recombination_birks_iter_2.root"
    filename = data_file
    if use_mc:
        filename = mc_file
    tree = uproot.open(filename)['output_tree']
    df = tree.pandas.df()
    # hack for sce off
    #df = df.rename(columns={'cali_dqdx':'cali_dqdx_sc_ON','cali_dqdx_sc_off':'cali_dqdx'})

    # use better rr to dEdx conversion
    conv = converter()
    # get the new residual range from Pip's fitting (needed to get resrange_new)
    # DOESN'T SEEM TO BE WORKING (LEAVE COMMENTED OUT)
    #df = recalculate_dedx_hyp(df)    
    #df['dedx_hyp'] = conv.get_dEdx(df['resrange_new'])


    # adjust the endpoint using Sungbin's shift
    df['resrange_new'] = df['resrange'] + df['rr_shift_sungbin']
    df['resrange'] = df['resrange_new']
    df['dedx_hyp'] = conv.get_dEdx(df['resrange'])
    return df


def apply_quality_cuts(df):
    # cut on end of tracks
    #df = df[df['resrange'] > 0.5]
    df = df[df['resrange'] > 1.0]

    # cut on angle - why does this work?
    df = df[df['pitch'] < 0.64]
  
    # cut whole event if suspect pion
    df = df[df['rr_shift_sungbin'] < 3.0] # in mc there is not much past 3cm of end

    # also to remove suspected pions
    df = df[df['chi2_sungbin'] < 10.0]

    # cut on dx
    df = df[df['dx'] > 0.5]
    df = df[df['dx'] < 0.7]

    # cut on number of hits (reco failure? hardly any here)
    df = df[df['nhits'] < 300]

    # try to deal with split tracks where the stopping part missing
    good_events = []
    for event in set(df['event']):
        hits = df[df['event'] == event]
        max_range = np.max(hits['resrange'])
        end_med = np.median(hits[hits['resrange'] < 5.0]['cali_dqdx'])
        start_med = np.median(hits[hits['resrange'] > max_range - 5.0]['cali_dqdx'])
        if end_med/start_med > 1.8:
            good_events.append(event)

    df = df[df['event'].isin(good_events)]

    return df



# ----- fit functions for the profiling --------
def floating_gaus(x, a, mu, sigma, ped=0.0):
    return ped + a*gaus(x,mu,sigma)

def landau(x, a, mu, sigma):
    #return ped + a*moyal(x,mu,sigma)     
    return a*moyal(x,mu,sigma)                                         

def both(x, a, mu_landau, sigma_landau, sigma_gaus=0.1):
    lnd = landau(x, a, mu_landau, sigma_landau)
    return gaussian_filter1d(lnd,sigma_gaus)

def moyal(x, mu=0, sigma=1):
    y = (x - mu)/sigma
    return np.exp(-(y + np.exp(-y))/2)/(sigma*np.sqrt(2*math.pi))

def gaus(x, mu=0, sigma=1):
    return (1/(sigma*np.sqrt(2*math.pi)))*np.exp(-0.5*((x-mu)/sigma)**2)
# ----------------------------------------------



# ------ recombination models to fit -----------


def birks(p, x):
    A_b = p[0]
    k_b = p[1]
    return (A_b/W)*(x/(1+(k_b/(eps*rho))*x))

def birks_wrapper(x, p0, p1):
    return birks(p=[p0,p1],x=x)

def modified_box(p, x):
    alpha = p[0]
    beta_etc = p[1]/(rho*eps)
    return (np.log((x)*beta_etc + alpha))/(beta_etc*W)

def modified_box_wrapper(x, p0, p1):
    return modified_box(p=[p0,p1],x=x)

def output_stats(x, y, sigma_y, model, result, cov, label):
    print(f"\nResults for {label} model: ")
    p0 = result[0]
    p0_error = np.sqrt(cov[0][0])
    p1 = result[1]
    p1_error = np.sqrt(cov[1][1])
    chi2 = np.sum(np.square((y - model(x, p0, p1))/sigma_y))
    ndof = len(x)
    print(f"  best fit parameters p0: {p0:.3f} +/-{p0_error:.3f}, p1: {p1:.3f} +/-{p1_error:.3f}") #, p1/rho*eps = {p1/(rho*eps):.3f}")
    print(f"  chi2/ndof = {chi2/ndof:.2f}")


def output_stats_new(result, sigma, chi2_ndof):
    p0 = result[0]
    p0_error = sigma[0]
    p1 = result[1]
    p1_error = sigma[1]
    print(f"  best fit parameters p0: {p0:.5f} +/-{p0_error:.3f}, p1: {p1:.5f} +/-{p1_error:.3f}") #, p1/rho*eps = {p1/(rho*eps):.3f}")
    print(f"  chi2/ndof = {chi2_ndof:.2f}")


# ----------------------------------------------





# --------- plotting helper functions -------------

def plot_2D(dEdx,dQdx,title='dQdx vs dEdx'):
    plt.figure(figsize=(10,6))
    plt.hist2d(dEdx,dQdx,bins=[np.arange(0.0,30.0,0.1),np.arange(0.0,400000.0,2000.0)],norm=LogNorm(vmin=1.0,vmax=2000.0),cmap='YlGnBu_r')    
    plt.colorbar()
    plt.xlabel("dE/dx [MeV/cm]",labelpad=10,fontsize=16)
    plt.ylabel("dQ/dx [e/cm]",labelpad=10,fontsize=16)
    plt.title(title,pad=20)

def plot_fits(dE_bin_centers,dQ_fitted,dQ_sigma,dE_sigma,birks_p0,birks_p1,box_p0,box_p1,label_1,label_2,title):
    plt.figure(figsize=(9,6))
    plt.errorbar(dE_bin_centers,dQ_fitted,yerr=dQ_sigma,xerr=dE_sigma,color='coral',linestyle='None',label="Fitted dQ/dx",linewidth=2)
    plt.plot(dE_bin_centers,birks_wrapper(dE_bin_centers,birks_p0,birks_p1),color='deeppink',label=label_1,linewidth=2)
    plt.plot(dE_bin_centers,modified_box_wrapper(dE_bin_centers,box_p0,box_p1),color='orange',label=label_2,linestyle='--',linewidth=2)
    plt.xlabel("dE/dx [MeV/cm]",labelpad=10,fontsize=16)
    plt.ylabel("dQ/dx [e/cm]",labelpad=10,fontsize=16)
    plt.legend(prop={'size':16})
    plt.title(title)

def overlay_fits(dEdx,dQdx,dE_bin_centers,dQ_fitted,dQ_sigma,dE_sigma,birks_p0,birks_p1,box_p0,box_p1,label_1,label_2,title):
    plot_2D(dEdx,dQdx,title)
    plt.errorbar(dE_bin_centers,dQ_fitted,yerr=dQ_sigma,xerr=dE_sigma,color='coral',linestyle='None',label="Fitted dQ/dx",linewidth=2)
    plt.plot(dE_bin_centers,birks_wrapper(dE_bin_centers,birks_p0,birks_p1),color='deeppink',label=label_1,linewidth=2)
    plt.plot(dE_bin_centers,modified_box_wrapper(dE_bin_centers,box_p0,box_p1),color='orange',label=label_2,linestyle='--',linewidth=2)
    plt.legend(prop={'size':16})

# -------------------------------------------------



# ------------ main part of the code ----------------

# read in the data frame (with all events)
df_all = get_df_from_file(USE_MC)

print(df_all.keys())

# sometimes we might want to look at the uncalibrated dq/dx values
dqdx_string = 'cali_dqdx'
if USE_UNCALIBRATED_DQDX:
    dqdx_string = 'dqdx'

# plot dQ/dx vs dE/dx pre cleaning
plot_2D(df_all['dedx_hyp'],df_all[dqdx_string],'dQdx vs dEdx (pre cleaning)')
plt.savefig(f"pre_cleaning_{TAIL}.png")

# apply quality cuts (e.g. on dx)
df = apply_quality_cuts(df_all)

# plot dQ/dx vs dE/dx post cleaning
plot_2D(df['dedx_hyp'],df[dqdx_string],'dQdx vs dEdx (post cleaning)')
plt.savefig(f"post_cleaning_{TAIL}.png")

# --------- bin dQ/dx in dE/dx and 1D fit -----------

# get dE/dx bins
dE_bin_edges, dE_bin_centers, dE_sigma_stat, dE_sigma_sys = fh.get_variable_dE_bins(SIGMA_RR,MIN_BIN_WIDTH)

if USE_FIXED_BINS:
    dE_bin_edges, dE_bin_centers, dE_sigma = fh.get_dE_bins(NBINS)



# fit dQ/dx in all of the dE/dx bins
profile_fit_function = both # floating_gaus # landau
profile_fit_function_hs = both # floating_gaus # landau
dQ_fitted, dQ_sigma, dQ_fit_okay, dQ_ns = fh.fit_all_dQ(df,dqdx_string,dE_bin_edges,dE_bin_centers,profile_fit_function,profile_fit_function_hs,PLOT_PROFILES) 

# update the dEdx stat systematic to be error on mean not rms
dE_sigma_stat = dE_sigma_stat/np.sqrt(dQ_ns)

# combine the dEdx systematic
# NOTE: basically no stat contribution!
dE_sigma = np.sqrt(dE_sigma_stat**2 + dE_sigma_sys**2)
#dE_sigma = dE_sigma_sys

# add the systematic onto the dQ/dx uncertainties
#print(dQ_sigma, SIGMA_DQDX*dQ_fitted)
#dQ_sigma = np.sqrt(dQ_sigma**2 + (SIGMA_DQDX*dQ_fitted)**2)

#print(f"stat = {dQ_sigma},  sys = {SIGMA_DQDX*dQ_fitted}")
#dQ_sigma = np.sqrt(dQ_sigma**2 + (SIGMA_DQDX*dQ_fitted)**2)



# remove any failed fits
if REMOVE_FAILED_FITS:
    dQ_fitted = dQ_fitted[dQ_fit_okay]
    dQ_sigma = dQ_sigma[dQ_fit_okay]
    dE_bin_centers = dE_bin_centers[dQ_fit_okay]
    dE_sigma = dE_sigma[dQ_fit_okay]


# make the covariance matrix
sigma2_dqdx_diag_array =  dQ_sigma**2 + sigma2_dqdx_tot*(dQ_fitted**2)
cov_dQ, sigma_dQ_diag = make_covariance(sigma2_dqdx_corr, sigma2_dqdx_diag_array, dQ_fitted)
cov_dE = make_covariance_dE(dE_sigma)
print(dE_sigma)
print(cov_dE)
print("*******************")


# ------------- 2D fit the dQ/dx vs dE/dx -------------------

birks_p0 = 0.8 # where were these from?
birks_p1 = 0.0486 # ?
box_p0 = 0.93 # argoneut                                      
box_p1 = 0.212 # argoneut  
model_birks = od.Model(birks)
model_box = od.Model(modified_box)

# plot the starting values of dQ/dx vs dE/dx and the input values
plot_fits(dE_bin_centers,dQ_fitted,sigma_dQ_diag,dE_sigma,birks_p0,birks_p1,box_p0,box_p1,"Birks Starting Values","Box Starting Values","Starting Parameters for the fits")
plt.savefig(f"fit_inputs_{TAIL}.pdf")
plt.savefig(f"fit_inputs_{TAIL}.png")


# do the 2D fit to Birks and Box
#data = od.RealData(dE_bin_centers,dQ_fitted,dE_sigma,dQ_sigma)
#data = od.RealData(dE_bin_centers,dQ_fitted,sx=dE_sigma,covy=cov_dQ)

#odr_birks = od.ODR(data, model_birks, beta0=[birks_p0,birks_p1])
#output_birks = odr_birks.run()
#result_birks, cov_birks = output_birks.beta, output_birks.cov_beta
#result_birks, cov_birks = curve_fit(birks_wrapper, dE_bin_centers, dQ_fitted, p0=[birks_p0, birks_p1], sigma=cov_dQ, absolute_sigma=True)
result_birks, sigma_birks, cov_birks, chi2_ndof_birks = xef.fit(birks, dE_bin_centers, np.identity(dE_bin_centers.shape[0])*dE_sigma**2, dQ_fitted, cov_dQ, p0=np.array([birks_p0, birks_p1]))

#odr_box = od.ODR(data, model_box, beta0=[box_p0,box_p1])
#output_box = odr_box.run()
#result_box, cov_box = output_box.beta, output_box.cov_beta
#result_box, cov_box = curve_fit(modified_box_wrapper, dE_bin_centers, dQ_fitted, p0=[box_p0, box_p1], sigma=cov_dQ, absolute_sigma=True)
result_box, sigma_box, cov_box, chi2_ndof_box = xef.fit(modified_box, dE_bin_centers, np.identity(dE_bin_centers.shape[0])*dE_sigma**2, dQ_fitted, cov_dQ, p0=np.array([box_p0, box_p1]))

#output_stats(dE_bin_centers, dQ_fitted, sigma_dQ_diag, modified_box_wrapper, result_box, cov_box, 'Modified Box')
output_stats_new(result_box,sigma_box,chi2_ndof_box)
print(f"ArgoNeuT values (Box): p0 = {box_p0}, p1 = {box_p1}")

#output_stats(dE_bin_centers, dQ_fitted, sigma_dQ_diag, birks_wrapper, result_birks, cov_birks, 'Birks')
output_stats_new(result_birks,sigma_birks,chi2_ndof_birks)
print(f"Icarus (?) values (Birks): p0 = {birks_p0}, p1 = {birks_p1}")

# plot the best fit dQ/dx vs dE/dx overlaid on input values
plot_fits(dE_bin_centers,dQ_fitted,sigma_dQ_diag,dE_sigma,result_birks[0],result_birks[1],result_box[0],result_box[1],"Birks Best Fit","Modified Box Best Fit","")
plt.savefig(f"fit_results_{TAIL}.pdf")

# plot the best fit dQ/dx vs dE/dx overlaid on input data
overlay_fits(df['dedx_hyp'],df[dqdx_string],dE_bin_centers,dQ_fitted,sigma_dQ_diag,dE_sigma,result_birks[0],result_birks[1],result_box[0],result_box[1],"Birks Best Fit","Modified Box Best Fit","")
plt.ylim((50000,350000))
plt.xlim((2.0,22.0))
plt.savefig(f"fit_results_overlay_{TAIL}.png")

# show the residual ranges (remember typical pitch = 0.6 cm with some variation (+/-0.1?))
plt.figure()
plt.hist(df['resrange'],bins=np.arange(0.0,3.0,0.01),color='mediumvioletred',histtype='step',label='Hit Distribution (Pandora)')
distance = 0.6
xs = np.arange(0.2,3.0,0.01)
cv = converter()
plt.plot(xs,cv.get_dEdx(xs),color='purple',label="dE/dx")
rainbow = ['salmon','orange','gold','lightgreen','darkturquoise']
i = 0
for x in np.arange(0.25,5.0*distance,distance):
    plt.axvline(x=x,color=rainbow[i],label=f"{x:.2f} cm, dE/dx = {cv.get_dEdx(x):.0f} MeV/cm")
    i += 1
plt.legend()
plt.xlabel("Residual Range [cm]")
plt.ylabel("Counts [per 0.01 cm] or dE/dx [Mev/cm], enjoy!")
plt.ylim((0.0,300.0))
plt.xlim((0.0,3.0))
plt.savefig(f"rainbow_{TAIL}.png")

# overlay
plot_2D(df['dedx_hyp'],df[dqdx_string],'dQdx vs dEdx (post cleaning)')
i = 1
for x in np.arange(0.25+distance,5.0*distance,distance):
    plt.axvline(x=cv.get_dEdx(x),color=rainbow[i],label=f"{x:.2f} cm, dE/dx = {cv.get_dEdx(x):.0f} MeV/cm")
    i += 1
plt.legend()
plt.ylim((0.0,5e5))
plt.text(21.0,1e5,"final-1 hits",color=rainbow[1])
plt.text(16.0,5e4,"final-2 hits",color=rainbow[2])
plt.savefig(f"rainbow_overlay_{TAIL}.png")

# show the pitch
plt.figure()
print(f"Average pitch = {np.mean(df['pitch']):.2f} cm")
plt.hist(df['pitch'],bins=np.arange(0.5,0.7,0.001),histtype='step',label='pitch')
plt.legend()
plt.xlabel("Pitch [cm]")
plt.savefig(f"pitch_{TAIL}.pdf")

plt.show()







