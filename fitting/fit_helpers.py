import numpy as np
import uproot
import math
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from rr_to_dEdx_converter import converter


def get_dE_bins(nbins):
    dE_min, dE_max = 0.0, 30.0
    bin_width = (dE_max - dE_min)/float(nbins)
    bin_edges = np.arange(dE_min,dE_max+bin_width,bin_width)
    bin_centers = np.arange(dE_min+0.5*bin_width,dE_max+0.5*bin_width,bin_width)
    return bin_edges, bin_centers, 0.5*bin_width*np.ones_like(bin_centers)


def dedx_hyp(rr):
    return 17.0*rr**(-0.42)


def get_variable_dE_bins(sigma,min_bin_size=0.2):
    # note: sigma is in resrange (e.g. wire spacing)
    #  however min_bin_size is in dE/dx (to not get too low stats)
    #  be careful of the logic in here, the ordering flips with the
    #  dE/dx transformation.

    # full range of the bins
    offset = 1.0
    #rr_min, rr_max = sigma+offset, 120.0 + sigma
    #rr_min, rr_max = offset, 120.0 + sigma
    rr_min, rr_max = offset, 95.0 + sigma


    conv = converter()

    # variable width bins (full range)
    rr_bin_edges = np.arange(rr_min,rr_max,2.0*sigma)
    #print(rr_bin_edges)
    var_bin_edges = conv.get_dEdx(rr_bin_edges)
    #print(var_bin_edges)
    var_bin_widths = np.array([var_bin_edges[i-1]-var_bin_edges[i] for i in range(1,var_bin_edges.shape[0])]) 
    

    # find out which are greater than min size and where they stop
    bins_to_switch = var_bin_widths[var_bin_widths > min_bin_size].shape[0]
    bin_edges = var_bin_edges[:bins_to_switch+1]

    # fixed bin width part
    fixed_start = bin_edges[bin_edges.shape[0]-1]
    n_fixed_bins = math.ceil((fixed_start - conv.get_dEdx(rr_max))/min_bin_size)
    fixed_end = fixed_start - min_bin_size*n_fixed_bins
    fixed_bin_edges = np.arange(fixed_start,fixed_end,-min_bin_size)

    # now choose variable bin width *unless* smaller than min bin width
    bin_edges = np.append(bin_edges,fixed_bin_edges)
    bin_widths = np.array([bin_edges[i-1]-bin_edges[i] for i in range(1,bin_edges.shape[0])])
    bin_centers = bin_edges[:-1] - 0.5*bin_widths


    # errors - approximate undertainty at bin center by looking at 
    #  the dedx at bin center and bin center - bin width
    bin_sigmas_stat = np.abs(bin_centers - conv.get_dEdx(conv.get_rr(bin_centers-0.28866*bin_widths)))
    bin_sigmas_sys = np.abs(bin_centers - conv.get_dEdx(conv.get_rr(bin_centers)-sigma))
    #bin_sigmas = np.sqrt(bin_sigmas_stat**2 + bin_sigmas_sys**2)
    #bin_sigmas = 0.5*bin_widths

    return np.flip(bin_edges), np.flip(bin_centers), np.flip(bin_sigmas_stat), np.flip(bin_sigmas_sys)



def fit_dQ(dQs, dQ_bin_centers, dQ_bin_values, fit_function, fit_function_hs, plot_profiles=False):
    # we'll only fit the peak region
    max_val = np.max(dQ_bin_values)
    threshold = 0.25*max_val
    if max_val > 1000.0:
        fit_function = fit_function_hs
    dQ_bin_centers_t = dQ_bin_centers[dQ_bin_values > threshold]
    dQ_bin_values_t = dQ_bin_values[dQ_bin_values > threshold]

    mean = np.average(dQ_bin_centers_t,weights=dQ_bin_values_t) #np.mean(dQs)
    rms = 10000 #np.sqrt(np.sum(np.square(dQs - mean)))
    a0 = max_val*rms*2.5

    mpv = mean
    mpv_error = -1
    fit_okay = False

    try:
        params, errors = curve_fit(fit_function,dQ_bin_centers_t,dQ_bin_values_t,p0=[a0,mean,rms])
        mpv = params[1]
        mpv_error = np.sqrt(errors[1][1])
        fit_okay = True

        if plot_profiles:
            plt.figure()
            plt.step(dQ_bin_centers, dQ_bin_values, where='mid')
            if params.shape[0] > 3:
                initial_fit = fit_function(dQ_bin_centers_t, a0, mean, rms, 0.1)
                plt.plot(dQ_bin_centers_t, initial_fit, label="pre fit")
                plt.plot(dQ_bin_centers_t, fit_function(dQ_bin_centers_t, params[0], params[1], params[2], params[3]), label="post fit")
            else:
                initial_fit = fit_function(dQ_bin_centers_t, a0, mean, rms)
                plt.plot(dQ_bin_centers_t, initial_fit, label="pre fit")
                plt.plot(dQ_bin_centers_t, fit_function(dQ_bin_centers_t, params[0], params[1], params[2]), label="post fit")
            plt.legend()
    except:
        print("fit failed, falling back!")

    return mpv, mpv_error, fit_okay



def fit_all_dQ(df,dqdx_string,bin_edges,bin_centers,fit_function,fit_function_hs,plot_profiles):
    plot_profiles_input = plot_profiles
    output, sigma = np.ones_like(bin_centers), np.ones_like(bin_centers)
    ns = np.zeros_like(bin_centers)
    fit_okay = np.array([False]*bin_centers.shape[0])
    # loop over dE bins (we have to fit each anyway)
    for i in range(bin_centers.shape[0]):

        # bin calibrated dQ in slices of dE        
        lb, ub = bin_edges[i], bin_edges[i+1]
        dQs = df[(df['dedx_hyp'] < ub) & (df['dedx_hyp'] > lb)][dqdx_string]
        # if there are not enough entries, skip
        #if dQs.shape[0] < 320:
         #   print(dQs.shape[0])
          #  continue
        if dQs.shape[0] < 100:
            continue

        # bin all the dQ in this dE bin in dQ bins
        lb_dQ, ub_dQ = 0.0, 400000.0
        nbins = 200
        bin_width_dQ = (ub_dQ - lb_dQ)/float(nbins)
        bin_centers_dQ = np.arange(lb_dQ + 0.5*bin_width_dQ,ub_dQ + 0.5*bin_width_dQ,bin_width_dQ)
        dQ_binned, _ = np.histogram(dQs,bins=nbins,range=[lb_dQ,ub_dQ])

        output[i], sigma[i], fit_okay[i] = fit_dQ(dQs,bin_centers_dQ,dQ_binned,fit_function,fit_function_hs,plot_profiles)
        ns[i] = dQs.shape[0]


        # if sigma is huge deem it failed
        if sigma[i] > 10*output[i]:
            fit_okay[i] = False

        # systematics for dQ/dx done elsewhere (since different for data/MC)

    return output, sigma, fit_okay, ns







