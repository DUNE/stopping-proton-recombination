import numpy as np
import uproot3 as uproot
import math
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LogNorm


base_path = "/vols/dune/awaldron/protoDUNE/sungbin_endpoint/"
#file_extensions = ["mc_proton_recombination_10mm.root","mc_proton_recombination_5mm.root","mc_proton_recombination_4mm.root","mc_proton_recombination_3mm.root","mc_proton_recombination_2mm.root","mc_proton_recombination_1mm.root"]
#file_extensions = ["mc_proton_recombination_4mm.root","mc_proton_recombination_3mm.root","mc_proton_recombination_2mm.root","mc_proton_recombination_1mm.root"]
file_extensions = ["mc_proton_recombination_5mm.root","mc_proton_recombination_4mm.root","mc_proton_recombination_3mm_full_stats.root","mc_proton_recombination_2mm_full_stats.root","mc_proton_recombination_1mm.root"]
labels = ["MC [5mm]","MC [4mm]","MC [3mm]","MC [2mm]","MC [1mm]"]
#labels = ["MC [5mm]","MC [1.5mm]","MC [1mm]","MC [0.5mm]"]
#bins = np.arange(-5.1,5.1,0.02)
bins = np.arange(-5.1,5.1,0.05)


def get_df_from_file(filename):
    tree = uproot.open(filename)['output_tree']
    df = tree.pandas.df()
    return df


def get_rms(stuff):
    # only want the peak region
    n = stuff.shape[0]
    stuff = stuff[np.abs(stuff) < 2.0]
    mu = np.mean(stuff)
    sigma = np.sqrt(np.mean((stuff-mu)**2))
    sigma_sigma = sigma/(np.sqrt(2*(n-1)))
    return sigma, sigma_sigma


plt.figure()

for i, f in enumerate(file_extensions):
    df = get_df_from_file(base_path + f)
    df = df.loc[(slice(None),0), :] # get only the first hit for each track
    shifts = df['rr_shift_sungbin']
    rms, sigma_rms = get_rms(shifts)
    this_label = labels[i] + f" rms = {rms:.3f} +/- {sigma_rms:.3f}"
    plt.hist(shifts, bins=bins, density=True, linewidth=2, histtype='step', label=this_label)
plt.xlabel("Endpoint Shift [cm]")
plt.ylabel("Density")
plt.ylim(0.0,2.5)
plt.legend()

plt.savefig("endpoint_resolution.pdf")
plt.savefig("endpoint_resolution.png")

plt.show()







