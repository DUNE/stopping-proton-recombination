import numpy as np
import uproot3 as uproot
import math
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LogNorm


base_path = "/vols/dune/awaldron/protoDUNE/sungbin_endpoint/"
file_extensions = ["data_proton_recombination_0_1mm.root","mc_proton_recombination_0_1mm_full_stats.root"]
labels = ["Data [0.1mm]","MC [0.1mm]"]
bins = np.arange(-20.1,20.1,0.2)
#bins = np.arange(-100.1,100.1,0.2)

def get_df_from_file(filename):
    tree = uproot.open(filename)['output_tree']
    df = tree.pandas.df()
    # cut on high chi2, suspected pions
    #df = df[df['chi2_sungbin'] < 10.0]
    return df



df = get_df_from_file(base_path + file_extensions[0])
plt.figure()
plt.hist(df['chi2_sungbin'], bins=np.arange(0.0,100.0,1.0))



plt.figure()

for i, f in enumerate(file_extensions):
    df = get_df_from_file(base_path + f)
    df = df.loc[(slice(None),0), :] # get only the first hit for each track
    shifts = df['rr_shift_sungbin']
    plt.hist(shifts, bins=bins, density=True, histtype='step', label=labels[i])

plt.xlabel("Endpoint Offset [cm]")
plt.ylabel("Density")

plt.legend()
plt.savefig("endpoint_offsets.pdf")
plt.savefig("endpoing_offsets.png")

plt.show()







