import numpy as np
import uproot3 as uproot
import math
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LogNorm
from rr_to_dEdx_converter import converter

def get_df_from_file():
    filename = "/vols/dune/awaldron/protons/root_files/mc_proton_recombination_xyz.root"
    tree = uproot.open(filename)['output_tree']
    df = tree.pandas.df()
    conv = converter()
    #df['dedx_hyp'] = conv.get_dEdx(df['resrange_new'])
    df['dedx_hyp'] = conv.get_dEdx(df['resrange'])
    return df


def apply_quality_cuts(df):
    # cut on end of tracks
    df = df[df['resrange'] > 1.0]

    # cut on angle - why does this work?
    df = df[df['pitch'] < 0.64]
  
    #df = df[df['rr_shift_sungbin'] < 3.0]

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





# ------------ main part of the code ----------------

# read in the data frame (with all events)
df_all = get_df_from_file()
print(df_all.keys())
print(df_all['x'])

# apply quality cuts (e.g. on dx)
df = apply_quality_cuts(df_all)


print(df['x'])

mu = np.mean(df['x'])
print(mu)
print(np.min(df['x']), np.max(df['x']))


sigma = np.sqrt(np.mean((df['x']-mu)**2))
print(sigma)
print(mu - sigma, mu + sigma)


x_bins = np.arange(-100.0,0.0,1.0)

plt.figure()
plt.hist(df_all['x'],bins=x_bins,histtype='step',label="pre cleaning")
plt.hist(df['x'],bins=x_bins,histtype='step',label="post cleaning")
plt.legend()


plt.show()







