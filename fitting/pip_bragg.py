import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d


USE_MC = False
prefix = 'data'
DEBUG_PLOTS = False

if USE_MC:
    prefix = 'mc'

# read in the bragg data
bragg = pd.read_csv("../bragg_data/apdata_lookup.txt"," ",index_col=False,names=['energy','energy_csda','resrange'])
LAr_density = 1.784e-3
bragg['resrange'] = bragg['resrange']/LAr_density
print(bragg['resrange'])
#bragg['energy'] = bragg['energy']/LAr_density
#bragg['energy_csda'] = bragg['energy_csda']/LAr_density
dedx = np.array([(bragg['energy'][i] - bragg['energy'][i-1])/(bragg['resrange'][i] - bragg['resrange'][i-1]) for i in range(1,len(bragg['energy']))]) 
dedx_csda = np.array([(bragg['energy_csda'][i] - bragg['energy_csda'][i-1])/(bragg['resrange'][i] - bragg['resrange'][i-1]) for i in range(1,len(bragg['energy_csda']))]) 
bragg_resrange_c = np.array([0.5*(bragg['resrange'][i] + bragg['resrange'][i-1])for i in range(1,len(bragg['resrange']))])


plt.figure()
plt.plot(bragg['resrange'][:-1],dedx)
plt.plot(bragg['resrange'][:-1],dedx_csda)


# Helper
def get_shift(x1,y1,z1,x2,y2,z2):
    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

# Floating Bragg Peak
bragg_func = interp1d(bragg['resrange'][:-1],dedx,kind='cubic',bounds_error=False,fill_value=0.0)
#bragg_func = interp1d(bragg_resrange_c,dedx,kind='cubic',bounds_error=False,fill_value=0.0)

def floating_bragg(x, p0, p1, p2):
    return p0*500.0*bragg_func(x + p1) + p2


#def floating_bragg(x, p1):
#    return 35.0*bragg_func(x + p1) + 1.0




# Bragg fitting function (just translated from the C++)
def extract_track_shift(resrange,dEdx):
    start_params = np.array([0.1,0.0,0.0])
    if USE_MC:
        start_params = np.array([10.0,0.0,0.0])    
    params, errors = start_params, np.zeros((3,3))
    #start_params = np.array([0.0])
    #params, errors = start_params, np.zeros((1,1))
    fit_status = 0
    
    # crop the low hits off the end of the track (for fitting)
    max_dedx = np.argmax(dEdx)
    #print(dEdx[0],dEdx[1],resrange[0])
    if USE_MC:
        dEdx = dEdx[:max_dedx+1]
        resrange = resrange[:max_dedx+1]
    else:
        # in data there are more spikes
        max_dedx = np.argmax(dEdx[resrange < 5.0])
        # also tracks are flipped
        dEdx = dEdx[max_dedx:]
        resrange = resrange[max_dedx:]

    #print(dEdx[0],resrange[0])

    try:
        params, errors = curve_fit(floating_bragg,resrange,dEdx,p0=start_params)
        fit_status = 1
        # also check the fit errors are not large?
    except:
        print('fitting failed!')
        fit_status = 0
    trackshift = params[1] # falls back to 0.0 if the fit fails
    return trackshift, params, fit_status



# Read in the data to shift
df = uproot.open(f"../output/{prefix}_proton_recombination.root")['output_tree'].pandas.df()
print(df.keys())

# Add the new endpoint column
df['resrange_bragg'] = np.ones_like(df['resrange'])
df['trackshift'] = np.zeros_like(df['resrange'])
df['bragg_okay'] = np.zeros_like(df['resrange'])
nplots = 0
for event in set(df['event']):
    df_e = df[df['event'] == event]
    trackshift, fit_params, fit_status = extract_track_shift(df_e['resrange'],df_e['dedx']) # think about it, is this the right dedx?
    #print(event, trackshift)
    df.loc[df['event'] == event,'resrange_bragg'] = df_e['resrange'] + trackshift # plus or minus?
    df.loc[df['event'] == event,'trackshift'] = trackshift
    df.loc[df['event'] == event,'bragg_okay'] = fit_status
    if nplots < 10 and DEBUG_PLOTS:
        plt.figure()
        if USE_MC:
            plt.scatter(df_e['hit_z'],df_e['hit_y'],label='reco hits')
            plt.scatter(df[df['event']==event]['true_endpos.fZ'],df[df['event']==event]['true_endpos.fY'],label='true end pos')
            plt.scatter(df[df['event']==event]['reco_endpos.fZ'],df[df['event']==event]['reco_endpos.fY'],label='reco end pos')
            plt.legend()
        plt.figure()
        plt.plot(df_e['resrange'],df_e['dedx'],label='dedx')        
        plt.plot(df_e['resrange'],df_e['dedx_hyp'],label='dedx_hyp')
        plt.plot(df_e['resrange'],floating_bragg(df_e['resrange'],fit_params[0],fit_params[1],fit_params[2]),label='fitted bragg')
        #plt.plot(df_e['resrange'],floating_bragg(df_e['resrange'],fit_params[0]),label='fitted bragg')
        plt.legend()
        plt.figure()
        plt.plot(df_e['resrange'],df_e['cali_dqdx'],label='cali_dqdx')
        plt.legend()
        nplots += 1

# analyse results
print(df['resrange_bragg'],df['resrange'],df['dedx'],df['bragg_okay'])
plt.figure()
plt.hist(df['trackshift'],bins=np.arange(-5.0,15.0,0.2))
plt.xlabel('trackshift (up track?) [mm]')


if USE_MC:
    LArsoft_shifts = get_shift(df['true_endpos.fX'],df['true_endpos.fY'],df['true_endpos.fZ'],df['reco_endpos.fX'],df['reco_endpos.fY'],df['reco_endpos.fZ'])
    #bragg_shifts = get_shift(df['true_endpos.fX'],df['true_endpos.fY'],df['true_endpos.fZ'],df['reco_endpos.fX'],df['reco_endpos.fY'],df['reco_endpos.fZ'])
    print(df['reco_endpos.fX'],df['true_endpos.fZ'])
    print(LArsoft_shifts, df['trackshift'])
    plt.figure()
    #plt.hist(LArsoft_shifts,bins=np.arange(-10.0,10.0,0.1),histtype='step',label='LArsoft')

    plt.hist(df['true_endpos.fZ'] - df['reco_endpos.fZ'],bins=np.arange(-10.0,10.0,0.1),histtype='step',label='LArsoft')
    # 0 = 20, A = 80 => H = 83, costheta = A/H = 0.96
    shift = df['true_endpos.fZ'] - df['reco_endpos.fZ']
    shift_bragg = df['true_endpos.fZ'] - (df['reco_endpos.fZ'] - 0.96*df['trackshift'])
    plt.hist(shift_bragg,bins=np.arange(-10.0,10.0,0.1),histtype='step',label='Bragg')
    #plt.hist(bragg_shifts,bins=np.arange(-10.0,10.0,0.1),histtype='step',label='Bragg')
    df['shift'] = shift
    df['shift_bragg'] = shift_bragg
    plt.xlabel('True - Reco (Z component)')
    plt.legend()

    #trunc = LArsoft_shifts[LArsoft_shifts < 2.0]
    trunc = df['shift_bragg']
    trunc = trunc[trunc < 1.5]
    trunc = trunc[trunc > -1.5]
    systematic = np.mean(trunc)
    ssigma = np.sqrt(np.mean((trunc - systematic)**2))
    print(systematic, ssigma)



# Save the dataframe
df.to_csv(f'../output/{prefix}_with_bragg_fit.csv',encoding='utf-8')

plt.show()
