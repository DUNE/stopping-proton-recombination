import numpy as np
from scipy import interpolate

class converter:
    # initialize
    def __init__(self):
        self.data = np.load("../proton_energy_loss/energy_loss_vs_rr.npz")
        self.fn = interpolate.interp1d(self.data['rr'],self.data['dEdx'],kind='cubic',bounds_error=False,fill_value=-1.0)
        self.inv_fn = interpolate.interp1d(self.data['dEdx'],self.data['rr'],kind='cubic',bounds_error=False,fill_value=-1.0)


    # convert residual range to dE/dx (one to one)
    def get_dEdx(self, rr):
        return self.fn(rr)


    def get_rr(self, dedx):
        return self.inv_fn(dedx)
