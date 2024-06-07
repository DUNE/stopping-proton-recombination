import numpy as np
import gvar as gv
import lsqfit


def get_values(var):
    values = var['p']
    return np.array([v for v in values])


# pass in birks or modified_box, not the wrapper
def fit(f,x,x_cov,y,y_cov,p0):
    x_data = make_data(x,x_cov)
    y_data = make_data(y,y_cov)
    prior = make_prior(x_data, p0)

    fit_function = lambda p: f(p['p'],p['x'])

    fit = lsqfit.nonlinear_fit(prior=prior, data=y_data, fcn=fit_function)
    return get_values(fit.pmean), get_values(fit.psdev), fit.cov, fit.chi2/fit.dof


def make_data(input_data, input_data_cov):
    return gv.gvar(input_data, input_data_cov)


def make_prior(x, p):
    prior = gv.BufferDict()
    n = p.shape[0]
    prior['p'] = gv.gvar(p,np.identity(n)*p*0.1) # does this affect results?
    prior['x'] = x
    return prior
