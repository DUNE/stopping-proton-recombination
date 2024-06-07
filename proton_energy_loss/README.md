Run on the lx machines (or change the path to link to Michael Mooney's root file)

$> python read_values.py

will make a .npz archive you can then read in and interpolate to convert residual range into dEdx.  Currently one to one - would there be a benefit in making it one to many?  It really depends on the binning of dEdx due to other systematics in the fit.

An example of using the conversion with interpolation can be found in test_read_in.py (also outputs some comparisons with the ArgoNeuT approximation for sanity checking)

Note correct use of histograms is to use the bin centers in both dEdx and residual range.