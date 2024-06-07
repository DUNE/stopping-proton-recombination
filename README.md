# stopping-protons

Contains code to extract reconstructed stopping proton-like events from the protoDUNE monte carlo (ProtonRecombination) and data (ProtonRecombination_run5387).  Note that there are some hard coded paths to data in there so it needs running on the dunegpvms.  The code currently only saves dQ/dx vs R, dE/dx vs R and dQ/dx vs dE/dx.  To run:

```
    $> root -l  
    root> .L BetheBloch.cxx+
    root> .L ProtonRecombination.C+  
    root> ProtonRecombination t;  
    root> t.Loop();  
```


Before running `ProtonRecombination`  or the data version, check it is using the right version of the calibration constants by checking which version of the `dedx_function` header it is using.  

In case calibration constants need to be regenerated, on the gpvms do `setup_protodune_calib` and follow the instructions here:  https://github.com/sungbinoh/protoduneana/blob/Calibration_constant/README.md#instruction-to-run-protodune-sp-calibration-codes.  Note that you only need to get the X and YZ maps, and copy the files back over (and the plane 2 normalisation: norm = 59.29/number from calibration).  You don't need to do the last dedx correction part.

Then make sure to rename the output root file something sensible and use that root file in `fitting/fit_data.py`.  Also check the options in this file are sensible.  Finally update the parameters you found in `fitting/all_sys.py` and run to get the final results for data with full systematic uncertainties.