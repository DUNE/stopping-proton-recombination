#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <string>
#include <TImage.h>
#include <iomanip>
#include <TSpline.h>
#include <TText.h>
#include <TFrame.h>
#include <TH3.h>
#include <TH2.h>
#include <TFile.h>
#include "TMinuit.h"


double calib_factor_sce_off =5.02e-3;
double normalisation_factor_sce_off=0.984;//for plane 2
TFile my_file0_sce_off("./cali_off/YZ_RITsceoff.root");//YZ correction factors
TFile my_file2_sce_off("./cali_off/X_RITsceoff.root"); //X correction factors


double Rho_sce_off = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
double betap_sce_off = 0.212;//(kV/cm)(g/cm^2)/MeV
double alpha_sce_off = 0.93;//parameter from ArgoNeuT experiment at 0.481kV/cm 
double Wion_sce_off = 23.6e-6;//parameter from ArgoNeuT experiment at 0.481kV/cm
double Efield_sce_off = 0.4867;//kV/cm protoDUNE electric filed



TH1F *X_correction_hist_sce_off = (TH1F*)my_file2_sce_off.Get("dqdx_X_correction_hist_2");
TH2F *YZ_correction_neg_hist_2_sce_off=(TH2F*)my_file0_sce_off.Get("correction_dqdx_ZvsY_negativeX_hist_2");
TH2F *YZ_correction_pos_hist_2_sce_off=(TH2F*)my_file0_sce_off.Get("correction_dqdx_ZvsY_positiveX_hist_2");


Double_t dedx_function_35ms_sce_off(double dqdx, double x, double y, double z){
  double Ef=0.4867;
  double Cx=X_correction_hist_sce_off->GetBinContent(X_correction_hist_sce_off->FindBin(x));
  double Cyz=0.0;
  if(x<0){
    Cyz=YZ_correction_neg_hist_2_sce_off->GetBinContent(YZ_correction_neg_hist_2_sce_off->FindBin(z,y));
  }
  if(x>=0){
    Cyz=YZ_correction_pos_hist_2_sce_off->GetBinContent(YZ_correction_pos_hist_2_sce_off->FindBin(z,y));
  }
 double corrected_dqdx=dqdx*Cx*Cyz*normalisation_factor/calib_factor;
 return (exp(corrected_dqdx*(betap_sce_off/(Rho_sce_off*Ef)*Wion))-alpha_sce_off)/(betap_sce_off/(Rho_sce_off*Ef));
}



Double_t dqdx_uniform_corr_sce_off(double dqdx, double x, double y, double z){
  double Ef=0.4867;//kV/cm//for SCE OFF constant Efield
  double Cx=X_correction_hist_sce_off->GetBinContent(X_correction_hist_sce_off->FindBin(x));
  double Cyz=0.0;
  if(x<0){
    Cyz=YZ_correction_neg_hist_2_sce_off->GetBinContent(YZ_correction_neg_hist_2_sce_off->FindBin(z,y));
  }
  if(x>=0){
    Cyz=YZ_correction_pos_hist_2_sce_off->GetBinContent(YZ_correction_pos_hist_2_sce_off->FindBin(z,y));
  }
 double corrected_dqdx=dqdx*Cx*Cyz*normalisation_factor/calib_factor;
 return corrected_dqdx;
}






