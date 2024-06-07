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

// MC Version!
double calib_factor = 1.0205e-3; //prod 4a
double normalisation_factor=1.0;//for plane 2 // prod 4 a

///change these parameters everytime you use new sample
//double calib_factor =5.02e-3;
//double normalisation_factor=1.02;//for plane 2
//TFile my_file0("./cali_oldp3/YZcalo_r6256030.root");//YZ correction factors
//TFile my_file2("./cali_oldp3/Xcalo_r6256030.root"); //X correction factors
TFile my_file0("./cali_p4/YZcalo_prod4_sceon.root");//YZ correction factors   
TFile my_file2("./cali_p4/Xcalo_prod4_sceon.root"); //X correction factors     
bool sceoff=0;//change this to 0 if you want sceon
////change parameters ends****************************************//

double Rho = 1.383;//g/cm^3 (liquid argon density at a pressure 18.0 psia) 
double betap = 0.212;//(kV/cm)(g/cm^2)/MeV
double alpha = 0.93;//parameter from ArgoNeuT experiment at 0.481kV/cm 
double Wion = 23.6e-6;//parameter from ArgoNeuT experiment at 0.481kV/cm
double Efield = 0.4867;//kV/cm protoDUNE electric filed

const int Z=18; //Atomic number of Argon
const double A=39.948; // g/mol Atomic mass of Argon
const double I=188.0e-6; // ev
const double K=0.307; // Mev.cm^2 / mol
const double Mmu=105.658; // Mev for Mu
const double Me=0.51; // Mev for electron
const double rho=1.396;//g/cm^3


////getting the variable Efield using data driven maps
TFile *ef=new TFile("./cali_p4/SCE_DataDriven_180kV_v4.root");
TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");
double tot_Ef(double xval,double yval,double zval){
  double E0v=0.4867;
  if(xval>=0){
    double ex=E0v+E0v*xpos->GetBinContent(xpos->FindBin(xval,yval,zval));
    double ey=E0v*ypos->GetBinContent(ypos->FindBin(xval,yval,zval));
    double ez=E0v*zpos->GetBinContent(zpos->FindBin(xval,yval,zval));
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
if(xval<0){
    double ex=E0v+E0v*xneg->GetBinContent(xneg->FindBin(xval,yval,zval));
    double ey=E0v*yneg->GetBinContent(yneg->FindBin(xval,yval,zval));
    double ez=E0v*zneg->GetBinContent(zneg->FindBin(xval,yval,zval));
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
 else
   return E0v;
}

/////////////////////////
//////////////////////Importing X fractional corrections//////////////////////
TH1F *X_correction_hist = (TH1F*)my_file2.Get("dqdx_X_correction_hist_2");
TH2F *YZ_correction_neg_hist_2=(TH2F*)my_file0.Get("correction_dqdx_ZvsY_negativeX_hist_2");
TH2F *YZ_correction_pos_hist_2=(TH2F*)my_file0.Get("correction_dqdx_ZvsY_positiveX_hist_2");
 ////////////////////////////////////////////////////////////////////////////////// 

Double_t dedx_function_35ms(double dqdx, double x, double y, double z){
  double Ef;
  Ef=tot_Ef(x,y,z);//for SCE ON
  if(sceoff) Ef=0.4867;//kV/cm//for SCE OFF constant Efield
  double Cx=X_correction_hist->GetBinContent(X_correction_hist->FindBin(x));
  double Cyz=0.0;
  if(x<0){
    Cyz=YZ_correction_neg_hist_2->GetBinContent(YZ_correction_neg_hist_2->FindBin(z,y));
  }
  if(x>=0){
    Cyz=YZ_correction_pos_hist_2->GetBinContent(YZ_correction_pos_hist_2->FindBin(z,y));
  }
 double corrected_dqdx=dqdx*Cx*Cyz*normalisation_factor/calib_factor;
 return (exp(corrected_dqdx*(betap/(Rho*Ef)*Wion))-alpha)/(betap/(Rho*Ef));
}


//HY: added the following func for the recombination study
Double_t dqdx_uniform_corr(double dqdx, double x, double y, double z, bool sc_off=false){
  double Ef;
  Ef=tot_Ef(x,y,z);//for SCE ON
  if(sc_off) Ef=0.4867;//kV/cm//for SCE OFF constant Efield
  double Cx=X_correction_hist->GetBinContent(X_correction_hist->FindBin(x));
  double Cyz=0.0;
  if(x<0){
    Cyz=YZ_correction_neg_hist_2->GetBinContent(YZ_correction_neg_hist_2->FindBin(z,y));
  }
  if(x>=0){
    Cyz=YZ_correction_pos_hist_2->GetBinContent(YZ_correction_pos_hist_2->FindBin(z,y));
  }
 double corrected_dqdx=dqdx*Cx*Cyz*normalisation_factor/calib_factor;
 return corrected_dqdx;
}

Double_t dqdx_uniform_corr_no_dedx(double dqdx, double x, double y, double z, bool sc_off=false){
  double Ef;
  Ef=tot_Ef(x,y,z);//for SCE ON
  if(sc_off) Ef=0.4867;//kV/cm//for SCE OFF constant Efield
  double Cx=X_correction_hist->GetBinContent(X_correction_hist->FindBin(x));
  double Cyz=0.0;
  if(x<0){
    Cyz=YZ_correction_neg_hist_2->GetBinContent(YZ_correction_neg_hist_2->FindBin(z,y));
  }
  if(x>=0){
    Cyz=YZ_correction_pos_hist_2->GetBinContent(YZ_correction_pos_hist_2->FindBin(z,y));
  }
  double corrected_dqdx=dqdx*Cx*Cyz*normalisation_factor/(1.0e-3);
 return corrected_dqdx;
}

////////////////////////////////////////////////////////////////////////////////
