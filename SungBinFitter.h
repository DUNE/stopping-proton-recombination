// Stolen from Sungbin!  Any mistakes my own.
#include "TMath.h"
#include "TVector3.h"
#include <iostream>
#include <map>
#include "BetheBloch.h"

std::map< int, BetheBloch* > map_BB;

void sf_initialize(){
  map_BB[2212] = new BetheBloch(2212);
  return;
}

std::pair<double, double> Fit_dEdx_Residual_Length(const vector<double> & dEdx, const vector<double> & ResRange, int PID=2212){

  std::pair<double, double> default_thing(-9999.,99999.);


  int N_max = 200; // == Maximum number of hits used for the Bethe-Bloch fitting


  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  int abs_PID = abs(PID);
  if(!(abs(PID) == 13 || PID == 2212 || abs(PID) == 211)){
    return default_thing;
  }


  // the fitter (is this slow?  maybe move)
 
  //map_BB[13] = new BetheBloch(13);
  //map_BB[211] = new BetheBloch(211);
  //map_BB[321] = new BetheBloch(321);




  double best_additional_res_length = -0.1;
  double best_chi2 = 99999.;
  double min_additional_res_length = -50.; // == [cm]
  double max_additional_res_length = 50.; // == [cm]
  //double min_additional_res_length = -10.; // == [cm]
  //double max_additional_res_length = 10.; // == [cm]

  double res_length_step = 0.01; //0.2; // == [cm]
  int N_skip = 3;
  int this_N_calo = dEdx.size();
  if(this_N_calo <= 15){
    std::cout << "sbf: too few hits" << std::endl;
    return default_thing;
    //return std::pair<double, double> thing(-9999., 99999.); // == Too small number of hits
  }

  int i_bestfit = -1;
  double dEdx_truncate_upper = 5.;
  double dEdx_truncate_bellow = 0.5;
  if(PID == 2212){
    max_additional_res_length = 50.;
    dEdx_truncate_upper = 20.;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;
  int this_N_hits = TMath::Min(this_N_calo, N_max); // == Use how many hits
  vector<double> chi2_vector;
  vector<double> additional_res_legnth_vector;

  double original_res_length = ResRange.at(this_N_calo - 1) - ResRange.at(this_N_calo - this_N_hits); // == [cm]
  original_res_length = ResRange.at(0) - ResRange.at(this_N_hits - 1); // == [cm]
  //std::cout << "INITIAL: " << original_res_length << std::endl; 

  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_chi2 = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip hits
      double this_res_length = ResRange.at(j) + this_additional_res_length;

      double this_KE = map_BB[abs_PID]->KEFromRangeSpline(this_res_length);
      //double dEdx_theory = dEdx_Bethe_Bloch(this_KE, this_mass);
      double dEdx_theory = map_BB[abs_PID]->meandEdx(this_KE);
      double dEdx_measured = dEdx.at(j);
      //std::cout << "looping " << i << " " << j << " " 
      //	<< dEdx_theory << " " << dEdx_measured << std::endl;
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate, it should be modified to consider protons

      // == Gaussian approx.
      //double dEdx_theory_err = dEdx_theory * 0.02;
      this_chi2 += pow(dEdx_measured - dEdx_theory, 2);
    }
    this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    //std::cout << "sbf: chi2/ndof = " << this_chi2 << "i = " << i << std::endl;
    if(this_chi2 < best_chi2){
      best_chi2 = this_chi2;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
  }




  double best_total_res_length = best_additional_res_length + original_res_length;
  
  double best_KE = map_BB[abs_PID]->KEFromRangeSpline(best_total_res_length);
  double best_mom = map_BB[abs_PID]->KEtoMomentum(best_KE);

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){
    cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : no mimumum" << endl;
    best_chi2 = 99999.;
    return default_thing;
    //return std::pair<double, double> thing(-9999., 999999.);
  }
  else if(best_chi2 > 99990.){
    cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 > 99990." << endl;
    best_chi2 = 99999.;
    return default_thing;
    //return std::pair<double, double> thing(-9999., 999999.);
  }
  else if(best_chi2 < 1.0e-11){
    cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 < 1.0e-11" << endl;
    best_chi2 = 99999.;
    return default_thing;
    //return std::pair<double, double> thing(-9999., 999999.);
  }
  //std::cout << "FINAL: " << best_total_res_length << std::endl;
  std::pair<double, double> thing(best_total_res_length, best_chi2);
  return thing; 
}
