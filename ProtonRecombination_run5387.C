#define ProtonRecombination_run5387_cxx
#include "ProtonRecombination_run5387.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <algorithm>
#include <math.h>
//#include "dedx_function_p4a_argoneut_box.h"
//#include "dedx_function_p4a_box_new.h"
//#include "dedx_function_p4a_box_iter_x.h"
#include "dedx_function_p4a_birks_iter_x.h"
//#include "dedx_function_p4a_box_iter_x_endpoint_argoneut.h"
//#include "dedx_function_p4a_birks.h"
#include "BraggFitter.h"
#include "SungBinFitter.h"


bool ProtonRecombination_run5387::IsNotEmptyCheck(Long64_t entry){
  bool IsbeamPosEmpty=beamPosz->empty();
  bool IsBeamEmpty=beamMomentum->empty();
  bool IsPrimtrk_hitEmpty=primtrk_hitz->empty();
  bool IsPrimtrk_start=primtrk_startz->empty();
  bool IsPrimtrk_dqdxEmpty=primtrk_dqdx->empty();
  bool IsPrimtrklen=primtrklen->empty();
  bool IsPrimtrk_range=primtrk_range->empty();
  bool IsPrimtrk_resrange=primtrk_resrange->empty();
  bool IsPrimtrk_pitch=primtrk_pitch->empty();
  bool IsEmpty=false;
  if (IsPrimtrk_start||IsPrimtrk_hitEmpty||IsBeamEmpty||IsbeamPosEmpty||IsPrimtrk_dqdxEmpty||IsPrimtrklen||IsPrimtrk_range||IsPrimtrk_resrange||IsPrimtrk_pitch) { //if any of these containers are empty                                                           
      IsEmpty=true;
  }
  return !IsEmpty;
}

bool ProtonRecombination_run5387::IsStoppingProton(Long64_t entry, TGraph* &graph){

  bool IsStoppingP = false;
  double primtrk_range_reco = primtrk_range->at(0);
  double csda_val_spec = graph->Eval(beamMomentum->at(0));
  //double mean_norm_trklen_csda=9.06809e-01; //prod3 spec           
  //double sigma_norm_trklen_csda=6.89073e-02; //prod3 spec          
  // WHY DIFFERENT?!?!?!------------
  //double mean_norm_trklen_csda=8.82596e-01; //new2   
  //double sigma_norm_trklen_csda=7.06573e-02; //new2 
  // -------------------------------
  double mean_norm_trklen_csda=8.84707e-01; //prod4               
  double sigma_norm_trklen_csda=7.02962e-02; //prod4      
  double min_norm_trklen_csda=mean_norm_trklen_csda-2.*sigma_norm_trklen_csda;
  double max_norm_trklen_csda=mean_norm_trklen_csda+3.*sigma_norm_trklen_csda;
  if ((primtrk_range_reco/csda_val_spec)>min_norm_trklen_csda&&(primtrk_range_reco/csda_val_spec)<max_norm_trklen_csda) { //csda from spec
    IsStoppingP=true;
  }
  return IsStoppingP;
}


bool ProtonRecombination_run5387::PassesBeamCuts(Long64_t entry){

  // the beam end position cut in mc is not needed

  // reco start position of track consistent with beam?            
  bool IsPos=false;
  //new p3                                    
  // double min1_x=-3.05368e+01-3.*5.13604e+00;
  //double min2_x=-3.05368e+01+3.*5.13604e+00;
  // double min1_y=4.22298e+02-3.*4.22119e+00;
  //double min2_y=4.22298e+02+3.*4.22119e+00;
  //double min1_z=5.72307e-02-3.*1.86602e-01;
  //double min2_z=5.72307e-02+3.*1.86602e-01;
  // different, why ?!?!?!?!
  /*
  double min1_z=3.67892e+00-3.*1.02152e+00; //new2 
  double min2_z=3.67892e+00+3.*1.02152e+00; //new2
  double min1_y=4.24055e+02-3.*4.54622e+00; //new2
  double min2_y=4.24055e+02+3.*4.54622e+00; //new2
  double min1_x=-2.82470e+01-3.*3.83924e+00; //new2
  double min2_x=-2.82470e+01+3.*3.83924e+00; //new2  
  */

  double min1_x=-2.82351e+01-3.*3.98756e+00; //prod4              
  double min2_x=-2.82351e+01+3.*3.98756e+00; //prod4              
  double min1_y=4.24140e+02-3.*4.68599e+00; //prod4               
  double min2_y=4.24140e+02+3.*4.68599e+00; //prod4               
  double min1_z=3.79565e+00-3.*1.01964e+00; //prod4               
  double min2_z=3.79565e+00+3.*1.01964e+00; //prod4  

  double primtrk_x_st_reco=primtrk_hitx->at(0);
  double primtrk_y_st_reco=primtrk_hity->at(0);
  double primtrk_z_st_reco=primtrk_hitz->at(0);
  
  // apparently the tracks can be flipped in data but not mc?
  if (primtrk_hitz->at(0)>primtrk_hitz->at(primtrk_hitz->size()-1)){
      primtrk_x_st_reco=primtrk_hitx->at(primtrk_hitx->size()-1);
      primtrk_y_st_reco=primtrk_hity->at(primtrk_hity->size()-1);
      primtrk_z_st_reco=primtrk_hitz->at(primtrk_hitz->size()-1);
  }


  if (primtrk_z_st_reco>min1_z&&primtrk_z_st_reco<min2_z) {
    if (primtrk_y_st_reco>min1_y&&primtrk_y_st_reco<min2_y) {
      if (primtrk_x_st_reco>min1_x&&primtrk_x_st_reco<min2_x) {
        IsPos=true;
      }
    }
  }

  // reco angle of track consistent with beam?                     
  bool IsCosine=false;
  //double cosine_beam_primtrk_min=(9.93062e-01)-4.*(3.53604e-03); //new p3 [spec]          
  // different, why?!?!?!
  //double cosine_beam_primtrk_min=9.78573e-01-4.*1.10895e-02;// new new p3
  double cosine_beam_primtrk_min=9.78550e-01-4.*1.10031e-02; // p4

  //cosine between beam_spec and primary trk direction   
  //double cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2];
  double cosine_beam_spec_primtrk = cosine_beam_primtrk;
  if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
  if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }

  //std::cout << IsCosine << " " << IsPos << std::endl;
  return (IsCosine and IsPos);
}




double ProtonRecombination_run5387::Calibrate_dQdx(Long64_t entry, int j, float prim_dqdx, float normalisation_factor, float calib_factor, TH1F* &X_correction_hist_2, TH2F* &YZ_correction_neg_hist_2){
  double prim_hitx=primtrk_hitx->at(j);
  double prim_hity=primtrk_hity->at(j);
  double prim_hitz=primtrk_hitz->at(j);
  double prim_pitch=primtrk_pitch->at(j);

  int x_bin = X_correction_hist_2->FindBin(prim_hitx);
  float Cx = X_correction_hist_2->GetBinContent(x_bin);
  double corrected_dq_dx = 0.0;

  if (prim_hitx>-360&&prim_hitx<0){ //negative X direction  
    float Cyzneg=YZ_correction_neg_hist_2->GetBinContent(YZ_correction_neg_hist_2->FindBin(prim_hitz,prim_hity));
    corrected_dq_dx=prim_dqdx*Cx*normalisation_factor*Cyzneg;
  }

  if(prim_hitx>0&&prim_hitx<360){ //positive X direction 
    float Cyzpos=YZ_correction_pos_hist_2->GetBinContent(YZ_correction_pos_hist_2->FindBin(prim_hitz,prim_hity));
    corrected_dq_dx=prim_dqdx*Cx*normalisation_factor*Cyzpos;
  }

  double scaled_corrected_dq_dx=float(corrected_dq_dx)/calib_factor;
  return scaled_corrected_dq_dx;
}



void ProtonRecombination_run5387::Loop()
{

  TFile *fmom_csda=new TFile("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/proton_mom_csda_converter.root");
  TGraph *csda_range_vs_mom_sm=(TGraph *)fmom_csda->Get("csda_range_vs_mom_sm");

  // to limit the number of entries                                         
  int max_entries = 100;//1000000;                                      
  bool limit_entries = false;

  // output histograms                                                          
  TH2D* dQdx_cal = new TH2D("dQdx_cal","dQdx_cal",240,0.,120.,350,0.,350000.); // vs remaining track length                                                    
  TH2D* dEdx_hyp = new TH2D("dEdx_hyp","dEdx_hyp",240,0.,120.,300,0.,30.); // vs remaining track length                                                        
  TH2D* dEdx_vs_dQdx = new TH2D("dEdx_vs_dQdx","dEdx_vs_dQdx",300,0.,30.,350,0.,350000.);

  // initialize Sungbin's fitter
  sf_initialize();


  // output tree                                                         
  //TFile* fout = new TFile("output/data_proton_recombination_birks.root", "RECREATE");
  TFile* fout = new TFile("output/data_proton_recombination.root", "RECREATE");
  TTree* output_tree = new TTree("output_tree","output_tree");
  const int kMaxHit = 1000;
  int event_number;
  int nhits;
  double rr_shift;
  double rr_shift_sungbin;
  double chi2_sungbin;
  int rr_fit_is_good;
  double rr_chi2;
  double prim_dqdx[kMaxHit];
  double prim_resrange[kMaxHit];
  double prim_resrange_reverse[kMaxHit];
  double cali_dqdx[kMaxHit];
  double cali_dqdx_with_dedx[kMaxHit];
  double dedx_hyp[kMaxHit];
  double dedx[kMaxHit];
  double dx[kMaxHit];
  double pitch[kMaxHit];
  output_tree->Branch("event",&event_number,"event_number/I");
  output_tree->Branch("nhits",&nhits,"nhits/I");
  output_tree->Branch("rr_shift",&rr_shift,"rr_shift/D");
  output_tree->Branch("rr_shift_sungbin",&rr_shift,"rr_shift_sungbin/D");
  output_tree->Branch("chi2_sungbin",&chi2_sungbin,"chi2_sungbin/D");
  output_tree->Branch("rr_fit_is_good",&rr_fit_is_good,"rr_fit_is_good/I");
  output_tree->Branch("rr_chi2",&rr_chi2,"rr_chi2/D");
  output_tree->Branch("dqdx",&prim_dqdx,"dqdx[nhits]/D");
  output_tree->Branch("resrange",&prim_resrange,"resrange[nhits]/D");
  output_tree->Branch("cali_dqdx",&cali_dqdx,"cali_dqdx[nhits]/D");
  output_tree->Branch("cali_dqdx_with_dedx",&cali_dqdx_with_dedx,"cali_dqdx_with_dedx[nhits]/D");
  output_tree->Branch("dedx_hyp",&dedx_hyp,"dedx_hyp[nhits]/D");
  output_tree->Branch("dedx",&dedx,"dedx[nhits]/D");
  output_tree->Branch("dx",&dx,"dx[nhits]/D");
  output_tree->Branch("pitch",&pitch,"pitch[nhits]/D");



  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;


  for (Long64_t jentry=0; jentry<nentries;jentry++) { //evt loop
    Long64_t ientry = LoadTree(jentry);
    event_number = jentry;
    if (ientry < 0) break;
    if ((ientry > max_entries) and limit_entries) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
		

    // make sure it is not empty
    if( !IsNotEmptyCheck(jentry) ) continue;
    
    // select protons
    // not needed for data, already all selected by beamline
    
		
    // select stopping protons
    if (!IsStoppingProton(ientry, csda_range_vs_mom_sm)) continue;

    // apply beam cuts on start position and angle
    if (!PassesBeamCuts(ientry)) continue;

    // these got past the cuts
    std::cout << ientry << std::endl;


    // loop over reco hits of the primary track                                                               
    int N = primtrk_dedx->size();
    nhits = N;   

    // loop over reco hits of the primary track                  
    for (size_t h=0; h<primtrk_dedx->size(); ++h) {
      

      // get the raw dQ/dx values                                                                             
      prim_dqdx[h] = primtrk_dqdx->at(h);
      prim_resrange[h] = primtrk_resrange->at(h);
      prim_resrange_reverse[N-h] = prim_resrange[h];

      // calibrate the dQ/dx values                                                                           
      cali_dqdx[h] = dqdx_uniform_corr_no_dedx(prim_dqdx[h],primtrk_hitx->at(h),primtrk_hity->at(h),primtrk_hitz->at(h));
      //cali_dqdx[h] = dqdx_uniform_corr(prim_dqdx[h],primtrk_hitx->at(h),primtrk_hity->at(h),primtrk_hitz->at(h)); // CAREFULL ABBEY
      cali_dqdx_with_dedx[h] = dqdx_uniform_corr(prim_dqdx[h],primtrk_hitx->at(h),primtrk_hity->at(h),primtrk_hitz->at(h));

      // "default" dedx (don't use for recombination study!)                          
      dedx[h] = dedx_function_p4(prim_dqdx[h],primtrk_hitx->at(h),primtrk_hity->at(h),primtrk_hitz->at(h));
    

      // get the dx values                                             
      if( h!=0 ){
	double dX = primtrk_hitx->at(h)-primtrk_hitx->at(h-1);
	double dY = primtrk_hity->at(h)-primtrk_hity->at(h-1);
	double dZ = primtrk_hitz->at(h)-primtrk_hitz->at(h-1);
	dx[h] = sqrt(dX*dX + dY*dY + dZ*dZ);
      }else{
	double dX = primtrk_hitx->at(h+1)-primtrk_hitx->at(h);
	double dY = primtrk_hity->at(h+1)-primtrk_hity->at(h);
	double dZ = primtrk_hitz->at(h+1)-primtrk_hitz->at(h);
	dx[h] = sqrt(dX*dX + dY*dY + dZ*dZ);
      }

      // get the pitch values                                          
      pitch[h] = primtrk_pitch->at(h);
    }

    /*
    // is the track flipped?     THIS DOESN'T SEEM TO BE NEEDED 
    if (primtrk_hitz->at(0)>primtrk_hitz->at(primtrk_hitz->size()-1)) {
      std::cout << "reversing flipped track!" << std::endl;
      for( size_t h=0; h<N; ++h){
	prim_resrange[h] = prim_resrange_reverse[h];
      }
    }
    */
    

    // do Pip's Bragg fitting                                            
    if (nhits < 100) std::cout << "WARNING: fewer than 100 hits!" << std::endl;
    double dedx_x[100];
    double dedx_y[100];

    int counter = 0;
    for (int j = nhits-1; j >= 0; j--)
      {
	if (counter > 99) break;
	dedx_x[counter] = prim_resrange[j];
	dedx_y[counter] = cali_dqdx[j];
	counter++;
      }

    // find minimum resrange 
    int min_rr_index = 99;
    for( int j=0; j < nhits; j++ ){
      if( prim_resrange[j] < prim_resrange[min_rr_index] ) min_rr_index = j;
    }

    // find maximum dedx         
    int max_dedx_index = 99;
    for( int j=0; j < nhits; j++ ){
      if (cali_dqdx[j] > cali_dqdx[max_dedx_index]) max_dedx_index = j;
    }

    // low end charge checking            
    bool low_end_charge = false;
    if (min_rr_index != max_dedx_index) low_end_charge = true;


    TGraph *dEdx_gr = new TGraph(100, dedx_x, dedx_y);
    FitResult result;
    double min_rr = dedx_x[max_dedx_index-1];
    result = ExtractTrackShift(dEdx_gr, min_rr);
    rr_shift = -1.0*result.trackshift;
    rr_fit_is_good = 0;
    if (!low_end_charge) rr_fit_is_good = 1;
    rr_chi2 = result.chi2;





    // Sungbin's endpoint finding                          
    const std::vector<double> vector_dedx(dedx, dedx + sizeof dedx / sizeof dedx[0]);
    // also flip for data at truncation stage
    std::vector<double> vector_dedx_trunc_flip = std::vector<double>(vector_dedx.begin(), vector_dedx.begin()+nhits);
    std::reverse(vector_dedx_trunc_flip.begin(), vector_dedx_trunc_flip.end());


    const std::vector<double> vector_prim_resrange(prim_resrange, prim_resrange + sizeof prim_resrange / sizeof prim_resrange[0]);
    // also flip for data at truncation stage
    std::vector<double> vector_prim_resrange_trunc_flip = std::vector<double>(vector_prim_resrange.begin(), vector_prim_resrange.begin()+nhits);
    std::reverse(vector_prim_resrange_trunc_flip.begin(), vector_prim_resrange_trunc_flip.end());


    std::pair<double, double> thingy = Fit_dEdx_Residual_Length(vector_dedx_trunc_flip, vector_prim_resrange_trunc_flip);
    std::cout << "Sungbin length: " << thingy.first << std::endl;
    double original_length = vector_prim_resrange_trunc_flip[0];
    std::cout << "Original length: " << original_length << std::endl;
    double sungbin_offset = thingy.first - original_length;
    rr_shift_sungbin = -sungbin_offset;
    chi2_sungbin = thingy.second;








    // fill the histograms                                                                                    
    for (size_t h=0; h<N; ++h){
      dQdx_cal->Fill(prim_resrange[h], cali_dqdx[h]);
      // get the (dE/dx)_{hyp} values                                                                         
      dedx_hyp[h] = 17.0*pow(prim_resrange[h],-0.42); // see ArgoNeuT paper                               
      dEdx_hyp->Fill(prim_resrange[h], dedx_hyp[h]);
      dEdx_vs_dQdx->Fill(dedx_hyp[h], cali_dqdx[h]);
    }
    
    output_tree->Fill();

  }
  
  // output the histograms                                        
  dQdx_cal->Write();
  dEdx_hyp->Write();
  dEdx_vs_dQdx->Write();

  // output the tree
  output_tree->Write();

  fout->Close();
}
