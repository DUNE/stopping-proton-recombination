#define ProtonRecombination_cxx
#include "ProtonRecombination.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
//#include "dedx_function_35ms_p4a.h"
#include "dedx_function_my_mc.h"
//#include "dedx_function_p4a.h" // SYSTEMATICS FOR SCE - use wrong map
#include "dedx_function_35ms_SCE_OFF.h"
#include "BraggFitter.h"
#include "SungBinFitter.h"


bool ProtonRecombination::IsNotEmptyCheck(Long64_t entry){

  bool is_data_notempty=false; //if the data containers not empty
  if (primtrk_hitz->empty()==false){
    if (primtrk_hitz->size()&&beamtrk_z->size()){
      //if (primtrk_hitz->at(0).size()>1){
	is_data_notempty=true;
	// }
    }
  }
  return is_data_notempty;
}


bool ProtonRecombination::IsStoppingProton(Long64_t entry, TGraph* &graph){

  bool IsStoppingP = false;
  double primtrk_range_reco = primtrk_range->at(0);
  double csda_val_spec = graph->Eval(beamMomentum_spec->at(0)); 
  //std::cout << primtrk_range_reco/csda_val_spec << std::endl;
  //double mean_norm_trklen_csda=9.06809e-01; //prod3 spec 
  //double sigma_norm_trklen_csda=6.89073e-02; //prod3 spec 
  double mean_norm_trklen_csda=9.32064e-01; //prod4 spec          
  double sigma_norm_trklen_csda=1.32800e-02; //prod4 spec   
  double min_norm_trklen_csda=mean_norm_trklen_csda-2.*sigma_norm_trklen_csda;
  double max_norm_trklen_csda=mean_norm_trklen_csda+3.*sigma_norm_trklen_csda;
  if ((primtrk_range_reco/csda_val_spec)>min_norm_trklen_csda&&(primtrk_range_reco/csda_val_spec)<max_norm_trklen_csda) { //csda from spec   
    IsStoppingP=true;
  }
  return IsStoppingP;
}


bool ProtonRecombination::PassesBeamCuts(Long64_t entry){

  // does beam instrumentation suggest it reaches TPC?
  bool is_beam_at_ff=false; //if the beam reach tpc                                                    
  try{
            
  int key_reach_tpc=-99;
  if (beamtrk_z){
    for (size_t kk=0; kk<beamtrk_z->size(); ++kk) {  //loop over all beam hits                               
      double zpos_beam=beamtrk_z->at(kk);
      //      if ((zpos_beam+0.49375)<0.01&&(zpos_beam+0.49375)>-0.01) { key_reach_tpc=(int)kk; } //kind key at ff
      if (zpos_beam>=0) {
	key_reach_tpc=(int)kk;
	break;
      } //kind key at ff   

  } //loop over all beam hits                                                                              
  }
  if (key_reach_tpc!=-99) { is_beam_at_ff=true; }


  // reco start position of track consistent with beam?
  bool IsPos=false;
  //new p3
  /*                                                                     
  double min1_x=-3.05368e+01-3.*5.13604e+00;
  double min2_x=-3.05368e+01+3.*5.13604e+00;
  double min1_y=4.22298e+02-3.*4.22119e+00;
  double min2_y=4.22298e+02+3.*4.22119e+00;
  double min1_z=5.72307e-02-3.*1.86602e-01;
  double min2_z=5.72307e-02+3.*1.86602e-01;
  */

  
  //new p4                                                         
  double min1_z=5.10816e-02-3.*2.13366e-01;
  double min2_z=5.10816e-02+3.*2.13366e-01;
  double min1_y=4.21863e+02-3.*4.11359e+00;
  double min2_y=4.21863e+02+3.*4.11359e+00;
  double min1_x=-3.05895e+01-3.*4.69242e+00;
  double min2_x=-3.05895e+01+3.*4.69242e+00;
  

  // hack p4 for delta study
  /*
  double min1_z=5.10816e-02-9.*2.13366e-01;
  double min2_z=5.10816e-02+9.*2.13366e-01;
  double min1_y=4.21863e+02-3.*4.11359e+00;
  double min2_y=4.21863e+02+3.*4.11359e+00;
  double min1_x=-3.05895e+01-3.*4.69242e+00;
  double min2_x=-3.05895e+01+3.*4.69242e+00;
  */

  double primtrk_x_st_reco=primtrk_hitx->at(0);//[0];
  double primtrk_y_st_reco=primtrk_hity->at(0);//[0];
  double primtrk_z_st_reco=primtrk_hitz->at(0);//[0];
  

  /* std::cout << "*****" << std::endl;
  std::cout << min1_x << " " << primtrk_x_st_reco << " " << min2_x << std::endl;
  std::cout << min1_y << " " << primtrk_y_st_reco << " " << min2_y << std::endl;
  std::cout << min1_z << " " << primtrk_z_st_reco << " " << min2_z << std::endl;
  std::cout << "*****" << std::endl;
  */

  // check if tracks are flipped
  if (primtrk_hitz->at(0)>primtrk_hitz->at(primtrk_hitz->size()-1)) {
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
  // double cosine_beam_primtrk_min=(9.93062e-01)-4.*(3.53604e-03); //new p3 [spec]
  double cosine_beam_primtrk_min=(9.92194e-01)-4.*(3.96921e-03);//new p4 [spec]    
  //cosine between beam_spec and primary trk direction                                                             
  double cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2];
  if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
  if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }

  //std::cout << IsCosine << " " << IsPos << " " << is_beam_at_ff << std::endl;

  return (IsCosine and IsPos and is_beam_at_ff);}
  catch (...){
    std::cout << "error encountered" << std::endl;
    return false;
  }


}


void ProtonRecombination::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ProtonRecombination.C
//      root> ProtonRecombination t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

  
// read in the proton reference dEdx curves? // IS THIS RIGHT??? SHOULD BE MC VERSION?
  // TFile *fmom_csda=new TFile("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/proton_mom_csda_converter.root");
  TFile *fmom_csda=new TFile("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/nosce/proton_mom_csda_converter.root");

  // sce off version (be careful)
  //TFile *fmom_csda=new TFile("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/nosce/proton_mom_csda_converter.root");
  TGraph *csda_range_vs_mom_sm=(TGraph *)fmom_csda->Get("csda_range_vs_mom_sm");
  //TGraph *mom_vs_csda_range_sm=(TGraph *)fmom_csda->Get("mom_vs_csda_range_sm");


  // to limit the number of entries
  int max_entries = 10000;//1000000;
  bool limit_entries = false;

  // output histograms
  TH2D* dQdx_cal = new TH2D("dQdx_cal","dQdx_cal",240,0.,120.,350,0.,350000.); // vs remaining track length
  TH2D* dEdx_hyp = new TH2D("dEdx_hyp","dEdx_hyp",240,0.,120.,300,0.,30.); // vs remaining track length
  TH2D* dEdx_vs_dQdx = new TH2D("dEdx_vs_dQdx","dEdx_vs_dQdx",300,0.,30.,350,0.,350000.);


  // output tree
  TFile* fout = new TFile("output/mc_proton_recombination.root", "RECREATE");
  // alternative for SCE off
  //TFile* fout = new TFile("output/mc_proton_recombination_sce_off.root", "RECREATE");
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
  double cali_dqdx_sc_off[kMaxHit];
  double dedx_hyp[kMaxHit];
  double dedx[kMaxHit];
  double dx[kMaxHit];
  double pitch[kMaxHit];
  double x[kMaxHit];
  double y[kMaxHit];
  double z[kMaxHit];
  output_tree->Branch("event",&event_number,"event_number/I");
  output_tree->Branch("nhits",&nhits,"nhits/I");
  output_tree->Branch("rr_shift",&rr_shift,"rr_shift/D");
  output_tree->Branch("rr_shift_sungbin",&rr_shift_sungbin,"rr_shift_sungbin/D");
  output_tree->Branch("chi2_sungbin",&chi2_sungbin,"chi2_sungbin/D");
  output_tree->Branch("rr_fit_is_good",&rr_fit_is_good,"rr_fit_is_good/I");  
  output_tree->Branch("rr_chi2",&rr_chi2,"rr_chi2/D");
  output_tree->Branch("dqdx",&prim_dqdx,"dqdx[nhits]/D");
  output_tree->Branch("resrange",&prim_resrange,"resrange[nhits]/D");
  output_tree->Branch("cali_dqdx",&cali_dqdx,"cali_dqdx[nhits]/D");
  output_tree->Branch("cali_dqdx_sc_off",&cali_dqdx_sc_off,"cali_dqdx_sc_off[nhits]/D");
  output_tree->Branch("dedx_hyp",&dedx_hyp,"dedx_hyp[nhits]/D");
  output_tree->Branch("dedx",&dedx,"dedx[nhits]/D");
  output_tree->Branch("dx",&dx,"dx[nhits]/D");
  output_tree->Branch("pitch",&pitch,"pitch[nhits]/D");
  output_tree->Branch("x",&x,"x[nhits]/D");
  output_tree->Branch("y",&y,"y[nhits]/D");
  output_tree->Branch("z",&z,"z[nhits]/D");


  // initialise Sungbin's fitter
  sf_initialize();

  // fail counters
  int quality_fails = 0;
  int proton_fails = 0;
  int beam_pos_fails = 0;
  int stopping_proton_fails = 0;

  // success counter
  int successes = 0;


   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      event_number = jentry;
      if (ientry < 0) break;
      if ((ientry > max_entries) and limit_entries) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      // make sure the entry is not empty
      if (!IsNotEmptyCheck(ientry)) {
	quality_fails += 1;  continue; }      

      // select protons with simulated beam instrumentation
      if (beamtrackPdg != 2212) {
	proton_fails += 1; continue; }

      // select stopping protons using reconstruction
      if (!IsStoppingProton(ientry, csda_range_vs_mom_sm)) {
      	stopping_proton_fails += 1; continue;} 



      try {
      // apply beam cuts on start position and angle
      if (!PassesBeamCuts(ientry)) {
	beam_pos_fails += 1; continue;}
      } catch (...) {
	std::cout << "caught an error" << std::endl;
	continue;
      }
      // these got past the cuts
      std::cout << jentry << " " << truthpdg << " " << beamtrackPdg << std::endl;
      successes += 1;


      // loop over reco hits of the primary track 
      int N = primtrk_dedx->size();
      nhits = N;


      if (primtrk_dqdx == NULL) continue; // Pip
      //std::ostringstream oss;
      //oss << primtrk_dqdx;
      //if (oss.str() > "0xffff0000") continue; // well, this is ugly.
      // I think there must have been some corruption in there - Abbey



      for (size_t h=0; h<N; ++h) {
	prim_dqdx[h] = primtrk_dqdx->at(h);
	prim_resrange[h] = primtrk_resrange->at(h);
	prim_resrange_reverse[N-h] = prim_resrange[h];


	// calibrate the dQ/dx values
	cali_dqdx[h] = dqdx_uniform_corr_no_dedx(prim_dqdx[h],primtrk_hitx->at(h),primtrk_hity->at(h),primtrk_hitz->at(h));
	//cali_dqdx[h] = dqdx_uniform_corr(prim_dqdx[h],primtrk_hitx->at(h),primtrk_hity->at(h),primtrk_hitz->at(h));

	// version of calibration with space charge OFF
	cali_dqdx_sc_off[h] = dqdx_uniform_corr_sce_off(prim_dqdx[h],primtrk_hitx->at(h),primtrk_hity->at(h),primtrk_hitz->at(h));

	// "default" dedx (don't use for recombination study!)
	dedx[h] = dedx_function_35ms(prim_dqdx[h],primtrk_hitx->at(h),primtrk_hity->at(h),primtrk_hitz->at(h));
	//dedx[h] = dedx_function_p4(prim_dqdx[h],primtrk_hitx->at(h),primtrk_hity->at(h),primtrk_hitz->at(h)); // ABBEY HACK


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


	// position information for systematics
	x[h] = primtrk_hitx->at(h);
	y[h] = primtrk_hity->at(h);
	z[h] = primtrk_hitz->at(h);


	// get the pitch values
	pitch[h] = primtrk_pitch->at(h);

      }


      /*
      // is the track flipped? THIS DOESN'T SEEM TO BE NEEDED
      if (primtrk_hitz->at(0)[0]>primtrk_hitz->at(0)[primtrk_hitz->at(0).size()-1]) {
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
      if( !low_end_charge ) rr_fit_is_good = 1;
      rr_chi2 = result.chi2;



      // Sungbin's endpoint finding
      const std::vector<double> vector_cali_dqdx(cali_dqdx, cali_dqdx + sizeof cali_dqdx / sizeof cali_dqdx[0]);
      const std::vector<double> vector_cali_dqdx_trunc = std::vector<double>(vector_cali_dqdx.begin(), vector_cali_dqdx.begin()+nhits);      
      const std::vector<double> vector_dedx(dedx, dedx + sizeof dedx / sizeof dedx[0]);      
      const std::vector<double> vector_dedx_trunc = std::vector<double>(vector_dedx.begin(), vector_dedx.begin()+nhits);            
      const std::vector<double> vector_prim_resrange(prim_resrange, prim_resrange + sizeof prim_resrange / sizeof prim_resrange[0]);      
      const std::vector<double> vector_prim_resrange_trunc = std::vector<double>(vector_prim_resrange.begin(), vector_prim_resrange.begin()+nhits);            

      std::pair<double, double> thingy = Fit_dEdx_Residual_Length(vector_dedx_trunc, vector_prim_resrange_trunc);
      std::cout << "Sungbin length: " << thingy.first << std::endl;
      std::cout << "Original length: " << vector_prim_resrange[0] << std::endl;
      double sungbin_offset = thingy.first - vector_prim_resrange[0];
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

      // fill the tree
      output_tree->Fill();
      
   }

   // how many failed at each stage
   std::cout << proton_fails << " " << beam_pos_fails << " "
	     << quality_fails << " " << stopping_proton_fails << " "
	     << successes << std::endl;


   // output the histograms
   dQdx_cal->Write();
   dEdx_hyp->Write();
   dEdx_vs_dQdx->Write();

   // output the tree
   output_tree->Write();

   fout->Close();

}
