#define ExtractTracks_cxx
#include "ExtractTracks.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include "dedx_function_35ms.h"
#include "dedx_function_35ms_SCE_OFF.h"

bool ExtractTracks::IsNotEmptyCheck(Long64_t entry){

  bool is_data_notempty=false; //if the data containers not empty
  if (primtrk_hitz->empty()==false){
    if (primtrk_hitz->size()&&beamtrk_z->size()){
      if (primtrk_hitz->at(0).size()>1){
	is_data_notempty=true;
      }
    }
  }
  return is_data_notempty;
}


bool ExtractTracks::IsStoppingProton(Long64_t entry, TGraph* &graph){

  bool IsStoppingP = false;
  double primtrk_range_reco = primtrk_range->at(0);
  double csda_val_spec = graph->Eval(beamMomentum_spec->at(0)); 
  double mean_norm_trklen_csda=9.06809e-01; //prod3 spec 
  double sigma_norm_trklen_csda=6.89073e-02; //prod3 spec 
  double min_norm_trklen_csda=mean_norm_trklen_csda-2.*sigma_norm_trklen_csda;
  double max_norm_trklen_csda=mean_norm_trklen_csda+3.*sigma_norm_trklen_csda;
  if ((primtrk_range_reco/csda_val_spec)>min_norm_trklen_csda&&(primtrk_range_reco/csda_val_spec)<max_norm_trklen_csda) { //csda from spec   
    IsStoppingP=true;
  }
  return IsStoppingP;
}


bool ExtractTracks::PassesBeamCuts(Long64_t entry){

  // does beam instrumentation suggest it reaches TPC?
  bool is_beam_at_ff=false; //if the beam reach tpc                                                                
  int key_reach_tpc=-99;
  if (beamtrk_z){
    for (size_t kk=0; kk<beamtrk_z->size(); ++kk) {  //loop over all beam hits                               
      double zpos_beam=beamtrk_z->at(kk);
      if ((zpos_beam+0.49375)<0.01&&(zpos_beam+0.49375)>-0.01) { key_reach_tpc=(int)kk; } //kind key at ff
  } //loop over all beam hits                                                                              
  }
  if (key_reach_tpc!=-99) { is_beam_at_ff=true; }


  // reco start position of track consistent with beam?
  bool IsPos=false;
  //new p3                                                                                                         
  double min1_x=-3.05368e+01-3.*5.13604e+00;
  double min2_x=-3.05368e+01+3.*5.13604e+00;
  double min1_y=4.22298e+02-3.*4.22119e+00;
  double min2_y=4.22298e+02+3.*4.22119e+00;
  double min1_z=5.72307e-02-3.*1.86602e-01;
  double min2_z=5.72307e-02+3.*1.86602e-01;

  double primtrk_x_st_reco=primtrk_hitx->at(0)[0];
  double primtrk_y_st_reco=primtrk_hity->at(0)[0];
  double primtrk_z_st_reco=primtrk_hitz->at(0)[0];


  // check if tracks are flipped
  if (primtrk_hitz->at(0)[0]>primtrk_hitz->at(0)[primtrk_hitz->at(0).size()-1]) {
    primtrk_x_st_reco=primtrk_hitx->at(0)[primtrk_hitx->at(0).size()-1];
    primtrk_y_st_reco=primtrk_hity->at(0)[primtrk_hity->at(0).size()-1];
    primtrk_z_st_reco=primtrk_hitz->at(0)[primtrk_hitz->at(0).size()-1];
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
  double cosine_beam_primtrk_min=(9.93062e-01)-4.*(3.53604e-03); //new p3 [spec]  
  //cosine between beam_spec and primary trk direction                                                             
  double cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2];
  if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
  if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }

  return (IsCosine and IsPos and is_beam_at_ff);
}


void ExtractTracks::Loop()
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
  TFile *fmom_csda=new TFile("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/proton_mom_csda_converter.root");
  // sce off version (be careful)
  //TFile *fmom_csda=new TFile("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/nosce/proton_mom_csda_converter.root");
  TGraph *csda_range_vs_mom_sm=(TGraph *)fmom_csda->Get("csda_range_vs_mom_sm");
  //TGraph *mom_vs_csda_range_sm=(TGraph *)fmom_csda->Get("mom_vs_csda_range_sm");


  // to limit the number of entries
  int max_entries = 10000;//1000000;
  bool limit_entries = false;


  // output tree
  TFile* fout = new TFile("output/mc_extract_tracks.root", "RECREATE");
  // alternative for SCE off
  //TFile* fout = new TFile("output/mc_proton_recombination_sce_off.root", "RECREATE");
  TTree* output_tree = new TTree("output_tree","output_tree");
  const int kMaxHit = 1000;
  int event_number;
  int nhits;
  int true_pdg;
  int beamline_pdg;
  double prim_dqdx[kMaxHit];
  double prim_resrange[kMaxHit];
  double prim_resrange_reverse[kMaxHit];
  double cali_dqdx[kMaxHit];
  double cali_dqdx_sc_off[kMaxHit];
  double dedx_hyp[kMaxHit];
  double dedx[kMaxHit];
  double dx[kMaxHit];
  double hit_x[kMaxHit];
  double hit_y[kMaxHit];
  double hit_z[kMaxHit];
  double pitch[kMaxHit];
  output_tree->Branch("event",&event_number,"event_number/I");
  output_tree->Branch("true_pdg",&true_pdg,"true_pdg/I");
  output_tree->Branch("beamline_pdg",&beamline_pdg,"beamline_pdg/I");
  output_tree->Branch("nhits",&nhits,"nhits/I");
  output_tree->Branch("dqdx",&prim_dqdx,"dqdx[nhits]/D");
  output_tree->Branch("resrange",&prim_resrange,"resrange[nhits]/D");
  output_tree->Branch("cali_dqdx",&cali_dqdx,"cali_dqdx[nhits]/D");
  output_tree->Branch("cali_dqdx_sc_off",&cali_dqdx_sc_off,"cali_dqdx_sc_off[nhits]/D");
  output_tree->Branch("dedx_hyp",&dedx_hyp,"dedx_hyp[nhits]/D");
  output_tree->Branch("dedx",&dedx,"dedx[nhits]/D");
  output_tree->Branch("dx",&dx,"dx[nhits]/D");
  output_tree->Branch("hit_x",&hit_x,"hit_x[nhits]/D");
  output_tree->Branch("hit_y",&hit_y,"hit_y[nhits]/D");
  output_tree->Branch("hit_z",&hit_z,"hit_z[nhits]/D");
  output_tree->Branch("pitch",&pitch,"pitch[nhits]/D");


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
      if (!IsNotEmptyCheck(ientry)) continue;      

      // select protons with simulated beam instrumentation
      //if (beamtrackPdg != 2212) continue;

      // select stopping protons using reconstruction
      //if (!IsStoppingProton(ientry, csda_range_vs_mom_sm)) continue; 

      // apply beam cuts on start position and angle
      if (!PassesBeamCuts(ientry)) continue;

      // these got past the cuts
      true_pdg = truthpdg;
      beamline_pdg = beamtrackPdg;
      std::cout << jentry << " " << truthpdg << " " << beamtrackPdg << std::endl;


      // loop over reco hits of the primary track 
      int N = primtrk_dedx->at(0).size();
      nhits = N;

      for (size_t h=0; h<N; ++h) {
	// get the raw dQ/dx values
	prim_dqdx[h] = primtrk_dqdx->at(0)[h];
	prim_resrange[h] = primtrk_resrange->at(0)[h];
	prim_resrange_reverse[N-h] = prim_resrange[h];

	// calibrate the dQ/dx values
	cali_dqdx[h] = dqdx_uniform_corr(prim_dqdx[h],primtrk_hitx->at(0)[h],primtrk_hity->at(0)[h],primtrk_hitz->at(0)[h]);

	// version of calibration with space charge OFF
	cali_dqdx_sc_off[h] = dqdx_uniform_corr_sce_off(prim_dqdx[h],primtrk_hitx->at(0)[h],primtrk_hity->at(0)[h],primtrk_hitz->at(0)[h]);


	// "default" dedx (don't use for recombination study!)
	dedx[h] = dedx_function_35ms(prim_dqdx[h],primtrk_hitx->at(0)[h],primtrk_hity->at(0)[h],primtrk_hitz->at(0)[h]);

	// get the dx values
	if( h!=0 ){
          double dX = primtrk_hitx->at(0)[h]-primtrk_hitx->at(0)[h-1];
          double dY = primtrk_hity->at(0)[h]-primtrk_hity->at(0)[h-1];
          double dZ = primtrk_hitz->at(0)[h]-primtrk_hitz->at(0)[h-1];
          dx[h] = sqrt(dX*dX + dY*dY + dZ*dZ);
        }else{
          double dX = primtrk_hitx->at(0)[h+1]-primtrk_hitx->at(0)[h];
          double dY = primtrk_hity->at(0)[h+1]-primtrk_hity->at(0)[h];
          double dZ = primtrk_hitz->at(0)[h+1]-primtrk_hitz->at(0)[h];
          dx[h] = sqrt(dX*dX + dY*dY + dZ*dZ);
        }

	// the hit positions
	hit_x[h] = primtrk_hitx->at(0)[h];
	hit_y[h] = primtrk_hity->at(0)[h];
	hit_z[h] = primtrk_hitz->at(0)[h];


	// get the pitch values
	pitch[h] = primtrk_pitch->at(0)[h];
      }


      // fill the tree
      output_tree->Fill();
      

   }


   // output the tree
   output_tree->Write();

   fout->Close();

}
