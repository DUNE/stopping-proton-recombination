//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 19 18:06:27 2021 by ROOT version 6.18/04
// from TChain protonbeamana/beamana/
//////////////////////////////////////////////////////////

#ifndef ProtonRecombination_run5387_h
#define ProtonRecombination_run5387_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class ProtonRecombination_run5387 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Double_t        evttime;
   Int_t           Nactivefembs[5];
   Bool_t          beam_inst_valid;
   Int_t           beam_inst_trigger;
   Int_t           beam_inst_nTracks;
   Int_t           beam_inst_nMomenta;
   vector<int>     *beam_inst_PDG_candidates;
   vector<double>  *beam_inst_TOF;
   vector<int>     *beam_inst_TOF_Chan;
   Int_t           isprimarytrack;
   Int_t           isprimaryshower;
   vector<double>  *beamPosx;
   vector<double>  *beamPosy;
   vector<double>  *beamPosz;
   vector<double>  *beamDirx;
   vector<double>  *beamDiry;
   vector<double>  *beamDirz;
   vector<double>  *beamMomentum;
   Double_t        tof;
   vector<double>  *tofs;
   vector<int>     *ch_tofs;
   Short_t         low_pressure_status;
   Short_t         high_pressure_status;
   Double_t        low_pressure;
   Double_t        high_pressure;
   Double_t        cosine_beam_primtrk;
   vector<double>  *primtrk_startx;
   vector<double>  *primtrk_starty;
   vector<double>  *primtrk_startz;
   vector<double>  *primtrk_endx;
   vector<double>  *primtrk_endy;
   vector<double>  *primtrk_endz;
   vector<double>  *primtrk_Dirx;
   vector<double>  *primtrk_Diry;
   vector<double>  *primtrk_Dirz;
   vector<double>  *primtrklen;
   vector<double>  *primtrkID;
   vector<double>  *primtrk_dqdx;
   vector<double>  *primtrk_dedx;
   vector<double>  *primtrk_resrange;
   vector<double>  *primtrk_range;
   vector<double>  *primtrk_hitx;
   vector<double>  *primtrk_hity;
   vector<double>  *primtrk_hitz;
   vector<double>  *primtrk_ke;
   vector<double>  *primtrk_pitch;
   vector<int>     *primtrk_wid;
   vector<double>  *primtrk_pt;
   vector<unsigned long> *primtrk_calo_hit_index;
   vector<int>     *primtrk_ch;
   vector<int>     *pdg_code;
   vector<int>     *n_beamparticle;
   vector<int>     *n_daughter;
   vector<int>     *isPrimary;
   vector<int>     *pfp_self;
   vector<int>     *pfp_daughter;
   vector<int>     *primtrk_trktag;
   vector<double>  *Zintersection;
   vector<double>  *Yintersection;
   vector<double>  *Xintersection;
   vector<double>  *timeintersection;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_evttime;   //!
   TBranch        *b_Nactivefembs;   //!
   TBranch        *b_beam_inst_valid;   //!
   TBranch        *b_beam_inst_trigger;   //!
   TBranch        *b_beam_inst_nTracks;   //!
   TBranch        *b_beam_inst_nMomenta;   //!
   TBranch        *b_beam_inst_PDG_candidates;   //!
   TBranch        *b_beam_inst_TOF;   //!
   TBranch        *b_beam_inst_TOF_Chan;   //!
   TBranch        *b_isprimarytrack;   //!
   TBranch        *b_isprimaryshower;   //!
   TBranch        *b_beamPosx;   //!
   TBranch        *b_beamPosy;   //!
   TBranch        *b_beamPosz;   //!
   TBranch        *b_beamDirx;   //!
   TBranch        *b_beamDiry;   //!
   TBranch        *b_beamDirz;   //!
   TBranch        *b_beamMomentum;   //!
   TBranch        *b_tof;   //!
   TBranch        *b_tofs;   //!
   TBranch        *b_ch_tofs;   //!
   TBranch        *b_low_pressure_status;   //!
   TBranch        *b_high_pressure_status;   //!
   TBranch        *b_low_pressure;   //!
   TBranch        *b_high_pressure;   //!
   TBranch        *b_cosine_beam_primtrk;   //!
   TBranch        *b_primtrk_startx;   //!
   TBranch        *b_primtrk_starty;   //!
   TBranch        *b_primtrk_startz;   //!
   TBranch        *b_primtrk_endx;   //!
   TBranch        *b_primtrk_endy;   //!
   TBranch        *b_primtrk_endz;   //!
   TBranch        *b_primtrk_Dirx;   //!
   TBranch        *b_primtrk_Diry;   //!
   TBranch        *b_primtrk_Dirz;   //!
   TBranch        *b_primtrklen;   //!
   TBranch        *b_primtrkID;   //!
   TBranch        *b_primtrk_dqdx;   //!
   TBranch        *b_primtrk_dedx;   //!
   TBranch        *b_primtrk_resrange;   //!
   TBranch        *b_primtrk_range;   //!
   TBranch        *b_primtrk_hitx;   //!
   TBranch        *b_primtrk_hity;   //!
   TBranch        *b_primtrk_hitz;   //!
   TBranch        *b_primtrk_ke;   //!
   TBranch        *b_primtrk_pitch;   //!
   TBranch        *b_primtrk_wid;   //!
   TBranch        *b_primtrk_pt;   //!
   TBranch        *b_primtrk_calo_hit_index;   //!
   TBranch        *b_primtrk_ch;   //!
   TBranch        *b_pdg_code;   //!
   TBranch        *b_n_beamparticle;   //!
   TBranch        *b_n_daughter;   //!
   TBranch        *b_isPrimary;   //!
   TBranch        *b_pfp_self;   //!
   TBranch        *b_pfp_daughter;   //!
   TBranch        *b_primtrk_trktag;   //!
   TBranch        *b_Zintersection;   //!
   TBranch        *b_Yintersection;   //!
   TBranch        *b_Xintersection;   //!
   TBranch        *b_timeintersection;   //!

   ProtonRecombination_run5387(TTree *tree=0);
   virtual ~ProtonRecombination_run5387();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     IsNotEmptyCheck(Long64_t entry);
   virtual bool     IsStoppingProton(Long64_t entry, TGraph* &graph);
   virtual bool     PassesBeamCuts(Long64_t entry);
   virtual double   Calibrate_dQdx(Long64_t entry, int j, float prim_dqdx, float normalisation_factor, float calib_factor, TH1F* &X_correction_hist_2,TH2F* &YZ_correction_neg_hist_2);

   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ProtonRecombination_run5387_cxx
ProtonRecombination_run5387::ProtonRecombination_run5387(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("protonbeamana/beamana",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("protonbeamana/beamana","");
      //chain->Add("/dune/data2/users/hyliao/protodune/proton/v09_26_00/ana/protodune-sp_runset_5387_reco2_new/Beam_all.root/protonbeamana/beamana");
      chain->Add("/pnfs/dune/persistent/users/hyliao/protodune/proton/v09_26_00/ana/protodune-sp_runset_5387_reco2_new/Beam_all.root/protonbeamana/beamana");
      tree = chain;
#endif // SINGLE_TREE

   }
if (tree->InheritsFrom("TChain")) ((TChain*)tree)->LoadTree(0);
   Init(tree);
}

ProtonRecombination_run5387::~ProtonRecombination_run5387()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ProtonRecombination_run5387::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ProtonRecombination_run5387::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ProtonRecombination_run5387::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   beam_inst_PDG_candidates = 0;
   beam_inst_TOF = 0;
   beam_inst_TOF_Chan = 0;
   beamPosx = 0;
   beamPosy = 0;
   beamPosz = 0;
   beamDirx = 0;
   beamDiry = 0;
   beamDirz = 0;
   beamMomentum = 0;
   tofs = 0;
   ch_tofs = 0;
   primtrk_startx = 0;
   primtrk_starty = 0;
   primtrk_startz = 0;
   primtrk_endx = 0;
   primtrk_endy = 0;
   primtrk_endz = 0;
   primtrk_Dirx = 0;
   primtrk_Diry = 0;
   primtrk_Dirz = 0;
   primtrklen = 0;
   primtrkID = 0;
   primtrk_dqdx = 0;
   primtrk_dedx = 0;
   primtrk_resrange = 0;
   primtrk_range = 0;
   primtrk_hitx = 0;
   primtrk_hity = 0;
   primtrk_hitz = 0;
   primtrk_ke = 0;
   primtrk_pitch = 0;
   primtrk_wid = 0;
   primtrk_pt = 0;
   primtrk_calo_hit_index = 0;
   primtrk_ch = 0;
   pdg_code = 0;
   n_beamparticle = 0;
   n_daughter = 0;
   isPrimary = 0;
   pfp_self = 0;
   pfp_daughter = 0;
   primtrk_trktag = 0;
   Zintersection = 0;
   Yintersection = 0;
   Xintersection = 0;
   timeintersection = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   fChain->SetBranchAddress("Nactivefembs", Nactivefembs, &b_Nactivefembs);
   fChain->SetBranchAddress("beam_inst_valid", &beam_inst_valid, &b_beam_inst_valid);
   fChain->SetBranchAddress("beam_inst_trigger", &beam_inst_trigger, &b_beam_inst_trigger);
   fChain->SetBranchAddress("beam_inst_nTracks", &beam_inst_nTracks, &b_beam_inst_nTracks);
   fChain->SetBranchAddress("beam_inst_nMomenta", &beam_inst_nMomenta, &b_beam_inst_nMomenta);
   fChain->SetBranchAddress("beam_inst_PDG_candidates", &beam_inst_PDG_candidates, &b_beam_inst_PDG_candidates);
   fChain->SetBranchAddress("beam_inst_TOF", &beam_inst_TOF, &b_beam_inst_TOF);
   fChain->SetBranchAddress("beam_inst_TOF_Chan", &beam_inst_TOF_Chan, &b_beam_inst_TOF_Chan);
   fChain->SetBranchAddress("isprimarytrack", &isprimarytrack, &b_isprimarytrack);
   fChain->SetBranchAddress("isprimaryshower", &isprimaryshower, &b_isprimaryshower);
   fChain->SetBranchAddress("beamPosx", &beamPosx, &b_beamPosx);
   fChain->SetBranchAddress("beamPosy", &beamPosy, &b_beamPosy);
   fChain->SetBranchAddress("beamPosz", &beamPosz, &b_beamPosz);
   fChain->SetBranchAddress("beamDirx", &beamDirx, &b_beamDirx);
   fChain->SetBranchAddress("beamDiry", &beamDiry, &b_beamDiry);
   fChain->SetBranchAddress("beamDirz", &beamDirz, &b_beamDirz);
   fChain->SetBranchAddress("beamMomentum", &beamMomentum, &b_beamMomentum);
   fChain->SetBranchAddress("tof", &tof, &b_tof);
   fChain->SetBranchAddress("tofs", &tofs, &b_tofs);
   fChain->SetBranchAddress("ch_tofs", &ch_tofs, &b_ch_tofs);
   fChain->SetBranchAddress("low_pressure_status", &low_pressure_status, &b_low_pressure_status);
   fChain->SetBranchAddress("high_pressure_status", &high_pressure_status, &b_high_pressure_status);
   fChain->SetBranchAddress("low_pressure", &low_pressure, &b_low_pressure);
   fChain->SetBranchAddress("high_pressure", &high_pressure, &b_high_pressure);
   fChain->SetBranchAddress("cosine_beam_primtrk", &cosine_beam_primtrk, &b_cosine_beam_primtrk);
   fChain->SetBranchAddress("primtrk_startx", &primtrk_startx, &b_primtrk_startx);
   fChain->SetBranchAddress("primtrk_starty", &primtrk_starty, &b_primtrk_starty);
   fChain->SetBranchAddress("primtrk_startz", &primtrk_startz, &b_primtrk_startz);
   fChain->SetBranchAddress("primtrk_endx", &primtrk_endx, &b_primtrk_endx);
   fChain->SetBranchAddress("primtrk_endy", &primtrk_endy, &b_primtrk_endy);
   fChain->SetBranchAddress("primtrk_endz", &primtrk_endz, &b_primtrk_endz);
   fChain->SetBranchAddress("primtrk_Dirx", &primtrk_Dirx, &b_primtrk_Dirx);
   fChain->SetBranchAddress("primtrk_Diry", &primtrk_Diry, &b_primtrk_Diry);
   fChain->SetBranchAddress("primtrk_Dirz", &primtrk_Dirz, &b_primtrk_Dirz);
   fChain->SetBranchAddress("primtrklen", &primtrklen, &b_primtrklen);
   fChain->SetBranchAddress("primtrkID", &primtrkID, &b_primtrkID);
   fChain->SetBranchAddress("primtrk_dqdx", &primtrk_dqdx, &b_primtrk_dqdx);
   fChain->SetBranchAddress("primtrk_dedx", &primtrk_dedx, &b_primtrk_dedx);
   fChain->SetBranchAddress("primtrk_resrange", &primtrk_resrange, &b_primtrk_resrange);
   fChain->SetBranchAddress("primtrk_range", &primtrk_range, &b_primtrk_range);
   fChain->SetBranchAddress("primtrk_hitx", &primtrk_hitx, &b_primtrk_hitx);
   fChain->SetBranchAddress("primtrk_hity", &primtrk_hity, &b_primtrk_hity);
   fChain->SetBranchAddress("primtrk_hitz", &primtrk_hitz, &b_primtrk_hitz);
   fChain->SetBranchAddress("primtrk_ke", &primtrk_ke, &b_primtrk_ke);
   fChain->SetBranchAddress("primtrk_pitch", &primtrk_pitch, &b_primtrk_pitch);
   fChain->SetBranchAddress("primtrk_wid", &primtrk_wid, &b_primtrk_wid);
   fChain->SetBranchAddress("primtrk_pt", &primtrk_pt, &b_primtrk_pt);
   fChain->SetBranchAddress("primtrk_calo_hit_index", &primtrk_calo_hit_index, &b_primtrk_calo_hit_index);
   fChain->SetBranchAddress("primtrk_ch", &primtrk_ch, &b_primtrk_ch);
   fChain->SetBranchAddress("pdg_code", &pdg_code, &b_pdg_code);
   fChain->SetBranchAddress("n_beamparticle", &n_beamparticle, &b_n_beamparticle);
   fChain->SetBranchAddress("n_daughter", &n_daughter, &b_n_daughter);
   fChain->SetBranchAddress("isPrimary", &isPrimary, &b_isPrimary);
   fChain->SetBranchAddress("pfp_self", &pfp_self, &b_pfp_self);
   fChain->SetBranchAddress("pfp_daughter", &pfp_daughter, &b_pfp_daughter);
   fChain->SetBranchAddress("primtrk_trktag", &primtrk_trktag, &b_primtrk_trktag);
   fChain->SetBranchAddress("Zintersection", &Zintersection, &b_Zintersection);
   fChain->SetBranchAddress("Yintersection", &Yintersection, &b_Yintersection);
   fChain->SetBranchAddress("Xintersection", &Xintersection, &b_Xintersection);
   fChain->SetBranchAddress("timeintersection", &timeintersection, &b_timeintersection);
   Notify();
}

Bool_t ProtonRecombination_run5387::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ProtonRecombination_run5387::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ProtonRecombination_run5387::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ProtonDemo_run5387_cxx
