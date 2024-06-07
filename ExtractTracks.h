//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 27 11:41:39 2020 by ROOT version 6.18/04
// from TChain protonmcnorw/PandoraBeam/
//////////////////////////////////////////////////////////

#ifndef ExtractTracks_h
#define ExtractTracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "string"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class ExtractTracks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   Double_t        timestamp;
   Int_t           Nactivefembs[5];
   Int_t           beamtrigger;
   Double_t        tof;
   Int_t           cerenkov1;
   Int_t           cerenkov2;
   Double_t        beamtrackMomentum;
   Double_t        beamtrackP[3];
   Double_t        beamtrackEnergy;
   Double_t        beamtrackPos[3];
   Double_t        beamtrackDir[3];
   Double_t        beamtrackTime;
   Int_t           beamtrackPdg;
   Int_t           beamtrackID;
   vector<double>  *beamtrk_x;
   vector<double>  *beamtrk_y;
   vector<double>  *beamtrk_z;
   vector<double>  *beamtrk_Px;
   vector<double>  *beamtrk_Py;
   vector<double>  *beamtrk_Pz;
   vector<double>  *beamtrk_Eng;
   vector<double>  *beamMomentum_spec;
   vector<double>  *beamPosx_spec;
   vector<double>  *beamPosy_spec;
   vector<double>  *beamPosz_spec;
   vector<double>  *beamDirx_spec;
   vector<double>  *beamDiry_spec;
   vector<double>  *beamDirz_spec;
   Double_t        ke_ff;
   vector<double>  *pos_ff;
   vector<double>  *primtrk_ke_true;
   vector<double>  *primtrk_ke_reco;
   string          *primtrk_end_g4process;
   vector<double>  *primtrk_range_true;
   Double_t        ke_ff_true;
   vector<vector<double> > *primtrk_hitx_true;
   vector<vector<double> > *primtrk_hity_true;
   vector<vector<double> > *primtrk_hitz_true;
   vector<vector<double> > *primtrk_trkid_true;
   vector<vector<double> > *primtrk_edept_true;
   vector<double>  *pfphit_peaktime_c;
   vector<double>  *pfphit_peaktime_u;
   vector<double>  *pfphit_peaktime_v;
   vector<double>  *pfphit_wireid_c;
   vector<double>  *pfphit_wireid_u;
   vector<double>  *pfphit_wireid_v;
   vector<vector<double> > *primtrk_true_x;
   vector<vector<double> > *primtrk_true_y;
   vector<vector<double> > *primtrk_true_z;
   vector<vector<double> > *primtrk_true_trkid;
   vector<vector<double> > *primtrk_true_edept;
   vector<vector<int> > *primtrk_true_wid;
   vector<vector<double> > *primtrj_true_x;
   vector<vector<double> > *primtrj_true_y;
   vector<vector<double> > *primtrj_true_z;
   vector<vector<double> > *primtrj_true_edept;
   Double_t        vertex[3];
   Double_t        secvertex[3];
   Int_t           isprimarytrack;
   Int_t           isprimaryshower;
   Double_t        primaryBDTScore;
   Int_t           primaryNHits;
   Double_t        primaryTheta;
   Double_t        primaryPhi;
   Double_t        primaryLength;
   Double_t        primaryMomentum;
   Double_t        primaryEndMomentum;
   Double_t        primaryEndPosition[3];
   Double_t        primaryStartPosition[3];
   Double_t        primaryEndDirection[3];
   Double_t        primaryStartDirection[3];
   Double_t        primaryOpeningAngle;
   Int_t           primaryID;
   Int_t           primaryShowerBestPlane;
   Int_t           primaryShowerEnergy;
   Int_t           primaryShowerMIPEnergy;
   Int_t           primaryShowerdEdx;
   Double_t        primaryMomentumByRangeProton;
   Double_t        primaryMomentumByRangeMuon;
   Double_t        primaryKineticEnergy[3];
   Double_t        primaryRange[3];
   Double_t        primaryT0;
   Int_t           primary_truth_TrackId;
   Int_t           primary_truth_Pdg;
   Int_t           truthpdg;
   Double_t        primary_truth_StartPosition[4];
   Double_t        primary_truth_StartPosition_MC[4];
   Double_t        primary_truth_EndPosition[4];
   Double_t        primary_truth_EndPosition_MC[4];
   Double_t        primary_truth_Momentum[4];
   Double_t        primary_truth_EndMomentum[4];
   Double_t        primary_truth_P;
   Double_t        primary_truth_Pt;
   Double_t        primary_truth_Mass;
   Double_t        primary_truth_Theta;
   Double_t        primary_truth_Phi;
   string          *primary_truth_EndProcess;
   Int_t           primary_truth_Isbeammatched;
   Int_t           primary_truth_NDaughters;
   vector<double>  *interactionX;
   vector<double>  *interactionY;
   vector<double>  *interactionZ;
   vector<double>  *interactionE;
   vector<string>  *interactionProcesslist;
   vector<double>  *interactionAngles;
   vector<double>  *Zintersection;
   vector<double>  *Yintersection;
   vector<double>  *Xintersection;
   vector<double>  *timeintersection;
   vector<double>  *interaction_wid_c;
   vector<double>  *interaction_tt_c;
   vector<double>  *interaction_wid_v;
   vector<double>  *interaction_tt_v;
   vector<double>  *interaction_wid_u;
   vector<double>  *interaction_tt_u;
   vector<float>   *x_c;
   vector<float>   *y_c;
   vector<float>   *z_c;
   vector<double>  *wid_c;
   vector<double>  *tt_c;
   vector<int>     *ch_c;
   vector<double>  *wirez_c;
   vector<double>  *pt_c;
   vector<double>  *q_c;
   vector<double>  *a_c;
   vector<int>     *wireno_c;
   vector<double>  *wid_v;
   vector<double>  *tt_v;
   vector<int>     *ch_v;
   vector<double>  *wirez_v;
   vector<double>  *pt_v;
   vector<double>  *q_v;
   vector<double>  *a_v;
   vector<int>     *wireno_v;
   vector<double>  *wid_u;
   vector<double>  *tt_u;
   vector<int>     *ch_u;
   vector<double>  *wirez_u;
   vector<double>  *pt_u;
   vector<double>  *q_u;
   vector<double>  *a_u;
   vector<int>     *wireno_u;
   Int_t           NDAUGHTERS;
   Int_t           isdaughtertrack[14];   //[NDAUGHTERS]
   Int_t           isdaughtershower[14];   //[NDAUGHTERS]
   Int_t           daughterNHits[14];   //[NDAUGHTERS]
   Double_t        daughterTheta[14];   //[NDAUGHTERS]
   Double_t        daughterPhi[14];   //[NDAUGHTERS]
   Double_t        daughterLength[14];   //[NDAUGHTERS]
   Double_t        daughterMomentum[14];   //[NDAUGHTERS]
   Double_t        daughterEndMomentum[14];   //[NDAUGHTERS]
   Double_t        daughterEndPosition[14][3];   //[NDAUGHTERS]
   Double_t        daughterStartPosition[14][3];   //[NDAUGHTERS]
   Double_t        daughterStartDirection[14][3];   //[NDAUGHTERS]
   Double_t        daughterEndDirection[14][3];   //[NDAUGHTERS]
   Double_t        daughterOpeningAngle[14];   //[NDAUGHTERS]
   Double_t        daughterShowerBestPlane[14];   //[NDAUGHTERS]
   Double_t        daughterShowerEnergy[14];   //[NDAUGHTERS]
   Double_t        daughterShowerMIPEnergy[14];   //[NDAUGHTERS]
   Double_t        daughterShowerdEdx[14];   //[NDAUGHTERS]
   Double_t        daughterMomentumByRangeProton[14];   //[NDAUGHTERS]
   Double_t        daughterMomentumByRangeMuon[14];   //[NDAUGHTERS]
   Double_t        daughterKineticEnergy[14][3];   //[NDAUGHTERS]
   Double_t        daughterRange[14][3];   //[NDAUGHTERS]
   Int_t           daughterID[14];   //[NDAUGHTERS]
   Int_t           daughter_truth_TrackId[14];   //[NDAUGHTERS]
   Int_t           daughter_truth_Pdg[14];   //[NDAUGHTERS]
   Double_t        daughter_truth_StartPosition[14][4];   //[NDAUGHTERS]
   Double_t        daughter_truth_EndPosition[14][4];   //[NDAUGHTERS]
   Double_t        daughter_truth_Momentum[14][4];   //[NDAUGHTERS]
   Double_t        daughter_truth_EndMomentum[14][4];   //[NDAUGHTERS]
   Double_t        daughter_truth_P[14];   //[NDAUGHTERS]
   Double_t        daughter_truth_Pt[14];   //[NDAUGHTERS]
   Double_t        daughter_truth_Mass[14];   //[NDAUGHTERS]
   Double_t        daughter_truth_Theta[14];   //[NDAUGHTERS]
   Double_t        daughter_truth_Phi[14];   //[NDAUGHTERS]
   Int_t           daughter_truth_Process[14];   //[NDAUGHTERS]
   vector<vector<double> > *primtrk_dqdx;
   vector<vector<double> > *primtrk_dedx;
   vector<vector<double> > *primtrk_resrange;
   vector<double>  *primtrk_range;
   vector<vector<double> > *primtrk_hitx;
   vector<vector<double> > *primtrk_hity;
   vector<vector<double> > *primtrk_hitz;
   vector<vector<double> > *primtrk_pitch;
   vector<vector<int> > *primtrk_wid;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_timestamp;   //!
   TBranch        *b_Nactivefembs;   //!
   TBranch        *b_beamtrigger;   //!
   TBranch        *b_tof;   //!
   TBranch        *b_cerenkov1;   //!
   TBranch        *b_cerenkov2;   //!
   TBranch        *b_beamtrackMomentum;   //!
   TBranch        *b_beamtrackP;   //!
   TBranch        *b_beamtrackEnergy;   //!
   TBranch        *b_beamtrackPos;   //!
   TBranch        *b_beamtrackDir;   //!
   TBranch        *b_beamtrackTime;   //!
   TBranch        *b_beamtrackPdg;   //!
   TBranch        *b_beamtrackID;   //!
   TBranch        *b_beamtrk_x;   //!
   TBranch        *b_beamtrk_y;   //!
   TBranch        *b_beamtrk_z;   //!
   TBranch        *b_beamtrk_Px;   //!
   TBranch        *b_beamtrk_Py;   //!
   TBranch        *b_beamtrk_Pz;   //!
   TBranch        *b_beamtrk_Eng;   //!
   TBranch        *b_beamMomentum_spec;   //!
   TBranch        *b_beamPosx_spec;   //!
   TBranch        *b_beamPosy_spec;   //!
   TBranch        *b_beamPosz_spec;   //!
   TBranch        *b_beamDirx_spec;   //!
   TBranch        *b_beamDiry_spec;   //!
   TBranch        *b_beamDirz_spec;   //!
   TBranch        *b_ke_ff;   //!
   TBranch        *b_pos_ff;   //!
   TBranch        *b_primtrk_ke_true;   //!
   TBranch        *b_primtrk_ke_reco;   //!
   TBranch        *b_primtrk_end_g4process;   //!
   TBranch        *b_primtrk_range_true;   //!
   TBranch        *b_ke_ff_true;   //!
   TBranch        *b_primtrk_hitx_true;   //!
   TBranch        *b_primtrk_hity_true;   //!
   TBranch        *b_primtrk_hitz_true;   //!
   TBranch        *b_primtrk_trkid_true;   //!
   TBranch        *b_primtrk_edept_true;   //!
   TBranch        *b_pfphit_peaktime_c;   //!
   TBranch        *b_pfphit_peaktime_u;   //!
   TBranch        *b_pfphit_peaktime_v;   //!
   TBranch        *b_pfphit_wireid_c;   //!
   TBranch        *b_pfphit_wireid_u;   //!
   TBranch        *b_pfphit_wireid_v;   //!
   TBranch        *b_primtrk_true_x;   //!
   TBranch        *b_primtrk_true_y;   //!
   TBranch        *b_primtrk_true_z;   //!
   TBranch        *b_primtrk_true_trkid;   //!
   TBranch        *b_primtrk_true_edept;   //!
   TBranch        *b_primtrk_true_wid;   //!
   TBranch        *b_primtrj_true_x;   //!
   TBranch        *b_primtrj_true_y;   //!
   TBranch        *b_primtrj_true_z;   //!
   TBranch        *b_primtrj_true_edept;   //!
   TBranch        *b_vertex;   //!
   TBranch        *b_secvertex;   //!
   TBranch        *b_isprimarytrack;   //!
   TBranch        *b_isprimaryshower;   //!
   TBranch        *b_primaryBDTScore;   //!
   TBranch        *b_fprimaryNHits;   //!
   TBranch        *b_primaryTheta;   //!
   TBranch        *b_primaryPhi;   //!
   TBranch        *b_primaryLength;   //!
   TBranch        *b_primaryMomentum;   //!
   TBranch        *b_primaryEndMomentum;   //!
   TBranch        *b_primaryEndPosition;   //!
   TBranch        *b_primaryStartPosition;   //!
   TBranch        *b_primaryEndDirection;   //!
   TBranch        *b_primaryStartDirection;   //!
   TBranch        *b_primaryOpeningAngle;   //!
   TBranch        *b_primaryID;   //!
   TBranch        *b_primaryShowerBestPlane;   //!
   TBranch        *b_primaryShowerEnergy;   //!
   TBranch        *b_primaryShowerMIPEnergy;   //!
   TBranch        *b_primaryShowerdEdx;   //!
   TBranch        *b_primaryMomentumByRangeProton;   //!
   TBranch        *b_primaryMomentumByRangeMuon;   //!
   TBranch        *b_primaryKineticEnergy;   //!
   TBranch        *b_primaryRange;   //!
   TBranch        *b_primaryT0;   //!
   TBranch        *b_primary_truth_TrackId;   //!
   TBranch        *b_primary_truth_Pdg;   //!
   TBranch        *b_truthpdg;   //!
   TBranch        *b_primary_truth_StartPosition;   //!
   TBranch        *b_primary_truth_StartPosition_MC;   //!
   TBranch        *b_primary_truth_EndPosition;   //!
   TBranch        *b_primary_truth_EndPosition_MC;   //!
   TBranch        *b_primary_truth_Momentum;   //!
   TBranch        *b_primary_truth_EndMomentum;   //!
   TBranch        *b_primary_truth_P;   //!
   TBranch        *b_primary_truth_Pt;   //!
   TBranch        *b_primary_truth_Mass;   //!
   TBranch        *b_primary_truth_Theta;   //!
   TBranch        *b_primary_truth_Phi;   //!
   TBranch        *b_primary_truth_EndProcess;   //!
   TBranch        *b_primary_truth_Isbeammatched;   //!
   TBranch        *b_primary_truth_NDaughters;   //!
   TBranch        *b_interactionX;   //!
   TBranch        *b_interactionY;   //!
   TBranch        *b_interactionZ;   //!
   TBranch        *b_interactionE;   //!
   TBranch        *b_interactionProcesslist;   //!
   TBranch        *b_interactionAngles;   //!
   TBranch        *b_Zintersection;   //!
   TBranch        *b_Yintersection;   //!
   TBranch        *b_Xintersection;   //!
   TBranch        *b_timeintersection;   //!
   TBranch        *b_interaction_wid_c;   //!
   TBranch        *b_interaction_tt_c;   //!
   TBranch        *b_interaction_wid_v;   //!
   TBranch        *b_interaction_tt_v;   //!
   TBranch        *b_interaction_wid_u;   //!
   TBranch        *b_interaction_tt_u;   //!
   TBranch        *b_x_c;   //!
   TBranch        *b_y_c;   //!
   TBranch        *b_z_c;   //!
   TBranch        *b_wid_c;   //!
   TBranch        *b_tt_c;   //!
   TBranch        *b_ch_c;   //!
   TBranch        *b_wirez_c;   //!
   TBranch        *b_pt_c;   //!
   TBranch        *b_q_c;   //!
   TBranch        *b_a_c;   //!
   TBranch        *b_wireno_c;   //!
   TBranch        *b_wid_v;   //!
   TBranch        *b_tt_v;   //!
   TBranch        *b_ch_v;   //!
   TBranch        *b_wirez_v;   //!
   TBranch        *b_pt_v;   //!
   TBranch        *b_q_v;   //!
   TBranch        *b_a_v;   //!
   TBranch        *b_wireno_v;   //!
   TBranch        *b_wid_u;   //!
   TBranch        *b_tt_u;   //!
   TBranch        *b_ch_u;   //!
   TBranch        *b_wirez_u;   //!
   TBranch        *b_pt_u;   //!
   TBranch        *b_q_u;   //!
   TBranch        *b_a_u;   //!
   TBranch        *b_wireno_u;   //!
   TBranch        *b_NDAUGHTERS;   //!
   TBranch        *b_isdaughtertrack;   //!
   TBranch        *b_isdaughtershower;   //!
   TBranch        *b_daughterNHits;   //!
   TBranch        *b_daughterTheta;   //!
   TBranch        *b_daughterPhi;   //!
   TBranch        *b_daughterLength;   //!
   TBranch        *b_daughterMomentum;   //!
   TBranch        *b_daughterEndMomentum;   //!
   TBranch        *b_daughterEndPosition;   //!
   TBranch        *b_daughterStartPosition;   //!
   TBranch        *b_daughterStartDirection;   //!
   TBranch        *b_daughterEndDirection;   //!
   TBranch        *b_daughterOpeningAngle;   //!
   TBranch        *b_daughterShowerBestPlane;   //!
   TBranch        *b_daughterShowerEnergy;   //!
   TBranch        *b_daughterShowerMIPEnergy;   //!
   TBranch        *b_daughterShowerdEdx;   //!
   TBranch        *b_daughterMomentumByRangeProton;   //!
   TBranch        *b_daughterMomentumByRangeMuon;   //!
   TBranch        *b_daughterKineticEnergy;   //!
   TBranch        *b_daughterRange;   //!
   TBranch        *b_daughterID;   //!
   TBranch        *b_daughter_truth_TrackId;   //!
   TBranch        *b_daughter_truth_Pdg;   //!
   TBranch        *b_daughter_truth_StartPosition;   //!
   TBranch        *b_daughter_truth_EndPosition;   //!
   TBranch        *b_daughter_truth_Momentum;   //!
   TBranch        *b_daughter_truth_EndMomentum;   //!
   TBranch        *b_daughter_truth_P;   //!
   TBranch        *b_daughter_truth_Pt;   //!
   TBranch        *b_daughter_truth_Mass;   //!
   TBranch        *b_daughter_truth_Theta;   //!
   TBranch        *b_daughter_truth_Phi;   //!
   TBranch        *b_daughter_truth_Process;   //!
   TBranch        *b_primtrk_dqdx;   //!
   TBranch        *b_primtrk_dedx;   //!
   TBranch        *b_primtrk_resrange;   //!
   TBranch        *b_primtrk_range;   //!
   TBranch        *b_primtrk_hitx;   //!
   TBranch        *b_primtrk_hity;   //!
   TBranch        *b_primtrk_hitz;   //!
   TBranch        *b_primtrk_pitch;   //!
   TBranch        *b_primtrk_wid;   //!

   ExtractTracks(TTree *tree=0);
   virtual ~ExtractTracks();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     IsNotEmptyCheck(Long64_t entry);
   virtual bool     IsStoppingProton(Long64_t entry, TGraph* &graph);
   virtual bool     PassesBeamCuts(Long64_t entry);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ExtractTracks_cxx
ExtractTracks::ExtractTracks(TTree *tree) : fChain(0) 
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
      f->GetObject("protonmcnorw/PandoraBeam",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("protonmcnorw/PandoraBeam","");
      chain->Add("/pnfs/dune/persistent/users/hyliao/protodune/proton/v09_04_01/ana/RITM0986948_MC_1GeV_reco_sce_datadriven_with_lifetime_corr_nog4rw/mcproton_all.root/protonmcnorw/PandoraBeam");
      // sce off version
      //chain->Add("/pnfs/dune/persistent/users/hyliao/protodune/proton/v09_04_01/ana/RITM0986948_MC_1GeV_reco_sce_off_nog4rw/mcproton_all.root/protonmcnorw/PandoraBeam");
      tree = chain;
#endif // SINGLE_TREE

   }
if (tree->InheritsFrom("TChain")) ((TChain*)tree)->LoadTree(0);
   Init(tree);
}

ExtractTracks::~ExtractTracks()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ExtractTracks::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ExtractTracks::LoadTree(Long64_t entry)
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

void ExtractTracks::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   beamtrk_x = 0;
   beamtrk_y = 0;
   beamtrk_z = 0;
   beamtrk_Px = 0;
   beamtrk_Py = 0;
   beamtrk_Pz = 0;
   beamtrk_Eng = 0;
   beamMomentum_spec = 0;
   beamPosx_spec = 0;
   beamPosy_spec = 0;
   beamPosz_spec = 0;
   beamDirx_spec = 0;
   beamDiry_spec = 0;
   beamDirz_spec = 0;
   pos_ff = 0;
   primtrk_ke_true = 0;
   primtrk_ke_reco = 0;
   primtrk_end_g4process = 0;
   primtrk_range_true = 0;
   primtrk_hitx_true = 0;
   primtrk_hity_true = 0;
   primtrk_hitz_true = 0;
   primtrk_trkid_true = 0;
   primtrk_edept_true = 0;
   pfphit_peaktime_c = 0;
   pfphit_peaktime_u = 0;
   pfphit_peaktime_v = 0;
   pfphit_wireid_c = 0;
   pfphit_wireid_u = 0;
   pfphit_wireid_v = 0;
   primtrk_true_x = 0;
   primtrk_true_y = 0;
   primtrk_true_z = 0;
   primtrk_true_trkid = 0;
   primtrk_true_edept = 0;
   primtrk_true_wid = 0;
   primtrj_true_x = 0;
   primtrj_true_y = 0;
   primtrj_true_z = 0;
   primtrj_true_edept = 0;
   primary_truth_EndProcess = 0;
   interactionX = 0;
   interactionY = 0;
   interactionZ = 0;
   interactionE = 0;
   interactionProcesslist = 0;
   interactionAngles = 0;
   Zintersection = 0;
   Yintersection = 0;
   Xintersection = 0;
   timeintersection = 0;
   interaction_wid_c = 0;
   interaction_tt_c = 0;
   interaction_wid_v = 0;
   interaction_tt_v = 0;
   interaction_wid_u = 0;
   interaction_tt_u = 0;
   x_c = 0;
   y_c = 0;
   z_c = 0;
   wid_c = 0;
   tt_c = 0;
   ch_c = 0;
   wirez_c = 0;
   pt_c = 0;
   q_c = 0;
   a_c = 0;
   wireno_c = 0;
   wid_v = 0;
   tt_v = 0;
   ch_v = 0;
   wirez_v = 0;
   pt_v = 0;
   q_v = 0;
   a_v = 0;
   wireno_v = 0;
   wid_u = 0;
   tt_u = 0;
   ch_u = 0;
   wirez_u = 0;
   pt_u = 0;
   q_u = 0;
   a_u = 0;
   wireno_u = 0;
   primtrk_dqdx = 0;
   primtrk_dedx = 0;
   primtrk_resrange = 0;
   primtrk_range = 0;
   primtrk_hitx = 0;
   primtrk_hity = 0;
   primtrk_hitz = 0;
   primtrk_pitch = 0;
   primtrk_wid = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   fChain->SetBranchAddress("Nactivefembs", Nactivefembs, &b_Nactivefembs);
   fChain->SetBranchAddress("beamtrigger", &beamtrigger, &b_beamtrigger);
   fChain->SetBranchAddress("tof", &tof, &b_tof);
   fChain->SetBranchAddress("cerenkov1", &cerenkov1, &b_cerenkov1);
   fChain->SetBranchAddress("cerenkov2", &cerenkov2, &b_cerenkov2);
   fChain->SetBranchAddress("beamtrackMomentum", &beamtrackMomentum, &b_beamtrackMomentum);
   fChain->SetBranchAddress("beamtrackP", beamtrackP, &b_beamtrackP);
   fChain->SetBranchAddress("beamtrackEnergy", &beamtrackEnergy, &b_beamtrackEnergy);
   fChain->SetBranchAddress("beamtrackPos", beamtrackPos, &b_beamtrackPos);
   fChain->SetBranchAddress("beamtrackDir", beamtrackDir, &b_beamtrackDir);
   fChain->SetBranchAddress("beamtrackTime", &beamtrackTime, &b_beamtrackTime);
   fChain->SetBranchAddress("beamtrackPdg", &beamtrackPdg, &b_beamtrackPdg);
   fChain->SetBranchAddress("beamtrackID", &beamtrackID, &b_beamtrackID);
   fChain->SetBranchAddress("beamtrk_x", &beamtrk_x, &b_beamtrk_x);
   fChain->SetBranchAddress("beamtrk_y", &beamtrk_y, &b_beamtrk_y);
   fChain->SetBranchAddress("beamtrk_z", &beamtrk_z, &b_beamtrk_z);
   fChain->SetBranchAddress("beamtrk_Px", &beamtrk_Px, &b_beamtrk_Px);
   fChain->SetBranchAddress("beamtrk_Py", &beamtrk_Py, &b_beamtrk_Py);
   fChain->SetBranchAddress("beamtrk_Pz", &beamtrk_Pz, &b_beamtrk_Pz);
   fChain->SetBranchAddress("beamtrk_Eng", &beamtrk_Eng, &b_beamtrk_Eng);
   fChain->SetBranchAddress("beamMomentum_spec", &beamMomentum_spec, &b_beamMomentum_spec);
   fChain->SetBranchAddress("beamPosx_spec", &beamPosx_spec, &b_beamPosx_spec);
   fChain->SetBranchAddress("beamPosy_spec", &beamPosy_spec, &b_beamPosy_spec);
   fChain->SetBranchAddress("beamPosz_spec", &beamPosz_spec, &b_beamPosz_spec);
   fChain->SetBranchAddress("beamDirx_spec", &beamDirx_spec, &b_beamDirx_spec);
   fChain->SetBranchAddress("beamDiry_spec", &beamDiry_spec, &b_beamDiry_spec);
   fChain->SetBranchAddress("beamDirz_spec", &beamDirz_spec, &b_beamDirz_spec);
   fChain->SetBranchAddress("ke_ff", &ke_ff, &b_ke_ff);
   fChain->SetBranchAddress("pos_ff", &pos_ff, &b_pos_ff);
   fChain->SetBranchAddress("primtrk_ke_true", &primtrk_ke_true, &b_primtrk_ke_true);
   fChain->SetBranchAddress("primtrk_ke_reco", &primtrk_ke_reco, &b_primtrk_ke_reco);
   fChain->SetBranchAddress("primtrk_end_g4process", &primtrk_end_g4process, &b_primtrk_end_g4process);
   fChain->SetBranchAddress("primtrk_range_true", &primtrk_range_true, &b_primtrk_range_true);
   fChain->SetBranchAddress("ke_ff_true", &ke_ff_true, &b_ke_ff_true);
   fChain->SetBranchAddress("primtrk_hitx_true", &primtrk_hitx_true, &b_primtrk_hitx_true);
   fChain->SetBranchAddress("primtrk_hity_true", &primtrk_hity_true, &b_primtrk_hity_true);
   fChain->SetBranchAddress("primtrk_hitz_true", &primtrk_hitz_true, &b_primtrk_hitz_true);
   fChain->SetBranchAddress("primtrk_trkid_true", &primtrk_trkid_true, &b_primtrk_trkid_true);
   fChain->SetBranchAddress("primtrk_edept_true", &primtrk_edept_true, &b_primtrk_edept_true);
   fChain->SetBranchAddress("pfphit_peaktime_c", &pfphit_peaktime_c, &b_pfphit_peaktime_c);
   fChain->SetBranchAddress("pfphit_peaktime_u", &pfphit_peaktime_u, &b_pfphit_peaktime_u);
   fChain->SetBranchAddress("pfphit_peaktime_v", &pfphit_peaktime_v, &b_pfphit_peaktime_v);
   fChain->SetBranchAddress("pfphit_wireid_c", &pfphit_wireid_c, &b_pfphit_wireid_c);
   fChain->SetBranchAddress("pfphit_wireid_u", &pfphit_wireid_u, &b_pfphit_wireid_u);
   fChain->SetBranchAddress("pfphit_wireid_v", &pfphit_wireid_v, &b_pfphit_wireid_v);
   fChain->SetBranchAddress("primtrk_true_x", &primtrk_true_x, &b_primtrk_true_x);
   fChain->SetBranchAddress("primtrk_true_y", &primtrk_true_y, &b_primtrk_true_y);
   fChain->SetBranchAddress("primtrk_true_z", &primtrk_true_z, &b_primtrk_true_z);
   fChain->SetBranchAddress("primtrk_true_trkid", &primtrk_true_trkid, &b_primtrk_true_trkid);
   fChain->SetBranchAddress("primtrk_true_edept", &primtrk_true_edept, &b_primtrk_true_edept);
   fChain->SetBranchAddress("primtrk_true_wid", &primtrk_true_wid, &b_primtrk_true_wid);
   fChain->SetBranchAddress("primtrj_true_x", &primtrj_true_x, &b_primtrj_true_x);
   fChain->SetBranchAddress("primtrj_true_y", &primtrj_true_y, &b_primtrj_true_y);
   fChain->SetBranchAddress("primtrj_true_z", &primtrj_true_z, &b_primtrj_true_z);
   fChain->SetBranchAddress("primtrj_true_edept", &primtrj_true_edept, &b_primtrj_true_edept);
   fChain->SetBranchAddress("vertex", vertex, &b_vertex);
   fChain->SetBranchAddress("secvertex", secvertex, &b_secvertex);
   fChain->SetBranchAddress("isprimarytrack", &isprimarytrack, &b_isprimarytrack);
   fChain->SetBranchAddress("isprimaryshower", &isprimaryshower, &b_isprimaryshower);
   fChain->SetBranchAddress("primaryBDTScore", &primaryBDTScore, &b_primaryBDTScore);
   fChain->SetBranchAddress("primaryNHits", &primaryNHits, &b_fprimaryNHits);
   fChain->SetBranchAddress("primaryTheta", &primaryTheta, &b_primaryTheta);
   fChain->SetBranchAddress("primaryPhi", &primaryPhi, &b_primaryPhi);
   fChain->SetBranchAddress("primaryLength", &primaryLength, &b_primaryLength);
   fChain->SetBranchAddress("primaryMomentum", &primaryMomentum, &b_primaryMomentum);
   fChain->SetBranchAddress("primaryEndMomentum", &primaryEndMomentum, &b_primaryEndMomentum);
   fChain->SetBranchAddress("primaryEndPosition", primaryEndPosition, &b_primaryEndPosition);
   fChain->SetBranchAddress("primaryStartPosition", primaryStartPosition, &b_primaryStartPosition);
   fChain->SetBranchAddress("primaryEndDirection", primaryEndDirection, &b_primaryEndDirection);
   fChain->SetBranchAddress("primaryStartDirection", primaryStartDirection, &b_primaryStartDirection);
   fChain->SetBranchAddress("primaryOpeningAngle", &primaryOpeningAngle, &b_primaryOpeningAngle);
   fChain->SetBranchAddress("primaryID", &primaryID, &b_primaryID);
   fChain->SetBranchAddress("primaryShowerBestPlane", &primaryShowerBestPlane, &b_primaryShowerBestPlane);
   fChain->SetBranchAddress("primaryShowerEnergy", &primaryShowerEnergy, &b_primaryShowerEnergy);
   fChain->SetBranchAddress("primaryShowerMIPEnergy", &primaryShowerMIPEnergy, &b_primaryShowerMIPEnergy);
   fChain->SetBranchAddress("primaryShowerdEdx", &primaryShowerdEdx, &b_primaryShowerdEdx);
   fChain->SetBranchAddress("primaryMomentumByRangeProton", &primaryMomentumByRangeProton, &b_primaryMomentumByRangeProton);
   fChain->SetBranchAddress("primaryMomentumByRangeMuon", &primaryMomentumByRangeMuon, &b_primaryMomentumByRangeMuon);
   fChain->SetBranchAddress("primaryKineticEnergy", primaryKineticEnergy, &b_primaryKineticEnergy);
   fChain->SetBranchAddress("primaryRange", primaryRange, &b_primaryRange);
   fChain->SetBranchAddress("primaryT0", &primaryT0, &b_primaryT0);
   fChain->SetBranchAddress("primary_truth_TrackId", &primary_truth_TrackId, &b_primary_truth_TrackId);
   fChain->SetBranchAddress("primary_truth_Pdg", &primary_truth_Pdg, &b_primary_truth_Pdg);
   fChain->SetBranchAddress("truthpdg", &truthpdg, &b_truthpdg);
   fChain->SetBranchAddress("primary_truth_StartPosition", primary_truth_StartPosition, &b_primary_truth_StartPosition);
   fChain->SetBranchAddress("primary_truth_StartPosition_MC", primary_truth_StartPosition_MC, &b_primary_truth_StartPosition_MC);
   fChain->SetBranchAddress("primary_truth_EndPosition", primary_truth_EndPosition, &b_primary_truth_EndPosition);
   fChain->SetBranchAddress("primary_truth_EndPosition_MC", primary_truth_EndPosition_MC, &b_primary_truth_EndPosition_MC);
   fChain->SetBranchAddress("primary_truth_Momentum", primary_truth_Momentum, &b_primary_truth_Momentum);
   fChain->SetBranchAddress("primary_truth_EndMomentum", primary_truth_EndMomentum, &b_primary_truth_EndMomentum);
   fChain->SetBranchAddress("primary_truth_P", &primary_truth_P, &b_primary_truth_P);
   fChain->SetBranchAddress("primary_truth_Pt", &primary_truth_Pt, &b_primary_truth_Pt);
   fChain->SetBranchAddress("primary_truth_Mass", &primary_truth_Mass, &b_primary_truth_Mass);
   fChain->SetBranchAddress("primary_truth_Theta", &primary_truth_Theta, &b_primary_truth_Theta);
   fChain->SetBranchAddress("primary_truth_Phi", &primary_truth_Phi, &b_primary_truth_Phi);
   fChain->SetBranchAddress("primary_truth_EndProcess", &primary_truth_EndProcess, &b_primary_truth_EndProcess);
   fChain->SetBranchAddress("primary_truth_Isbeammatched", &primary_truth_Isbeammatched, &b_primary_truth_Isbeammatched);
   fChain->SetBranchAddress("primary_truth_NDaughters", &primary_truth_NDaughters, &b_primary_truth_NDaughters);
   fChain->SetBranchAddress("interactionX", &interactionX, &b_interactionX);
   fChain->SetBranchAddress("interactionY", &interactionY, &b_interactionY);
   fChain->SetBranchAddress("interactionZ", &interactionZ, &b_interactionZ);
   fChain->SetBranchAddress("interactionE", &interactionE, &b_interactionE);
   fChain->SetBranchAddress("interactionProcesslist", &interactionProcesslist, &b_interactionProcesslist);
   fChain->SetBranchAddress("interactionAngles", &interactionAngles, &b_interactionAngles);
   fChain->SetBranchAddress("Zintersection", &Zintersection, &b_Zintersection);
   fChain->SetBranchAddress("Yintersection", &Yintersection, &b_Yintersection);
   fChain->SetBranchAddress("Xintersection", &Xintersection, &b_Xintersection);
   fChain->SetBranchAddress("timeintersection", &timeintersection, &b_timeintersection);
   fChain->SetBranchAddress("interaction_wid_c", &interaction_wid_c, &b_interaction_wid_c);
   fChain->SetBranchAddress("interaction_tt_c", &interaction_tt_c, &b_interaction_tt_c);
   fChain->SetBranchAddress("interaction_wid_v", &interaction_wid_v, &b_interaction_wid_v);
   fChain->SetBranchAddress("interaction_tt_v", &interaction_tt_v, &b_interaction_tt_v);
   fChain->SetBranchAddress("interaction_wid_u", &interaction_wid_u, &b_interaction_wid_u);
   fChain->SetBranchAddress("interaction_tt_u", &interaction_tt_u, &b_interaction_tt_u);
   fChain->SetBranchAddress("x_c", &x_c, &b_x_c);
   fChain->SetBranchAddress("y_c", &y_c, &b_y_c);
   fChain->SetBranchAddress("z_c", &z_c, &b_z_c);
   fChain->SetBranchAddress("wid_c", &wid_c, &b_wid_c);
   fChain->SetBranchAddress("tt_c", &tt_c, &b_tt_c);
   fChain->SetBranchAddress("ch_c", &ch_c, &b_ch_c);
   fChain->SetBranchAddress("wirez_c", &wirez_c, &b_wirez_c);
   fChain->SetBranchAddress("pt_c", &pt_c, &b_pt_c);
   fChain->SetBranchAddress("q_c", &q_c, &b_q_c);
   fChain->SetBranchAddress("a_c", &a_c, &b_a_c);
   fChain->SetBranchAddress("wireno_c", &wireno_c, &b_wireno_c);
   fChain->SetBranchAddress("wid_v", &wid_v, &b_wid_v);
   fChain->SetBranchAddress("tt_v", &tt_v, &b_tt_v);
   fChain->SetBranchAddress("ch_v", &ch_v, &b_ch_v);
   fChain->SetBranchAddress("wirez_v", &wirez_v, &b_wirez_v);
   fChain->SetBranchAddress("pt_v", &pt_v, &b_pt_v);
   fChain->SetBranchAddress("q_v", &q_v, &b_q_v);
   fChain->SetBranchAddress("a_v", &a_v, &b_a_v);
   fChain->SetBranchAddress("wireno_v", &wireno_v, &b_wireno_v);
   fChain->SetBranchAddress("wid_u", &wid_u, &b_wid_u);
   fChain->SetBranchAddress("tt_u", &tt_u, &b_tt_u);
   fChain->SetBranchAddress("ch_u", &ch_u, &b_ch_u);
   fChain->SetBranchAddress("wirez_u", &wirez_u, &b_wirez_u);
   fChain->SetBranchAddress("pt_u", &pt_u, &b_pt_u);
   fChain->SetBranchAddress("q_u", &q_u, &b_q_u);
   fChain->SetBranchAddress("a_u", &a_u, &b_a_u);
   fChain->SetBranchAddress("wireno_u", &wireno_u, &b_wireno_u);
   fChain->SetBranchAddress("NDAUGHTERS", &NDAUGHTERS, &b_NDAUGHTERS);
   fChain->SetBranchAddress("isdaughtertrack", isdaughtertrack, &b_isdaughtertrack);
   fChain->SetBranchAddress("isdaughtershower", isdaughtershower, &b_isdaughtershower);
   fChain->SetBranchAddress("daughterNHits", daughterNHits, &b_daughterNHits);
   fChain->SetBranchAddress("daughterTheta", daughterTheta, &b_daughterTheta);
   fChain->SetBranchAddress("daughterPhi", daughterPhi, &b_daughterPhi);
   fChain->SetBranchAddress("daughterLength", daughterLength, &b_daughterLength);
   fChain->SetBranchAddress("daughterMomentum", daughterMomentum, &b_daughterMomentum);
   fChain->SetBranchAddress("daughterEndMomentum", daughterEndMomentum, &b_daughterEndMomentum);
   fChain->SetBranchAddress("daughterEndPosition", daughterEndPosition, &b_daughterEndPosition);
   fChain->SetBranchAddress("daughterStartPosition", daughterStartPosition, &b_daughterStartPosition);
   fChain->SetBranchAddress("daughterStartDirection", daughterStartDirection, &b_daughterStartDirection);
   fChain->SetBranchAddress("daughterEndDirection", daughterEndDirection, &b_daughterEndDirection);
   fChain->SetBranchAddress("daughterOpeningAngle", daughterOpeningAngle, &b_daughterOpeningAngle);
   fChain->SetBranchAddress("daughterShowerBestPlane", daughterShowerBestPlane, &b_daughterShowerBestPlane);
   fChain->SetBranchAddress("daughterShowerEnergy", daughterShowerEnergy, &b_daughterShowerEnergy);
   fChain->SetBranchAddress("daughterShowerMIPEnergy", daughterShowerMIPEnergy, &b_daughterShowerMIPEnergy);
   fChain->SetBranchAddress("daughterShowerdEdx", daughterShowerdEdx, &b_daughterShowerdEdx);
   fChain->SetBranchAddress("daughterMomentumByRangeProton", daughterMomentumByRangeProton, &b_daughterMomentumByRangeProton);
   fChain->SetBranchAddress("daughterMomentumByRangeMuon", daughterMomentumByRangeMuon, &b_daughterMomentumByRangeMuon);
   fChain->SetBranchAddress("daughterKineticEnergy", daughterKineticEnergy, &b_daughterKineticEnergy);
   fChain->SetBranchAddress("daughterRange", daughterRange, &b_daughterRange);
   fChain->SetBranchAddress("daughterID", daughterID, &b_daughterID);
   fChain->SetBranchAddress("daughter_truth_TrackId", daughter_truth_TrackId, &b_daughter_truth_TrackId);
   fChain->SetBranchAddress("daughter_truth_Pdg", daughter_truth_Pdg, &b_daughter_truth_Pdg);
   fChain->SetBranchAddress("daughter_truth_StartPosition", daughter_truth_StartPosition, &b_daughter_truth_StartPosition);
   fChain->SetBranchAddress("daughter_truth_EndPosition", daughter_truth_EndPosition, &b_daughter_truth_EndPosition);
   fChain->SetBranchAddress("daughter_truth_Momentum", daughter_truth_Momentum, &b_daughter_truth_Momentum);
   fChain->SetBranchAddress("daughter_truth_EndMomentum", daughter_truth_EndMomentum, &b_daughter_truth_EndMomentum);
   fChain->SetBranchAddress("daughter_truth_P", daughter_truth_P, &b_daughter_truth_P);
   fChain->SetBranchAddress("daughter_truth_Pt", daughter_truth_Pt, &b_daughter_truth_Pt);
   fChain->SetBranchAddress("daughter_truth_Mass", daughter_truth_Mass, &b_daughter_truth_Mass);
   fChain->SetBranchAddress("daughter_truth_Theta", daughter_truth_Theta, &b_daughter_truth_Theta);
   fChain->SetBranchAddress("daughter_truth_Phi", daughter_truth_Phi, &b_daughter_truth_Phi);
   fChain->SetBranchAddress("daughter_truth_Process", daughter_truth_Process, &b_daughter_truth_Process);
   fChain->SetBranchAddress("primtrk_dqdx", &primtrk_dqdx, &b_primtrk_dqdx);
   fChain->SetBranchAddress("primtrk_dedx", &primtrk_dedx, &b_primtrk_dedx);
   fChain->SetBranchAddress("primtrk_resrange", &primtrk_resrange, &b_primtrk_resrange);
   fChain->SetBranchAddress("primtrk_range", &primtrk_range, &b_primtrk_range);
   fChain->SetBranchAddress("primtrk_hitx", &primtrk_hitx, &b_primtrk_hitx);
   fChain->SetBranchAddress("primtrk_hity", &primtrk_hity, &b_primtrk_hity);
   fChain->SetBranchAddress("primtrk_hitz", &primtrk_hitz, &b_primtrk_hitz);
   fChain->SetBranchAddress("primtrk_pitch", &primtrk_pitch, &b_primtrk_pitch);
   fChain->SetBranchAddress("primtrk_wid", &primtrk_wid, &b_primtrk_wid);
   Notify();
}

Bool_t ExtractTracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ExtractTracks::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}



Int_t ExtractTracks::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ProtonRecombination_cxx
