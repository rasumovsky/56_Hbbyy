//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 30 13:44:00 2016 by ROOT version 6.04/14
// from TTree mini/mini
// found on file: MGPy8_X400tohh_yybb.root
//////////////////////////////////////////////////////////

#ifndef MiniNtup_h
#define MiniNtup_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
using namespace std;
using std::vector;

class MiniNtup {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Bool_t          isMC;
   Int_t           runNumber;
   Int_t           eventNumber;
   Float_t         mu;
   Int_t           npv;
   Float_t         hgam_weightInitial;
   Float_t         hgam_weight;
   Float_t         hgam_weightSF;
   Float_t         mcWeight;
   Float_t         crossSectionBRfilterEff;
   Bool_t          isDalitz;
   Float_t         pileupWeight;
   Float_t         vertexWeight;
   Float_t         weightJvt;
   Float_t         yybb_weight;
   Int_t           yybb_btagCat;
   Float_t         m_yybb_constrained;
   Bool_t          hgam_isPassed;
   Int_t           hgam_cutList;
   Int_t           yybb_cutList;
   Float_t         m_yy;
   Float_t         truth_m_yy;
   Float_t         m_jj;
   Float_t         dy_jj;
   Float_t         dphi_jj;
   Float_t         dphi_yyjj;
   Float_t         dy_yyjj;
   Float_t         dR_yy;
   Float_t         dy_yy;
   Int_t           jet_n;
   Float_t         jet_pt[14];   //[jet_n]
   Float_t         jet_eta[14];   //[jet_n]
   Float_t         jet_phi[14];   //[jet_n]
   Float_t         jet_m[14];   //[jet_n]
   Float_t         jet_Jvt[14];   //[jet_n]
   Float_t         jet_SF_jvt[14];   //[jet_n]
   Int_t           jet_label[14];   //[jet_n]
   Float_t         jet_SF_MV2c10_FixedCutBEff_60[14];   //[jet_n]
   Float_t         jet_SF_MV2c10_FixedCutBEff_70[14];   //[jet_n]
   Float_t         jet_SF_MV2c10_FixedCutBEff_77[14];   //[jet_n]
   Float_t         jet_SF_MV2c10_FixedCutBEff_85[14];   //[jet_n]
   Bool_t          jet_MV2c10_FixedCutBEff_60[14];   //[jet_n]
   Bool_t          jet_MV2c10_FixedCutBEff_70[14];   //[jet_n]
   Bool_t          jet_MV2c10_FixedCutBEff_77[14];   //[jet_n]
   Bool_t          jet_MV2c10_FixedCutBEff_85[14];   //[jet_n]
   Float_t         met;
   Float_t         met_phi;
   Float_t         met_sumet;
   Int_t           photon_n;
   Float_t         photon_m[4];   //[photon_n]
   Float_t         photon_eta[4];   //[photon_n]
   Float_t         photon_phi[4];   //[photon_n]
   Float_t         photon_pt[4];   //[photon_n]
   Bool_t          photon_isTight[4];   //[photon_n]
   Bool_t          photon_iso_Loose[4];   //[photon_n]
   Bool_t          photon_iso_LooseCaloOnly[4];   //[photon_n]
   Bool_t          photon_iso_Tight[4];   //[photon_n]
   Bool_t          photon_iso_TightCaloOnly[4];   //[photon_n]
   Float_t         photon_ptcone20[4];   //[photon_n]
   Float_t         photon_topoEtcone40[4];   //[photon_n]
   Int_t           photon_truthOrigin[4];   //[photon_n]
   Int_t           photon_truthType[4];   //[photon_n]
   Float_t         photon_scaleFactor[4];   //[photon_n]
   Int_t           muon_n;
   Float_t         muon_pt;
   Float_t         muon_eta;
   Float_t         muon_phi;
   Float_t         muon_charge;
   Float_t         muon_scaleFactor;
   Float_t         muon_topoetcone20;
   Float_t         muon_ptvarcone20;
   Int_t           electron_n;
   Float_t         electron_pt;
   Float_t         electron_eta;
   Float_t         electron_phi;
   Float_t         electron_m;
   Float_t         electron_charge;
   Float_t         electron_scaleFactor;
   Float_t         electron_topoetcone20;
   Float_t         electron_ptvarcone20;
   Bool_t          electron_isTight;

   // List of branches
   TBranch        *b_isMC;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_hgam_weightInitial;   //!
   TBranch        *b_hgam_weight;   //!
   TBranch        *b_hgam_weightSF;   //!
   TBranch        *b_mcWeight;   //!
   TBranch        *b_crossSectionBRfilterEff;   //!
   TBranch        *b_isDalitz;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_vertexWeight;   //!
   TBranch        *b_weightJvt;   //!
   TBranch        *b_yybb_weight;   //!
   TBranch        *b_yybb_btagCat;   //!
   TBranch        *b_m_yybb_constrained;   //!
   TBranch        *b_hgam_isPassed;   //!
   TBranch        *b_hgam_cutList;   //!
   TBranch        *b_yybb_cutList;   //!
   TBranch        *b_m_yy;   //!
   TBranch        *b_truth_m_yy;   //!
   TBranch        *b_m_jj;   //!
   TBranch        *b_dy_jj;   //!
   TBranch        *b_dphi_jj;   //!
   TBranch        *b_dphi_yyjj;   //!
   TBranch        *b_dy_yyjj;   //!
   TBranch        *b_dR_yy;   //!
   TBranch        *b_dy_yy;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_jet_Jvt;   //!
   TBranch        *b_jet_SF_jvt;   //!
   TBranch        *b_jet_label;   //!
   TBranch        *b_jet_SF_MV2c10_FixedCutBEff_60;   //!
   TBranch        *b_jet_SF_MV2c10_FixedCutBEff_70;   //!
   TBranch        *b_jet_SF_MV2c10_FixedCutBEff_77;   //!
   TBranch        *b_jet_SF_MV2c10_FixedCutBEff_85;   //!
   TBranch        *b_jet_MV2c10_FixedCutBEff_60;   //!
   TBranch        *b_jet_MV2c10_FixedCutBEff_70;   //!
   TBranch        *b_jet_MV2c10_FixedCutBEff_77;   //!
   TBranch        *b_jet_MV2c10_FixedCutBEff_85;   //!
   TBranch        *b_met;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_photon_n;   //!
   TBranch        *b_photon_m;   //!
   TBranch        *b_photon_eta;   //!
   TBranch        *b_photon_phi;   //!
   TBranch        *b_photon_pt;   //!
   TBranch        *b_photon_isTight;   //!
   TBranch        *b_photon_iso_Loose;   //!
   TBranch        *b_photon_iso_LooseCaloOnly;   //!
   TBranch        *b_photon_iso_Tight;   //!
   TBranch        *b_photon_iso_TightCaloOnly;   //!
   TBranch        *b_photon_ptcone20;   //!
   TBranch        *b_photon_topoEtcone40;   //!
   TBranch        *b_photon_truthOrigin;   //!
   TBranch        *b_photon_truthType;   //!
   TBranch        *b_photon_scaleFactor;   //!
   TBranch        *b_muon_n;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_scaleFactor;   //!
   TBranch        *b_muon_topoetcone20;   //!
   TBranch        *b_muon_ptvarcone20;   //!
   TBranch        *b_electron_n;   //!
   TBranch        *b_electron_pt;   //!
   TBranch        *b_electron_eta;   //!
   TBranch        *b_electron_phi;   //!
   TBranch        *b_electron_m;   //!
   TBranch        *b_electron_charge;   //!
   TBranch        *b_electron_scaleFactor;   //!
   TBranch        *b_electron_topoetcone20;   //!
   TBranch        *b_electron_ptvarcone20;   //!
   TBranch        *b_electron_isTight;   //!

   MiniNtup(TTree *tree=0);
   virtual ~MiniNtup();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MiniNtup_cxx
MiniNtup::MiniNtup(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MGPy8_X400tohh_yybb.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MGPy8_X400tohh_yybb.root");
      }
      f->GetObject("mini",tree);

   }
   Init(tree);
}

MiniNtup::~MiniNtup()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MiniNtup::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MiniNtup::LoadTree(Long64_t entry)
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

void MiniNtup::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isMC", &isMC, &b_isMC);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("hgam_weightInitial", &hgam_weightInitial, &b_hgam_weightInitial);
   fChain->SetBranchAddress("hgam_weight", &hgam_weight, &b_hgam_weight);
   fChain->SetBranchAddress("hgam_weightSF", &hgam_weightSF, &b_hgam_weightSF);
   fChain->SetBranchAddress("mcWeight", &mcWeight, &b_mcWeight);
   fChain->SetBranchAddress("crossSectionBRfilterEff", &crossSectionBRfilterEff, &b_crossSectionBRfilterEff);
   fChain->SetBranchAddress("isDalitz", &isDalitz, &b_isDalitz);
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("vertexWeight", &vertexWeight, &b_vertexWeight);
   fChain->SetBranchAddress("weightJvt", &weightJvt, &b_weightJvt);
   fChain->SetBranchAddress("yybb_weight", &yybb_weight, &b_yybb_weight);
   fChain->SetBranchAddress("yybb_btagCat", &yybb_btagCat, &b_yybb_btagCat);
   fChain->SetBranchAddress("m_yybb_constrained", &m_yybb_constrained, &b_m_yybb_constrained);
   fChain->SetBranchAddress("hgam_isPassed", &hgam_isPassed, &b_hgam_isPassed);
   fChain->SetBranchAddress("hgam_cutList", &hgam_cutList, &b_hgam_cutList);
   fChain->SetBranchAddress("yybb_cutList", &yybb_cutList, &b_yybb_cutList);
   fChain->SetBranchAddress("m_yy", &m_yy, &b_m_yy);
   fChain->SetBranchAddress("truth_m_yy", &truth_m_yy, &b_truth_m_yy);
   fChain->SetBranchAddress("m_jj", &m_jj, &b_m_jj);
   fChain->SetBranchAddress("dy_jj", &dy_jj, &b_dy_jj);
   fChain->SetBranchAddress("dphi_jj", &dphi_jj, &b_dphi_jj);
   fChain->SetBranchAddress("dphi_yyjj", &dphi_yyjj, &b_dphi_yyjj);
   fChain->SetBranchAddress("dy_yyjj", &dy_yyjj, &b_dy_yyjj);
   fChain->SetBranchAddress("dR_yy", &dR_yy, &b_dR_yy);
   fChain->SetBranchAddress("dy_yy", &dy_yy, &b_dy_yy);
   fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_m", jet_m, &b_jet_m);
   fChain->SetBranchAddress("jet_Jvt", jet_Jvt, &b_jet_Jvt);
   fChain->SetBranchAddress("jet_SF_jvt", jet_SF_jvt, &b_jet_SF_jvt);
   fChain->SetBranchAddress("jet_label", jet_label, &b_jet_label);
   fChain->SetBranchAddress("jet_SF_MV2c10_FixedCutBEff_60", jet_SF_MV2c10_FixedCutBEff_60, &b_jet_SF_MV2c10_FixedCutBEff_60);
   fChain->SetBranchAddress("jet_SF_MV2c10_FixedCutBEff_70", jet_SF_MV2c10_FixedCutBEff_70, &b_jet_SF_MV2c10_FixedCutBEff_70);
   fChain->SetBranchAddress("jet_SF_MV2c10_FixedCutBEff_77", jet_SF_MV2c10_FixedCutBEff_77, &b_jet_SF_MV2c10_FixedCutBEff_77);
   fChain->SetBranchAddress("jet_SF_MV2c10_FixedCutBEff_85", jet_SF_MV2c10_FixedCutBEff_85, &b_jet_SF_MV2c10_FixedCutBEff_85);
   fChain->SetBranchAddress("jet_MV2c10_FixedCutBEff_60", jet_MV2c10_FixedCutBEff_60, &b_jet_MV2c10_FixedCutBEff_60);
   fChain->SetBranchAddress("jet_MV2c10_FixedCutBEff_70", jet_MV2c10_FixedCutBEff_70, &b_jet_MV2c10_FixedCutBEff_70);
   fChain->SetBranchAddress("jet_MV2c10_FixedCutBEff_77", jet_MV2c10_FixedCutBEff_77, &b_jet_MV2c10_FixedCutBEff_77);
   fChain->SetBranchAddress("jet_MV2c10_FixedCutBEff_85", jet_MV2c10_FixedCutBEff_85, &b_jet_MV2c10_FixedCutBEff_85);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_sumet", &met_sumet, &b_met_sumet);
   fChain->SetBranchAddress("photon_n", &photon_n, &b_photon_n);
   fChain->SetBranchAddress("photon_m", photon_m, &b_photon_m);
   fChain->SetBranchAddress("photon_eta", photon_eta, &b_photon_eta);
   fChain->SetBranchAddress("photon_phi", photon_phi, &b_photon_phi);
   fChain->SetBranchAddress("photon_pt", photon_pt, &b_photon_pt);
   fChain->SetBranchAddress("photon_isTight", photon_isTight, &b_photon_isTight);
   fChain->SetBranchAddress("photon_iso_Loose", photon_iso_Loose, &b_photon_iso_Loose);
   fChain->SetBranchAddress("photon_iso_LooseCaloOnly", photon_iso_LooseCaloOnly, &b_photon_iso_LooseCaloOnly);
   fChain->SetBranchAddress("photon_iso_Tight", photon_iso_Tight, &b_photon_iso_Tight);
   fChain->SetBranchAddress("photon_iso_TightCaloOnly", photon_iso_TightCaloOnly, &b_photon_iso_TightCaloOnly);
   fChain->SetBranchAddress("photon_ptcone20", photon_ptcone20, &b_photon_ptcone20);
   fChain->SetBranchAddress("photon_topoEtcone40", photon_topoEtcone40, &b_photon_topoEtcone40);
   fChain->SetBranchAddress("photon_truthOrigin", photon_truthOrigin, &b_photon_truthOrigin);
   fChain->SetBranchAddress("photon_truthType", photon_truthType, &b_photon_truthType);
   fChain->SetBranchAddress("photon_scaleFactor", photon_scaleFactor, &b_photon_scaleFactor);
   fChain->SetBranchAddress("muon_n", &muon_n, &b_muon_n);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_charge", &muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_scaleFactor", &muon_scaleFactor, &b_muon_scaleFactor);
   fChain->SetBranchAddress("muon_topoetcone20", &muon_topoetcone20, &b_muon_topoetcone20);
   fChain->SetBranchAddress("muon_ptvarcone20", &muon_ptvarcone20, &b_muon_ptvarcone20);
   fChain->SetBranchAddress("electron_n", &electron_n, &b_electron_n);
   fChain->SetBranchAddress("electron_pt", &electron_pt, &b_electron_pt);
   fChain->SetBranchAddress("electron_eta", &electron_eta, &b_electron_eta);
   fChain->SetBranchAddress("electron_phi", &electron_phi, &b_electron_phi);
   fChain->SetBranchAddress("electron_m", &electron_m, &b_electron_m);
   fChain->SetBranchAddress("electron_charge", &electron_charge, &b_electron_charge);
   fChain->SetBranchAddress("electron_scaleFactor", &electron_scaleFactor, &b_electron_scaleFactor);
   fChain->SetBranchAddress("electron_topoetcone20", &electron_topoetcone20, &b_electron_topoetcone20);
   fChain->SetBranchAddress("electron_ptvarcone20", &electron_ptvarcone20, &b_electron_ptvarcone20);
   fChain->SetBranchAddress("electron_isTight", &electron_isTight, &b_electron_isTight);
   Notify();
}

Bool_t MiniNtup::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MiniNtup::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MiniNtup::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MiniNtup_cxx
