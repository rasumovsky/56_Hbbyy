//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 24 10:28:49 2015 by ROOT version 5.34/09
// from TTree toy/toy
// found on file: toy_mu0_1987.root
//////////////////////////////////////////////////////////

#ifndef DMToyTree_h
#define DMToyTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;
using std::vector;

// Fixed size dimensions of array or collections stored in the TTree if any.

class DMToyTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           seed;
   Double_t        numEvents;
   Double_t        muDMVal;
   Bool_t          convergedMu0;
   Bool_t          convergedMu1;
   Bool_t          convergedMuFree;
   Double_t        nllMu0;
   Double_t        nllMu1;
   Double_t        nllMuFree;
   Double_t        llrL1L0;
   Double_t        llrL0Lfree;
   Double_t        llrL1Lfree;
   vector<double>  *numEventsPerCate;
   vector<string>  *namesNP;
   vector<double>  *valuesNPMu0;
   vector<double>  *valuesNPMu1;
   vector<double>  *valuesNPMuFree;
   vector<string>  *namesGlobs;
   vector<double>  *valuesGlobsMu1;
   vector<double>  *valuesGlobsMu0;
   vector<double>  *valuesGlobsMuFree;

   // List of branches
   TBranch        *b_seed;   //!
   TBranch        *b_numEvents;   //!
   TBranch        *b_muDMVal;   //!
   TBranch        *b_convergedMu0;   //!
   TBranch        *b_convergedMu1;   //!
   TBranch        *b_convergedMuFree;   //!
   TBranch        *b_nllMu0;   //!
   TBranch        *b_nllMu1;   //!
   TBranch        *b_nllMuFree;   //!
   TBranch        *b_llrL1L0;   //!
   TBranch        *b_llrL0Lfree;   //!
   TBranch        *b_llrL1Lfree;   //!
   TBranch        *b_numEventsPerCate;   //!
   TBranch        *b_namesNP;   //!
   TBranch        *b_valuesNPMu0;   //!
   TBranch        *b_valuesNPMu1;   //!
   TBranch        *b_valuesNPMuFree;   //!
   TBranch        *b_namesGlobs;   //!
   TBranch        *b_valuesGlobsMu1;   //!
   TBranch        *b_valuesGlobsMu0;   //!
   TBranch        *b_valuesGlobsMuFree;   //!

   DMToyTree(TTree *tree=0);
   virtual ~DMToyTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DMToyTree_cxx
DMToyTree::DMToyTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("toy_mu0_1987.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("toy_mu0_1987.root");
      }
      f->GetObject("toy",tree);

   }
   Init(tree);
}

DMToyTree::~DMToyTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DMToyTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DMToyTree::LoadTree(Long64_t entry)
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

void DMToyTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   numEventsPerCate = 0;
   namesNP = 0;
   valuesNPMu0 = 0;
   valuesNPMu1 = 0;
   valuesNPMuFree = 0;
   namesGlobs = 0;
   valuesGlobsMu1 = 0;
   valuesGlobsMu0 = 0;
   valuesGlobsMuFree = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("seed", &seed, &b_seed);
   fChain->SetBranchAddress("numEvents", &numEvents, &b_numEvents);
   fChain->SetBranchAddress("muDMVal", &muDMVal, &b_muDMVal);
   fChain->SetBranchAddress("convergedMu0", &convergedMu0, &b_convergedMu0);
   fChain->SetBranchAddress("convergedMu1", &convergedMu1, &b_convergedMu1);
   fChain->SetBranchAddress("convergedMuFree", &convergedMuFree, &b_convergedMuFree);
   fChain->SetBranchAddress("nllMu0", &nllMu0, &b_nllMu0);
   fChain->SetBranchAddress("nllMu1", &nllMu1, &b_nllMu1);
   fChain->SetBranchAddress("nllMuFree", &nllMuFree, &b_nllMuFree);
   fChain->SetBranchAddress("llrL1L0", &llrL1L0, &b_llrL1L0);
   fChain->SetBranchAddress("llrL0Lfree", &llrL0Lfree, &b_llrL0Lfree);
   fChain->SetBranchAddress("llrL1Lfree", &llrL1Lfree, &b_llrL1Lfree);
   fChain->SetBranchAddress("numEventsPerCate", &numEventsPerCate, &b_numEventsPerCate);
   fChain->SetBranchAddress("namesNP", &namesNP, &b_namesNP);
   fChain->SetBranchAddress("valuesNPMu0", &valuesNPMu0, &b_valuesNPMu0);
   fChain->SetBranchAddress("valuesNPMu1", &valuesNPMu1, &b_valuesNPMu1);
   fChain->SetBranchAddress("valuesNPMuFree", &valuesNPMuFree, &b_valuesNPMuFree);
   fChain->SetBranchAddress("namesGlobs", &namesGlobs, &b_namesGlobs);
   fChain->SetBranchAddress("valuesGlobsMu1", &valuesGlobsMu1, &b_valuesGlobsMu1);
   fChain->SetBranchAddress("valuesGlobsMu0", &valuesGlobsMu0, &b_valuesGlobsMu0);
   fChain->SetBranchAddress("valuesGlobsMuFree", &valuesGlobsMuFree, &b_valuesGlobsMuFree);
   Notify();
}

Bool_t DMToyTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DMToyTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DMToyTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DMToyTree_cxx
