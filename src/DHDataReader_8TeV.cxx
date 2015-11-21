////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHDataReader.cxx                                                    //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 09/07/2015                                                          //
//                                                                            //
//  Retrieves data and MC sets for workspace creation.                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHDataReader.h"

/**
   -----------------------------------------------------------------------------
   Initialize the class for checking job statuses.
   @param configFile - The analysis configuration file.
   @param observable - The mass observable for the datasets.
*/
DHDataReader::DHDataReader(TString configFile, RooRealVar *observable) {

  // Load the analysis configuration file:
  m_config = new Config(configFile);
  
  if (observable) setMassObservable(observable);
  else {// construct diphoton invariant mass observable by default.
    RooRealVar *m_yy = new RooRealVar("m_yy", "m_yy", 
				      m_config->getNum("DHMyyRangeLo"),
				      m_config->getNum("DHMyyRangeHi"));
    setMassObservable(m_yy);
  }
}

/**
   -----------------------------------------------------------------------------
   Load the non-resonant data in the signal or control region into a RooDataSet.
   @param cateName - The name of the category (NonResCR or NonResSR):
*/
RooDataSet* DHDataReader::loadNonResData(TString cateName) {
  
  int cateOfInterest = 0;
  std::vector<TString> cateNameList = m_config->getStrV("cateNamesNonRes");
  for (cateOfInterest = 0; cateOfInterest < (int)cateNameList.size(); 
       cateOfInterest++) {
    if (cateName.EqualTo(cateNameList[cateOfInterest])) break;
  }
  
  // Create a RooDataSet to return:
  RooDataSet *currData = new RooDataSet(Form("obsData_%s", cateName.Data()),
					Form("obsData_%s", cateName.Data()),
					RooArgSet(*m_observable));
  // Load the TTree from file:
  TFile inputFile(m_config->getStr("nonResInput"));
  TTree *tree = (TTree*)inputFile.Get("nonres");
  double gg_mass; int catNonRes_idx;
  tree->SetBranchAddress("gg_mass", &gg_mass);
  tree->SetBranchAddress("catNonRes_idx", &catNonRes_idx);
  
  // Loop over TTree contents:
  for (int index = 0; index < tree->GetEntries(); index++) {
    tree->GetEntry(index);
    
    // Cut to make sure we are using the proper category:
    if (catNonRes_idx != cateOfInterest) continue;
    
    // Add points to the dataset:
    m_observable->setVal(gg_mass);
    currData->add(*m_observable);
  }
  return currData;
}

/**
   -----------------------------------------------------------------------------
   Load the resonant data in the signal or control region into a RooDataSet.
   @param cateName - The name of the category (ResonantCR or ResonantSR):
*/
RooDataSet* DHDataReader::loadResData(TString cateName) {
  
  // Create a RooDataSet to return:
  RooDataSet *currData = new RooDataSet(Form("obsData_%s", cateName.Data()),
					Form("obsData_%s", cateName.Data()),
					RooArgSet(*m_observable));

  if (cateName.EqualTo("ResonantCR")) {
    // Load the TTree from file:
    TFile inputFile(m_config->getStr("resData"));
    TTree *tree = (TTree*)inputFile.Get("photon");
    Float_t gg_jj_mass; Int_t cut_jj; Int_t blind; Float_t inv_mass;
    //tree->SetBranchAddress("gg_jj_mass", &gg_jj_mass);
    tree->SetBranchAddress("gg_jj_mass_constr", &gg_jj_mass);
    tree->SetBranchAddress("cut_jj", &cut_jj);
    tree->SetBranchAddress("blind", &blind);
    tree->SetBranchAddress("inv_mass", &inv_mass);
    
    // Loop over TTree contents:
    for (int index = 0; index < tree->GetEntries(); index++) {
      tree->GetEntry(index);
      
      // Apply a few cuts:
      if (cut_jj != 1 || blind == 1) continue;
      
      // Also apply the mass window cut:
      if (fabs(inv_mass - 125.4) > 3.2) continue;
      
      // Add points to the dataset:
      m_observable->setVal(gg_jj_mass/1000.0);
      currData->add(*m_observable);
    }
  }
  else if (cateName.EqualTo("ResonantSR")) {
    TFile inputFile(m_config->getStr("resInput"));
    TTree *tree = (TTree*)inputFile.Get("nonres");
    double mgg; int catNonRes_idx; double gg_bb_mass_constr;
    tree->SetBranchAddress("mgg", &mgg);
    tree->SetBranchAddress("catNonRes_idx", &catNonRes_idx);
    tree->SetBranchAddress("gg_bb_mass_constr", &gg_bb_mass_constr);
    // Loop over TTree contents:
    for (int index = 0; index < tree->GetEntries(); index++) {
      tree->GetEntry(index);
      if (!(fabs(mgg - 125.4) < 3.2)) continue;
      if (catNonRes_idx != 1) continue;
      m_observable->setVal(gg_bb_mass_constr);
      currData->add(*m_observable);
    }
  }
  else {
    std::cout << "DHDataReader: Undefined source for " << cateName << " data."
	      << std::endl;
    exit(0);
  }
  return currData;
}

/**
   -----------------------------------------------------------------------------
   Load the single Higgs MC in the resonant analysis signal region from file,
   and return a pointer to a RooDataSet object.
   @param cateName - The name of the category (should be ResonantSR):
*/
RooDataSet* DHDataReader::loadSingleHiggs(TString cateName) {
  RooRealVar *wt = new RooRealVar("wt", "wt", 1);
  // Create a RooDataSet to return:
  RooDataSet *currData = new RooDataSet(Form("resSHData_%s", cateName.Data()),
					Form("resSHData_%s", cateName.Data()),
					RooArgSet(*m_observable, *wt),
					RooFit::WeightVar(*wt));
  TFile inputFile(m_config->getStr("resSHData"));
  TH1F *histSH = (TH1F*)inputFile.Get("h");
  for (int i_b = 1; i_b < histSH->GetNbinsX(); i_b++) {
    double currWeight = histSH->GetBinContent(i_b);
    m_observable->setVal(histSH->GetBinCenter(i_b));
    wt->setVal(currWeight);
    currData->add(RooArgSet(*m_observable,*wt), currWeight);
  }
  return currData;
}

/**
   -----------------------------------------------------------------------------
   Provide a pointer to the RooRealVar used as the observable.
   @param observable - The RooRealVar used as the data observable.
*/
void DHDataReader::setMassObservable(RooRealVar *observable) {
  m_observable = observable;
}
