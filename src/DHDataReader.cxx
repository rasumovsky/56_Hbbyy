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
   @param observable - the mass observable for the datasets.
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
*/
RooDataSet* DHDataReader::loadResData(TString cateName) {
  // Create a RooDataSet to return:
  RooDataSet *currData = new RooDataSet(Form("obsData_%s", cateName.Data()),
					Form("obsData_%s", cateName.Data()),
					RooArgSet(*m_observable));
  // Load the TTree from file:
  TFile inputFile(m_config->getStr("resData"));
  TTree *tree = (TTree*)inputFile.Get("photon");
  double gg_jj_mass; bool cut_jj; bool blind;
  tree->SetBranchAddress("gg_jj_mass", &gg_jj_mass);
  tree->SetBranchAddress("cut_jj", &cut_jj);
  tree->SetBranchAddress("blind", &blind);
  
  // Loop over TTree contents:
  for (int index = 0; index < tree->GetEntries(); index++) {
    tree->GetEntry(index);
    
    // Apply a few cuts:
    if (!cut_jj || blind) continue;
    
    // Add points to the dataset:
    m_observable->setVal(gg_jj_mass);
    currData->add(*m_observable);
  }
  return currData;
}

/**
   -----------------------------------------------------------------------------
*/
void DHDataReader::setMassObservable(RooRealVar *observable) {
  m_observable = observable;
}
