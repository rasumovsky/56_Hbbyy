////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHDataReader.cxx                                                    //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 01/13/2016                                                          //
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
DHDataReader::DHDataReader(TString configFile) {
  // Load the analysis configuration file:
  m_config = new Config(configFile);
  m_storedMxAODTrees.clear();
}

/**
   -----------------------------------------------------------------------------
   Load the non-resonant data in the signal or control region into a RooDataSet.
   @param sampleName - The name of the MxAOD sample.
*/
void DHDataReader::loadMxAOD(TString sampleName) {
  TString fileName = Form("%s/%s.root",
			  (m_config->getStr("MxAODDirectory")).Data(), 
			  sampleName.Data());
  TFile inputFile(fileName);
  if (inputFile.IsOpen()) {
    m_storedMxAODTrees[sampleName] = (TTree*)inputFile.Get("CollectionTree");
  }
  else {
    std::cout << "DHDataReader: ERROR Opening file " << sampleName << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   @param sampleName - The name of the MxAOD sample.
   @param cateName - The name of the category.
*/
RooDataSet* DHDataReader::getDataSet(TString sampleName, TString cateName) {
  
  // Get the name of the observable from the config file:
  TString obsForm = m_config->getStr(Form("OBS_%s",cateName.Data()));
  TString obsName = ""; double obsMin = 0.0; double obsMax = 0.0;
  varData(obsForm, obsName, obsMin, obsMax);
  RooRealVar obs(obsName, obsName, obsMin, obsMax);
  RooRealVar wt("wt", "wt", 1.0);
  RooDataSet* data 
    = new RooDataSet(Form("data_%s_%s",sampleName.Data(),cateName.Data()),
		     Form("data_%s_%s",sampleName.Data(),cateName.Data()),
		     RooArgSet(obs,wt), RooFit::WeightVar(wt));
  
  // Load the MxAOD file:
  loadMxAOD(sampleName);
  
  // Need to distinguish resonant, nonresonant, etc...

  // Load relevant branches:
  float b_weight; float b_mass; TString b_category;
  m_storedMxAODTrees[sampleName]->SetBranchAddress("weightname", &b_weight);
  m_storedMxAODTrees[sampleName]->SetBranchAddress("massname", &b_mass);
  m_storedMxAODTrees[sampleName]->SetBranchAddress("catename", &b_category);
  
  // Loop over the contents of the TTree:
  for (int index = 0; index < m_storedMxAODTrees[sampleName]->GetEntries(); 
       index++) {
    m_storedMxAODTrees[sampleName]->GetEntry(index);
    
    // Check that we are in the proper category:
    
    // Perform necesssary cuts:
    
    // Add point to the dataset:
    obs.setVal(b_mass);
    wt.setVal(b_weight);
    data->add(RooArgSet(obs,wt), b_weight);
  }
  return data;
}

/*
   -----------------------------------------------------------------------------
   Load the resonant data in the signal or control region into a RooDataSet.
   @param cateName - The name of the category (ResonantCR or ResonantSR):

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
*/

/**
   -----------------------------------------------------------------------------
   Takes in a variable declaration such as "name[1,0,5]" and returns "name".
   @param varForm - The form of the variable declared.
   @param varName - The extracted name of the variable.
   @param varMin - The extracted minimum of the variable.
   @param varMax - The extracted maximum of the variable.
   @return - The name of the variable without the rest of the expression.
*/
void DHDataReader::varData(TString varForm, TString& varName, double& varMin,
			   double& varMax) {
  varName = varForm;
  varName.Remove(varName.First("["));
  varForm.ReplaceAll("[","");
  varForm.ReplaceAll("]","");
  TString strVarMin = varForm;
  strVarMin.Remove(strVarMin.First(","));
  varMin = strVarMin.Atof();
  TString strVarMax = varForm;
  strVarMax.Remove(0,strVarMax.First(",")+1);
  varMax = strVarMax.Atof();
}
