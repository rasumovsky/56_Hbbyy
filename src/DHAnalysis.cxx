////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHAnalysis.cxx                                                      //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/06/2015                                                          //
//                                                                            //
//  This namespace stores all of the global information for the H->bb +yy     //
//  search with 13 TeV data in 2015. It also has all of the includes that are //
//  necessary for the analysis.                                               //
//                                                                            //
//  cateScheme = "note8TeV"                                                   //
//    cate 0: nonresonant CR (<2b)                                            //
//    cate 1: nonresonant SR (=2b)                                            //
//    cate 2: resonant CR (<2b)                                               //
//    cate 3: resonant SR (=2b)                                               //
//                                                                            //
//      TString sig  = "Signal";                                              //
//      TString bkg1 = "BkgSingleHiggs";                                      //
//      TString bkg2 = "BkgNonResonant";                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHAnalysis.h"

/**
   Get the analysis type based on the name of the signal.
   @param config - The configuration file for the analysis.
   @param sampleName - The name of the sample being used.
   @returns - The name of the analysis type.
*/
TString DHAnalysis::getAnalysisType(Config *config, TString sampleName) {
  std::vector<TString> anaTypes = config->getStrV("analysisTypes");
  for (int i_a = 0; i_a < (int)anaTypes.size(); i_a++) {
    if (sampleName.Contains(anaTypes[i_a])) return anaTypes[i_a];
  }
  std::cout << "DHAnalysis: Error! Analysis type unknown for " << sampleName 
	    << std::endl;
  exit(0);
}

/**
   Check if the given component should be incorporated in the fit in the
   specified category. 
   @param config - The configuration file for the analysis.
   @param cateName - The category name.
   @param component - The fit component.
   @returns - True iff the fit is included in the current category fit.
*/
bool DHAnalysis::cateHasComponent(Config *config, TString cateName,
				  TString component) {
  // A list of fit components in this category:
  std::vector<TString> fitComponentList
    = config->getStrV(Form("components%s", cateName.Data()));
  
  for (int i_c = 0; i_c < (int)fitComponentList.size(); i_c++) {
    if (component.EqualTo(fitComponentList[i_c])) return true;
  }
  return false;
}

/**
   Check if the sample is among those listed as a SM or DH signal. 
   @param config - The configuration file for the analysis.
   @param sampleName - the name of the sample being used.
   @returns - True iff the sample is a signal sample.
*/
bool DHAnalysis::isSMSample(Config *config, TString sampleName) {
  std::vector<TString> sigSMModes = config->getStrV("sigSMModes");
  for (int i_SM = 0; i_SM < (int)sigSMModes.size(); i_SM++) {
    if (sampleName.EqualTo(sigSMModes[i_SM])) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Check if the sample is among those listed as a di-Higgs signal. 
   @param config - The configuration file for the analysis.
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample is a DM signal sample.
*/
bool DHAnalysis::isDHSample(Config *config, TString sampleName) {
  std::vector<TString> sigDHModes = config->getStrV("sigDHModes");
  for (int i_DH = 0; i_DH < (int)sigDHModes.size(); i_DH++) {
    if (sampleName.EqualTo(sigDHModes[i_DH])) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Check if the sample is among those listed as a SM or DH signal. 
   @param config - The configuration file for the analysis.
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample is a SM or DM signal sample.
*/
bool DHAnalysis::isSignalSample(Config *config, TString sampleName) {
  return (isSMSample(config, sampleName) || isDHSample(config, sampleName));
}

/**
   Check whether a sample should be weighted.
   @param config - The configuration file for the analysis.
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample has associated event weights.
*/
bool DHAnalysis::isWeightedSample(Config *config, TString sampleName) {
  // First check if it is a SM or DM signal process:
  if (isSignalSample(config, sampleName)) return true;
  
  // Finally, check if it is one of the other MC processes:
  std::vector<TString> MCProcesses = config->getStrV("MCProcesses");
  for (int i_MC = 0; i_MC < (int)MCProcesses.size(); i_MC++) {
    if (sampleName.EqualTo(MCProcesses[i_MC])) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Get the mediator mass for the resonant analysis, else return zero. 
   @param sampleName - the name of the sample being used.
   @returns - the integer mass of the mediator, if resonant analysis. 
*/
int DHAnalysis::getMediatorMass(TString sampleName) {
  if (sampleName.Contains("ResMx") && !sampleName.Contains("NonRes")) {
    TString currSampleName = sampleName;
    currSampleName.ReplaceAll("ResMx", "");
    return currSampleName.Atoi();
  }
  else {
    std::cout << "DHAnalysis: Request for mediator mass failed. sampleName="
	      << sampleName << std::endl;
    return 0;
  }
}
