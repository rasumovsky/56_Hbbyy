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
   -----------------------------------------------------------------------------
   Check if a particular fit component is used for the given category.
   @param cateName - The name of the category.
   @param component - The fit component.
   @returns - True iff specified component is found in the specified category.
*/
bool DHAnalysis::cateHasComponent(TString cateName, TString component) {
  std::vector<TString> componentList = getFitComponents(cateName);
  for (std::vector<TString>::iterator iter = componentList.begin();
       iter != componentList.end(); iter++) {
    if ((*iter).EqualTo(component)) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Convert the category index to the category name. 
   @param cateScheme - The name of the event categorization scheme.
   @param cateIndex - The category index.
   @returns - The category name.
*/
TString DHAnalysis::cateIndexToName(TString cateScheme, TString analysisType,
				    int cateIndex) {
  if (cateScheme.EqualTo("note8TeV")) {
    if (cateIndex == 0) {
      return Form("%s_%sCR", analysisType.Data(), cateScheme.Data());
    }
    else if (cateIndex == 1) {
      return Form("%s_%sSR", analysisType.Data(), cateScheme.Data());
    }
  }
  std::cout << "DHAnalysis: Error! no name for scheme=" << cateScheme 
	    << ", analysis=" << analysisType << ", category=" 
	    << cateIndex << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Convert the category index to the category name. 
   @param cateName - The name of the category.
   @returns - The category index.
*/
int DHAnalysis::cateNameToIndex(TString cateName) {
  if (cateName.Contains("note8TeV")) {
    if (cateName.Contains("CR")) return 0;
    else if (cateName.Contains("SR")) return 1;
  }
  std::cout << "DHAnalysis: Error! no index for " << cateName << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Determine the background PDF based on the category.
   @param cateScheme - the name of the event categorization scheme.
   @param cateIndex - the category index.
   @returns - the name of the background PDF.
*/
TString DHAnalysis::cateToBkgFunc(TString cateName) {
  //Possibilities are "Landau", "BernO1",..., "ExppolO1",... "Landau";
  TString result = "ExppolO1";
  if (cateName.Contains("note8TeV")) {
    if (cateName.Contains("NonRes")) return "ExppolO1";
    else if (cateName.Contains("Res")) return "Landau";
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Get a list of the fit components for a given category.
   @param cateName - The name of the category.
   @returns - Vector of template types (Signal, BkgSingleHiggs, BkgNonResonant).
*/
std::vector<TString> DHAnalysis::getFitComponents(TString cateName) {
  int cateIndex = cateNameToIndex(cateName);
  TString sig = "Signal";
  TString bkg1 = "BkgSingleHiggs";
  TString bkg2 = "BkgNonResonant";
  std::vector<TString> result; result.clear();
  if (cateName.Contains("note8TeV")) {
    if (cateIndex == 0) {
      result.push_back(bkg1);
      result.push_back(bkg2);
    }
    else if (cateIndex == 1) {
      result.push_back(sig);
      result.push_back(bkg1);
      result.push_back(bkg2);
    }
    else if (cateIndex == 2) {
      result.push_back(bkg2);
    }
    else if (cateIndex == 3) {
      result.push_back(sig);
      result.push_back(bkg1);
      result.push_back(bkg2);
    }
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Get number of categories for a particular categorization scheme.
   @param cateScheme - The name of the categorization scheme.
   @returns - The number of categories.
*/
int DHAnalysis::getNumCategories(TString cateScheme, TString analysisType) {
  if (cateScheme.EqualTo("note8TeV")) {
    if (analysisType.EqualTo("Res")) return 2;
    else if (analysisType.EqualTo("NonRes")) return 2;
  }
  std::cout << "DHAnalysis: Error! Categorization " << cateScheme << " and "
	    << analysisType << " has not been defined!." << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Check if the sample is among those listed as a SM signal. 
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample is a SM signal sample.
*/
bool DHAnalysis::isSMSample(TString sampleName) {
  for (int i_SM = 0; i_SM < nSMModes; i_SM++) {
    if (sampleName.Contains(sigSMModes[i_SM])) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Check if the sample is among those listed as a di-Higgs signal. 
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample is a DM signal sample.
*/
bool DHAnalysis::isDHSample(TString sampleName) {
  for (int i_DH = 0; i_DH < nDHModes; i_DH++) {
    if (sampleName.EqualTo(sigDHModes[i_DH])) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Check if the sample is among those listed as a SM or DM signal. 
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample is a SM or DM signal sample.
*/
bool DHAnalysis::isSignalSample(TString sampleName) {
  return (isSMSample(sampleName) || isDHSample(sampleName));
}

/**
   -----------------------------------------------------------------------------
   Check whether a sample should be weighted.
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample has associated event weights.
*/
bool DHAnalysis::isWeightedSample(TString sampleName) {
  // First check if it is a SM or DM signal process:
  if (isSignalSample(sampleName)) return true;
  
  // Finally, check if it is one of the other MC processes:
  for (int i_MC = 0; i_MC < nMCProcesses; i_MC++) {
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

/**
   -----------------------------------------------------------------------------
   Get the analysis type (Res, NonRes) using the name of the signal that is
   under investigation.
   @param signalName - the name of the signal being analyzed. 
   @returns "Res" or "NonRes" for resonant or non-resonant search. 
*/
TString DHAnalysis::getAnalysisType(TString signalName) {
  if (signalName.Contains("NonRes")) return "NonRes";
  else if (signalName.Contains("ResMx")) return "Res";
  else {
    std::cout << "DHAnalysis: analysis type not defined for " << signalName 
	      << std::endl;
    exit(0);
  }
}
