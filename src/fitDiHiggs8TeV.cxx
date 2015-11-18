////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  createSignalParameterization.cxx                                          //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 22/07/2015                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This main method provides a tool for performing individual fits to the    //
//  resonance Monte Carlo. Settings for the utility are provided in           //
//  singleSigFitExample.cfg.                                                  //
//                                                                            //
//  NOTE: When using other ntuples, the user will have to modify a few things.//
//  First, make sure the TTree branch types are properly assigned. Second,    //
//  use the proper units for mass. The SigParam tool always uses GeV. Third,  //
//  modify the luminosity and event weight. All of these items are marked in  //
//  the example below by: // USER MODIFICATION NECESSARY.                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamTools/HggTwoSidedCBPdf.h"
#include "HGamTools/SigParam.h"
#include "HGamTools/AtlasStyle.h"

#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"

/**
   -----------------------------------------------------------------------------
   The main method for this utility. Provide 1 argument - the location of the 
   config (.cfg) file, which should be stored in the data/ directory. The main()
   method runs over the samples provided, performs the fits requests, and gives
   comparisons of parameterized and non-parameterized fits. 
*/
int main(int argc, char *argv[])
{
  // Check that the config file location is provided.
  if (argc < 2) HG::fatal("No arguemnts provided");
  //HG::Config settings(TString(argv[1]));
  HG::Config *settings = new HG::Config(TString(argv[1]));

  // Print configuration for benefit of user:
  std::cout << "createSignalParameterization will run with parameters:"
	    << std::endl;
  settings->printDB();
  
  // Set the function type:
  TString function = settings->getStr("SignalFunctionalForm");
  
  // Check that output directory exists:
  TString outputDir = settings->getStr("OutputDir");
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Set the ATLAS Style for plots:
  SetAtlasStyle();
  
  // Instantiate SigParam class for individual & parameterized fits:
  SigParam *sps
    = new SigParam(settings->getStr("SampleName"), outputDir+"/Individual");
  sps->setLogYAxis(settings->getBool("MakeLogPlots"));
  //sps->doBinnedFit(settings->getBool("DoBinnedFit"), 1);
  
  // Use a hard-coded resonance mass if a branch is not included in MxAODs:
  int resMass = settings->getInt("ResonanceMass");
  double resMassToUse = (double)resMass;
  double massToUse; double weightToUse;
  double cateToUse = 0;
  
  //--------------------------------------//
  // Loop over the input file:
  TString fileName = settings->getStr("InputFile");
  std::ifstream inputFile;
  inputFile.open(fileName);
  while (!inputFile.eof()) {
    inputFile >> massToUse >> weightToUse;
    std::cout << massToUse << " " << weightToUse << std::endl;
    // Add the mass and weight values to the datasets for fitting:
    sps->addMassPoint(resMassToUse, cateToUse, massToUse, weightToUse);
  }
  inputFile.close();
  
  //--------------------------------------//
  // Now fit and plot the resonance shapes!
  std::cout << "createSignalParameterization: Start fitting and plotting!" 
	    << std::endl;
  
  // Loop over the analysis categories:
  if (sps->makeSingleResonance(resMassToUse, cateToUse, function)) {
    sps->plotSingleResonance(resMassToUse, cateToUse);
  }
  else {
    std::cout << "createSignalParameterization: Fit at mRes=" << resMassToUse
	      << ", cate=" << cateToUse << " did not converge :(" << std::endl;
  }
  sps->saveAll();
  
  //--------------------------------------//
  // Print some values for resonance mass of 125 GeV and category 0:
  std::cout << "\nfitResonantDiHiggs: Printing parameters for cate 0, M_Res = "
	    << resMass << " GeV" << std::endl;
  std::vector<TString> names = sps->getVariableNames(resMass, cateToUse);
  for (int i_n = 0; i_n < (int)names.size(); i_n++) {
    std::cout << "\t" << names[i_n] << " = " 
	      << sps->getParameterValue(names[i_n], resMass, cateToUse)
	      << std::endl;
  }
  return 0;
}
