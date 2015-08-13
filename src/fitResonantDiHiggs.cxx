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
   Checks whether the given mass value is contained in the vector, and adds it
   if it is not found.
   @param currList - The current vector of mass points.
   @param newMass - The mass point to check for membership in the list.
   @returns - A list of unique resonance mass points.
*/
std::vector<double> checkMResList(std::vector<double> currList,double newMass) {
  bool matched = false;
  for (int i_m = 0; i_m < (int)currList.size(); i_m++) {
    if (fabs(currList[i_m] - newMass) <= 0.1) {
      matched = true;
      break;
    }
  }
  if (!matched) currList.push_back(newMass);
  return currList;
}

/**
   -----------------------------------------------------------------------------
   Prints a progress bar to screen to provide elapsed time and remaining time
   information to the user. This is useful when processing large datasets. 
   @param index - The current event index.
   @param total - The total number of events.
*/
void PrintProgressBar(int index, int total) {
  if (index%10000 == 0) {
    TString print_bar = " [";
    for (int bar = 0; bar < 20; bar++) {
      double current_fraction = double(bar) / 20.0;
      if (double(index)/double(total) > current_fraction) print_bar.Append("/");
      else print_bar.Append(".");
    }
    print_bar.Append("] ");
    double percent = 100.0 * (double(index) / double(total));
    TString text = Form("%s %2.2f ", print_bar.Data(), percent);
    std::cout << text << "%\r" << std::flush; 
  }
} 

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
   
  // Prepare for loop over input MxAOD/TTree:
  std::vector<TString> fileNames = settings->getStrV("InputFile");
  TChain *chain = new TChain(settings->getStr("TreeName"));
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    chain->AddFile(fileNames[i_f]);
  }
  
  // USER MODIFICATION NECESSARY:
  // Assign the MxAOD/TTree branches to variables:
  int category = 0;
  double massVar; float weightVar; int resMass; 
  chain->SetBranchAddress(settings->getStr("MassBranchName"), &massVar);
  chain->SetBranchAddress(settings->getStr("WeightBranchName"), &weightVar);
  if (settings->isDefined("ResMassBranchName")) {
    chain->SetBranchAddress(settings->getStr("ResMassBranchName"), &resMass);
  }
  if (settings->isDefined("CateBranchName")) {
    chain->SetBranchAddress(settings->getStr("CateBranchName"), &category);
  }
  int nEvents = chain->GetEntries();
  
  int nCategories = 0;
  std::vector<double> mResList; mResList.clear();
  
  // Use a hard-coded resonance mass if a branch is not included in MxAODs:
  if (settings->isDefined("ResonanceMass") && 
      !settings->isDefined("ResMassBranchName")) {
    resMass = settings->getInt("ResonanceMass");
  }
  
  //--------------------------------------//
  // Loop over the Trees:
  // Loop over events to get parameterization:
  std::cout << "There are "<< nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {
    
    chain->GetEntry(index);
    PrintProgressBar(index, nEvents);
    
    // The category index fed into the SigParam tool should start at 0.
    int currCate = 0;// category;
    
    // The observed mass fed into the SigParam tool should always be in GeV.
    double massToUse = (double)massVar;
    if (TString(settings->getStr("MassBranchUnits")).EqualTo("MeV")) {
      massToUse = massVar / 1000.0;
    }
    
    if (massToUse < 0.0) continue;
    
    // The resonance mass fed into the SigParam tool should always be in GeV.
    // At the moment, the MxAODs do not have a branch for the truth Higgs mass.
    double resMassToUse = (double)resMass;
    if (TString(settings->getStr("ResMassBranchUnits")).EqualTo("MeV")) {
      resMassToUse = resMass / 1000.0;
    }
    
    // Counter for categories:
    if (currCate >= nCategories) nCategories = currCate+1;
    
    // USER MODIFICATION NECESSARY:
    double analysisLuminosity = 20300;
    double weightToUse = 1.0;//weightVar * analysisLuminosity;
    
    // Add the mass and weight values to the datasets for fitting:
    sps->addMassPoint(resMassToUse, currCate, massToUse, weightToUse);
    
    mResList = checkMResList(mResList, resMassToUse);
  }
  
  //--------------------------------------//
  // Now fit and plot the resonance shapes!
  std::cout << "createSignalParameterization: Start fitting and plotting!" << std::endl;
  
  // Loop over the analysis categories:
  for (int i_c = 0; i_c < nCategories; i_c++) {
    
    std::cout << "createSignalParameterization: # of masses: "
	      << mResList.size() << std::endl;
    for (int i_m = 0; i_m < (int)mResList.size(); i_m++) {
      if (sps->makeSingleResonance(mResList[i_m], i_c, function)) {
	sps->plotSingleResonance(mResList[i_m], i_c);
      }
      else {
	std::cout << "createSignalParameterization: Fit at mRes="
		  << mResList[i_m] << ", cate=" << i_c 
		  << " did not converge :(" << std::endl;
      }
    }
    sps->saveAll();
  }
  
  //--------------------------------------//
  // Print some values for resonance mass of 125 GeV and category 0:
  std::cout << "\nfitResonantDiHiggs: Printing parameters for cate 0, M_Res = "
	    << resMass << " GeV" << std::endl;
  std::vector<TString> names = sps->getVariableNames(resMass,0);
  for (int i_n = 0; i_n < (int)names.size(); i_n++) {
    std::cout << "\t" << names[i_n] << " = " 
	      << sps->getParameterValue(names[i_n],resMass,0) << std::endl;
  }
  return 0;
}
