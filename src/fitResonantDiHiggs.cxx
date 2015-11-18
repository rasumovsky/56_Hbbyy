////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  createSignalParameterization.cxx                                          //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 19/08/2015                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This main method provides a tool for performing individual fits to the    //
//  resonance Monte Carlo for the resonant bb+yy analysis. Settings for the   //
//  utility are provided in resDiHiggsConfig.cfg.                             //
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
  
  sps->doBinnedFit(settings->getBool("DoBinnedFit"), 1);
  
  // Prepare for loop over input MxAOD/TTree:
  std::vector<TString> fileNames = settings->getStrV("InputFile");
  TChain *chain = new TChain(settings->getStr("TreeName"));
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    chain->AddFile(fileNames[i_f]);
  }
  
  // USER MODIFICATION NECESSARY:
  // Assign the MxAOD/TTree branches to variables:
  int category; double massVar; float weightVar; int resMass; 
  chain->SetBranchAddress(settings->getStr("MassBranchName"), &massVar);
  if (settings->isDefined("WeightBranchName")) {
    chain->SetBranchAddress(settings->getStr("WeightBranchName"), &weightVar);
  }
  else {
    weightVar = 1.0;
  }
  if (settings->isDefined("ResMassBranchName")) {
    chain->SetBranchAddress(settings->getStr("ResMassBranchName"), &resMass);
  }
  else {
    resMass = settings->getNum("ResonanceMass");
  }
  if (settings->isDefined("CateBranchName")) {
    chain->SetBranchAddress(settings->getStr("CateBranchName"), &category);
  }
  else {
    category = 0;
  }
  int nEvents = chain->GetEntries();
  int nCategories = 0;
      
  //--------------------------------------//
  // Loop over the Trees:
  std::cout << "There are "<< nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {
    chain->GetEntry(index);
    PrintProgressBar(index, nEvents);
    
    // The observed mass fed into the SigParam tool should always be in GeV.
    double massToUse = (double)massVar;
    if (TString(settings->getStr("MassBranchUnits")).EqualTo("MeV")) {
      massToUse = massVar / 1000.0;
    }
    if (massToUse < 0.0) continue;
    
    // The resonance mass fed into the SigParam tool should always be in GeV.
    double resMassToUse = (double)resMass;
    if (settings->isDefined("ResMassBranchName") && 
	TString(settings->getStr("ResMassBranchUnits")).EqualTo("MeV")) {
      resMassToUse = resMass / 1000.0;
    }
    
    // The category index fed into the SigParam tool should start at 0.
    int currCate = category;
    // Counter for categories:
    if (currCate >= nCategories) nCategories = currCate+1;
    
    double weightToUse = weightVar;
    
    // Add the mass and weight values to the datasets for fitting:
    sps->addMassPoint(resMassToUse, currCate, massToUse, weightToUse);
  }
  
  //--------------------------------------//
  // Now fit and plot the resonance shapes!
  std::cout << "createSignalParameterization: Start fitting and plotting!"
	    << std::endl;
  // Loop over the analysis categories:
  for (int i_c = 0; i_c < nCategories; i_c++) {
    if (sps->makeSingleResonance(resMass, i_c, function)) {
      sps->plotSingleResonance(resMass, i_c);
    }
    else {
      std::cout << "createSignalParameterization: Fit at cate=" << i_c 
		<< " did not converge :(" << std::endl;
    }
    sps->saveAll();
  }
  
  //--------------------------------------//
  // Print some values for resonance mass of 125 GeV and category 0:
  std::cout << "\nfitResonantDiHiggs: Printing parameters for cate 0, M_Res = "
	    << resMass << " GeV" << std::endl;
  std::vector<TString> names = sps->getVariableNames(resMass, 0);
  for (int i_n = 0; i_n < (int)names.size(); i_n++) {
    std::cout << "\t" << names[i_n] << " = " 
	      << sps->getParameterValue(names[i_n], resMass, 0) << std::endl;
  }
  return 0;
}
