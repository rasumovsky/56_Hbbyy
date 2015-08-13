////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHPseudoExp.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 23/06/2015                                                          //
//                                                                            //
//  This program creates pseudoexperiment ensembles for the di-Higgs bbyy     //
//  search with 13 TeV data at ATLAS.                                         //
//                                                                            //
//  Note: this program has been developed to work in conjunction with the     //
//  DHTestStat class, in order to reduce the code redundancy for calculation  //
//  of test statistics. Also, the input workspace is loaded in a loop because //
//  for some reason it is impossible to overwrite entries in a workspace.     //
//                                                                            //
//  options:                                                                  //
//      Binned, FixMu                                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Package includes:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "DHAnalysis.h"
#include "DHTestStat.h"
#include "RooBernsteinM.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "statistics.h"

/**
   -----------------------------------------------------------------------------
   The main method. 
   @param configFile - the name of the analysis config file.
   @param DHSignal - the name of the DH signal model.
   @param options - the options (see header note).
   @param seed - the random seed for pseudoexperiment creation.
   @param toysPerJob - the number of pseudoexperiments to create per job.
   @param muDHVal - the value of the DH signal strength to use.
*/
int main(int argc, char **argv) {
  if (argc < 7) {
    std::cout << "Usage: " << argv[0] << " <configFile> <DHSignal> <options> <seed> <toysPerJob> <mu_DH>" << std::endl;
    exit(0);
  }
  
  // Assign input parameters:
  TString configFile = argv[1];
  TString DHSignal = argv[2];
  TString options = argv[3];
  int seed = atoi(argv[4]);
  int nToysPerJob = atoi(argv[5]);
  int inputMuDH = atoi(argv[6]);
  
  // Load the analysis configurations from file:
  Config *config = new Config(configFile);
  TString anaType = DHAnalysis::getAnalysisType(config, DHSignal);
  // Copy the input workspace file locally:
  TString originFile = Form("%s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
			    (config->getStr("masterOutput")).Data(), 
			    (config->getStr("jobName")).Data(), anaType.Data());
  TString copiedFile = Form("workspaceDH_%s.root", anaType.Data());
  system(Form("cp %s %s", originFile.Data(), copiedFile.Data()));
  
  // Create output TTree:
  TString outputDir = Form("%s/%s/DHPseudoExp", 
			   (config->getStr("masterOutput")).Data(),
			   (config->getStr("jobName")).Data());
  TString tempOutputFileName = Form("%s/single_files/toy_mu%i_%i.root",
				    outputDir.Data(), inputMuDH, seed);
  
  // Construct the output directories:
  system(Form("mkdir -vp %s/err", outputDir.Data()));
  system(Form("mkdir -vp %s/log", outputDir.Data()));
  system(Form("mkdir -vp %s/single_files", outputDir.Data()));
  
  // Output file:
  TFile fOutputFile(tempOutputFileName, "recreate");
  TTree fOutputTree("toy", "toy");
  
  // Variables to store in the TTree:
  double numEvents;
  bool convergedMu1, convergedMu0, convergedMuFree;
  double muDHVal, nllMu0, nllMu1, nllMuFree, llrL1L0, llrL0Lfree, llrL1Lfree;
  
  std::vector<std::string> namesNP; namesNP.clear();
  std::vector<double> valuesNPMu0; valuesNPMu0.clear();
  std::vector<double> valuesNPMu1; valuesNPMu1.clear();
  std::vector<double> valuesNPMuFree; valuesNPMuFree.clear();
  std::vector<string> namesGlobs; namesGlobs.clear();
  std::vector<double> valuesGlobsMu0; valuesGlobsMu0.clear();
  std::vector<double> valuesGlobsMu1; valuesGlobsMu1.clear();
  std::vector<double> valuesGlobsMuFree; valuesGlobsMuFree.clear();
    
  fOutputTree.Branch("seed", &seed, "seed/I");
  fOutputTree.Branch("numEvents", &numEvents, "numEvents/D");
  fOutputTree.Branch("muDHVal", &muDHVal, "muDHVal/D");
  fOutputTree.Branch("convergedMu0", &convergedMu0, "convergedMu0/O");
  fOutputTree.Branch("convergedMu1", &convergedMu1, "convergedMu1/O");
  fOutputTree.Branch("convergedMuFree", &convergedMuFree, "convergedMuFree/O");
  fOutputTree.Branch("nllMu0", &nllMu0, "nllMu0/D");
  fOutputTree.Branch("nllMu1", &nllMu1, "nllMu1/D");
  fOutputTree.Branch("nllMuFree", &nllMuFree, "nllMuFree/D");
  fOutputTree.Branch("llrL1L0", &llrL1L0, "llrL1L0/D");
  fOutputTree.Branch("llrL0Lfree", &llrL0Lfree, "llrL0Lfree/D");
  fOutputTree.Branch("llrL1Lfree", &llrL1Lfree, "llrL1Lfree/D");
  //fOutputTree.Branch("numEventsPerCate", &numEventsPerCate);
  fOutputTree.Branch("namesNP", &namesNP);
  fOutputTree.Branch("valuesNPMu0", &valuesNPMu0);
  fOutputTree.Branch("valuesNPMu1", &valuesNPMu1);
  fOutputTree.Branch("valuesNPMuFree", &valuesNPMuFree);
  fOutputTree.Branch("namesGlobs", &namesGlobs);
  fOutputTree.Branch("valuesGlobsMu1", &valuesGlobsMu1);
  fOutputTree.Branch("valuesGlobsMu0", &valuesGlobsMu0);
  fOutputTree.Branch("valuesGlobsMuFree", &valuesGlobsMuFree);
  
  // Loop to generate pseudo experiments:
  std::cout << "DHPseudoExp: Generating " << nToysPerJob
	    << " toys with mu_DH = " << inputMuDH << endl;
  for (int i_t = 0; i_t < nToysPerJob; i_t++) {
    
    // Load model, data, etc. from workspace:
    TFile inputFile(copiedFile, "read");
    RooWorkspace *workspace = (RooWorkspace*)inputFile.Get("combinedWS");    
     
    DHTestStat *dhts = new DHTestStat(configFile, DHSignal, "new", workspace);
    
    RooDataSet *newToyData
      = dhts->createPseudoData(seed, inputMuDH, 1, options.Contains("FixMu"));
    numEvents = workspace->data("toyData")->sumEntries();
    
    // Mu = 0 fits:
    nllMu0 = dhts->getFitNLL("toyData", 0, true, muDHVal);
    convergedMu0 = dhts->fitsAllConverged();
    namesNP = dhts->getNPNames();
    valuesNPMu0 = dhts->getNPValues();
    namesGlobs = dhts->getGlobsNames();
    valuesGlobsMu0 = dhts->getGlobsValues();
    
    // Mu = 1 fits:
    nllMu1 = dhts->getFitNLL("toyData", 1, true, muDHVal);
    convergedMu1 = dhts->fitsAllConverged();
    valuesNPMu1 = dhts->getNPValues();
    valuesGlobsMu1 = dhts->getGlobsValues();
    
    // Mu free fits:
    nllMuFree = dhts->getFitNLL("toyData", 1, false, muDHVal);
    convergedMuFree = dhts->fitsAllConverged();
    valuesNPMuFree = dhts->getNPValues();
    valuesGlobsMuFree = dhts->getGlobsValues();
    
    // Calculate profile likelihood ratios:
    llrL1L0 = nllMu1 - nllMu0;
    llrL1Lfree = muDHVal > 1.0 ? 0.0 : (nllMu1 - nllMuFree);
    llrL0Lfree = muDHVal < 0.0 ? 0.0 : (nllMu0 - nllMuFree);
    
    // Fill the tree:
    fOutputTree.Fill();
    
    // Count the toys:
    seed++;
    fOutputTree.AutoSave("SaveSelf");
    
    // Close the input file before the loop repeats:
    inputFile.Close();
  }
  
  // Write the output file, delete local file copies:
  fOutputFile.cd();
  fOutputTree.Write();
  fOutputFile.Close();
  system(Form("rm %s",copiedFile.Data()));
  return 0;
}
