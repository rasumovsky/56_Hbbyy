////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHPseudoExp.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 01/21/2016                                                          //
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
//      Binned, FixMu, CL%d                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Package includes:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "DHTestStat.h"
#include "RooBernsteinM.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "statistics.h"

/**
   -----------------------------------------------------------------------------
   Convert key and value in a map to two separate vectors that are passed by 
   reference. Stupid, yeah...
   @param map - The input map with keys and values to be split.
   @names - A vector of names.
   @values - A vector of values.
*/
void mapToVectors(std::map<std::string,double> map, 
		  std::vector<std::string>& names, std::vector<double>& values){
  names.clear(); 
  values.clear();
  for (std::map<std::string,double>::iterator mapIter = map.begin(); 
       mapIter != map.end(); mapIter++) {
    names.push_back(mapIter->first);
    values.push_back(mapIter->second);
  }
}

/**
   -----------------------------------------------------------------------------
   The main method tosses toys and saves data in a TTree.
   @param configFile - The name of the analysis config file.
   @param options - The options (see header note).
   @param seed - The random seed for pseudoexperiment creation.
   @param toysPerJob - The number of pseudoexperiments to create per job.
   @param poiVal - The value of the parameter of interest to use.
   @param resonanceMass - The mass of the resonance (resonant analysis only)
*/
int main(int argc, char **argv) {
  if (argc < 7) {
    std::cout << "Usage: " << argv[0] 
	      << " <configFile> <options> <seed> <toysPerJob> <poiVal>" 
	      << std::endl;
    exit(0);
  }
  
  // Assign input parameters:
  TString configFile = argv[1];
  TString options = argv[2];
  int seed = atoi(argv[3]);
  int nToysPerJob = atoi(argv[4]);
  int inputPoIVal = atoi(argv[5]);
  int resonanceMass = atoi(argv[6]);
  
  // Load the analysis configurations from file:
  Config *config = new Config(configFile);
  TString anaType = config->getStr("AnalysisType");
  
  // Copy the input workspace file locally:
  TString originFile = Form("%s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
			    (config->getStr("MasterOutput")).Data(), 
			    (config->getStr("JobName")).Data(), anaType.Data());
  TString copiedFile = Form("workspaceDH_%s.root", anaType.Data());
  system(Form("cp %s %s", originFile.Data(), copiedFile.Data()));
  
  // Create output TTree:
  TString outputDir = Form("%s/%s/DHPseudoExp", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  
  if (options.Contains(config->getStr("CLScanToyOptions"))) {
    if (anaType.EqualTo("Resonant")) {
      outputDir = Form("%s/%s/DHPseudoExpForCLScan/MX%d", 
		       (config->getStr("MasterOutput")).Data(),
		       (config->getStr("JobName")).Data(), resonanceMass);
    }
    else{
      outputDir = Form("%s/%s/DHPseudoExpForCLScan", 
		       (config->getStr("MasterOutput")).Data(),
		       (config->getStr("JobName")).Data());
    }
  }
  
  // Get the index of the job (for CL Scan only):
  TString jobIndex = options;
  jobIndex.ReplaceAll(config->getStr("CLScanToyOptions"),"");
  
  TString tempOutputFileName = Form("%s/single_files/toy_mu%i_%i.root",
				    outputDir.Data(), inputPoIVal, seed);
  if (options.Contains(config->getStr("CLScanToyOptions"))) {
    tempOutputFileName = Form("%s/single_files/toy_mu%i_%i_%i.root",
			      outputDir.Data(), inputPoIVal, jobIndex.Atoi(),
			      seed);
  }
  
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
  double profiledPOIVal;
  double nllMu0, nllMu1, nllMuFree, llrL1L0, llrL0Lfree, llrL1Lfree;
  std::vector<std::string> namesNP; namesNP.clear();
  std::vector<double> valuesNPMu0; valuesNPMu0.clear();
  std::vector<double> valuesNPMu1; valuesNPMu1.clear();
  std::vector<double> valuesNPMuFree; valuesNPMuFree.clear();
  std::vector<string> namesGlobs; namesGlobs.clear();
  std::vector<double> valuesGlobsMu0; valuesGlobsMu0.clear();
  std::vector<double> valuesGlobsMu1; valuesGlobsMu1.clear();
  std::vector<double> valuesGlobsMuFree; valuesGlobsMuFree.clear();
  std::vector<string> namesPars; namesPars.clear();
  std::vector<double> valuesParsMu0; valuesParsMu0.clear();
  std::vector<double> valuesParsMu1; valuesParsMu1.clear();
  std::vector<double> valuesParsMuFree; valuesParsMuFree.clear();
  std::vector<double> numEventsPerCate; numEventsPerCate.clear();
  
  fOutputTree.Branch("seed", &seed, "seed/I");
  fOutputTree.Branch("numEvents", &numEvents, "numEvents/D");
  fOutputTree.Branch("numEventsPerCate", &numEventsPerCate);
  fOutputTree.Branch("profiledPOIVal", &profiledPOIVal, "profiledPOIVal/D");
  fOutputTree.Branch("convergedMu0", &convergedMu0, "convergedMu0/O");
  fOutputTree.Branch("convergedMu1", &convergedMu1, "convergedMu1/O");
  fOutputTree.Branch("convergedMuFree", &convergedMuFree, "convergedMuFree/O");
  fOutputTree.Branch("nllMu0", &nllMu0, "nllMu0/D");
  fOutputTree.Branch("nllMu1", &nllMu1, "nllMu1/D");
  fOutputTree.Branch("nllMuFree", &nllMuFree, "nllMuFree/D");
  fOutputTree.Branch("llrL1L0", &llrL1L0, "llrL1L0/D");
  fOutputTree.Branch("llrL0Lfree", &llrL0Lfree, "llrL0Lfree/D");
  fOutputTree.Branch("llrL1Lfree", &llrL1Lfree, "llrL1Lfree/D");
  fOutputTree.Branch("namesNP", &namesNP);
  fOutputTree.Branch("valuesNPMu0", &valuesNPMu0);
  fOutputTree.Branch("valuesNPMu1", &valuesNPMu1);
  fOutputTree.Branch("valuesNPMuFree", &valuesNPMuFree);
  fOutputTree.Branch("namesGlobs", &namesGlobs);
  fOutputTree.Branch("valuesGlobsMu1", &valuesGlobsMu1);
  fOutputTree.Branch("valuesGlobsMu0", &valuesGlobsMu0);
  fOutputTree.Branch("valuesGlobsMuFree", &valuesGlobsMuFree);
  fOutputTree.Branch("namesPars", &namesPars);
  fOutputTree.Branch("valuesParsMu1", &valuesParsMu1);
  fOutputTree.Branch("valuesParsMu0", &valuesParsMu0);
  fOutputTree.Branch("valuesParsMuFree", &valuesParsMuFree);
  
  // Loop to generate pseudo experiments:
  std::cout << "DHPseudoExp: Generating " << nToysPerJob
	    << " toys with mu_DH = " << inputPoIVal << endl;
  for (int i_t = 0; i_t < nToysPerJob; i_t++) {
    
    // Load model, data, etc. from workspace:
    TFile inputFile(copiedFile, "read");
    RooWorkspace *workspace = (RooWorkspace*)inputFile.Get("combinedWS");    
     
    DHTestStat *dhts = new DHTestStat(configFile, "new", workspace);
    
    if (options.Contains(config->getStr("CLScanToyOptions"))) {
      std::vector<double> scanValues = config->getNumV("CLScanValues");
      double crossSection = scanValues[jobIndex.Atoi()];
      dhts->setParam(config->getStr("CLScanVar"), crossSection, true);
    }
    
    // Be sure to set the value of the resonance mass for fitting!
    if ((config->getStr("AnalysisType")).EqualTo("Resonant")) {
      dhts->setParam("mResonance", resonanceMass, true);
    }

    // NEW! If background-only fit, set spurious signal constant:
    //if (config->getBool("UseSystematics") && inputPoIVal == 0) {
    //dhts->setParam("MyyMODELING", 0, true);
    //dhts->setParam("RNDM_MyyMODELING", 0, true);
    //}
    
    // Then create snapshot for Mu=1 or Mu=0 hypothesis! This must be re-done
    // for every toy job, since the signal hypothesis can change!
    dhts->saveSnapshots(true);
    TString dataToProf = (config->getBool("DoBlind")) ? 
      "asimovDataMu1" : "obsData";
    if (inputPoIVal == 0) dhts->getFitNLL(dataToProf, 0, true, profiledPOIVal);
    else dhts->getFitNLL(dataToProf, 1, true, profiledPOIVal);
    dhts->saveSnapshots(false);
    
    // Create the pseudo data:
    RooDataSet *newToyData
      = dhts->createPseudoData(seed, inputPoIVal, options.Contains("FixMu"));
    numEvents = workspace->data("toyData")->sumEntries();
    numEventsPerCate = dhts->getNEventsToys();

    // Globs are only set once for each dataset:
    mapToVectors(dhts->getGlobalObservables(), namesGlobs, valuesGlobsMu0);
    mapToVectors(dhts->getGlobalObservables(), namesGlobs, valuesGlobsMu1);
    mapToVectors(dhts->getGlobalObservables(), namesGlobs, valuesGlobsMuFree);
    
    // Mu = 0 fits:
    nllMu0 = dhts->getFitNLL("toyData", 0, true, profiledPOIVal, false);
    convergedMu0 = dhts->fitsAllConverged();
    mapToVectors(dhts->getNuisanceParameters(), namesNP, valuesNPMu0);
    mapToVectors(dhts->getParameters(), namesPars, valuesParsMu0);
    
    // Mu = 1 fits:
    nllMu1 = dhts->getFitNLL("toyData", 1, true, profiledPOIVal, false);
    convergedMu1 = dhts->fitsAllConverged();
    mapToVectors(dhts->getNuisanceParameters(), namesNP, valuesNPMu1);
    mapToVectors(dhts->getParameters(), namesPars, valuesParsMu1);
    
    // Mu free fits:
    nllMuFree = dhts->getFitNLL("toyData", 1, false, profiledPOIVal, false);
    convergedMuFree = dhts->fitsAllConverged();
    mapToVectors(dhts->getNuisanceParameters(), namesNP, valuesNPMuFree);
    mapToVectors(dhts->getParameters(), namesPars, valuesParsMuFree);
    
    // Calculate profile likelihood ratios:
    llrL1L0 = nllMu1 - nllMu0;
    llrL1Lfree = profiledPOIVal > 1.0 ? 0.0 : (nllMu1 - nllMuFree);
    llrL0Lfree = profiledPOIVal < 0.0 ? 0.0 : (nllMu0 - nllMuFree);
    
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
  system(Form("rm %s", copiedFile.Data()));
  return 0;
}
