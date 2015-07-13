////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMPseudoExp.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 23/06/2015                                                          //
//                                                                            //
//  This program creates pseudoexperiment ensembles for the H->yy + DM search //
//  with 13 TeV data at ATLAS.                                                //
//                                                                            //
//  Note: this program has been developed to work in conjunction with the     //
//  DMTestStat class, in order to reduce the code redundancy for calculation  //
//  of test statistics. Also, the input workspace is loaded in a loop because //
//  for some reason it is impossible to overwrite entries in a workspace.     //
//                                                                            //
//  options:                                                                  //
//      Binned, FixMu                                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMPseudoExp.h"

/**
   -----------------------------------------------------------------------------
   Generates a pseudo-experiment dataset. 
   @param w - input workspace containing the statistical model
   @param mc - the modelconfig object containing the statistical model
   @param seed - the random seed for data generation.
*/
void createPseudoData(RooWorkspace *w, ModelConfig *mc, int seed) {
  std::cout << "DMPseudoExp: Create pseudodata with seed " << seed << std::endl;
  
  RooSimultaneous* combPdf = (RooSimultaneous*)mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)mc->GetGlobalObservables();
  RooArgSet* observables = (RooArgSet*)mc->GetObservables();
  RooArgSet* originValsNP = (RooArgSet*)mc->GetNuisanceParameters()->snapshot();
  RooRealVar* firstPOI = (RooRealVar*)mc->GetParametersOfInterest()->first();
  
  RooRandom::randomGenerator()->SetSeed(seed);
  statistics::constSet(nuisanceParameters, true);
  statistics::constSet(globalObservables, false);
  
  map<string,RooDataSet*> toyDataMap; 
  RooCategory *categories = (RooCategory*)w->obj("categories");
  TIterator *cateIter = combPdf->indexCat().typeIterator();
  RooCatType *cateType = NULL;
  RooDataSet *dataTemp[20];
  
  // Loop over all channels
  int index = 0;
  // Previously this was commented and similar line below was uncommented
  statistics::randomizeSet(combPdf, globalObservables, seed); 
  statistics::constSet(globalObservables, true);
  
  numEventsPerCate.clear();
  
  // Iterate over the categories:
  while ((cateType=(RooCatType*) cateIter->Next())) {
    RooAbsPdf *currPDF = combPdf->getPdf(cateType->GetName());
    RooArgSet *currObs = currPDF->getObservables(observables);
    RooArgSet *currGlobs = currPDF->getObservables(globalObservables);
    RooRealVar *t = (RooRealVar*)currObs->first();
    
    //statistics::randomizeSet(currPDF, currGlobs, -1);
    //statistics::constSet(currGlobs, true);
    
    // Option for producing binned pseudo-data:
    if (options.Contains("Binned")) {
      currPDF->setAttribute("PleaseGenerateBinned");
      TIterator *iterObs = currObs->createIterator();
      RooRealVar *currObs = NULL;
      // Bin each of the observables:
      while ((currObs = (RooRealVar*)iterObs->Next())) currObs->setBins(120);
      dataTemp[index]
	= (RooDataSet*)currPDF->generate(*currObs, AutoBinned(true),
					 Extended(currPDF->canBeExtended()),
					 GenBinned("PleaseGenerateBinned"));
    }
    // Default: construct unbinned pseudo-data:
    else {
      dataTemp[index] = (RooDataSet*)currPDF->generate(*currObs,Extended(true));
    }
    
    toyDataMap[(std::string)cateType->GetName()] = dataTemp[index];
    numEventsPerCate.push_back((double)dataTemp[index]->sumEntries());
    index++;
  }
    
  // Import the new data into the workspace:
  RooDataSet* toyData = new RooDataSet("toyData", "", *observables, 
				       RooFit::Index(*categories),
				       RooFit::Import(toyDataMap));
  
  // release nuisance parameters:
  statistics::constSet(nuisanceParameters, false);
  
  // Import into the workspace:
  w->import(*toyData);
}

/**
   -----------------------------------------------------------------------------
   The main method. 
   @param jobName - the name of the analysis job.
   @param DMSignal - the name of the DM signal model.
   @param cateScheme - the category scheme for the analysis.
   @param options - the options (see header note).
   @param seed - the random seed for pseudoexperiment creation.
   @param toysPerJob - the number of pseudoexperiments to create per job.
   @param muDMVal - the value of the DM signal strength to use.
*/
int main(int argc, char **argv) {
  if (argc < 8) {
    std::cout << "Usage: " << argv[0] << " <jobName> <DMSignal> <cateScheme> <options> <seed> <toysPerJob> <mu_DM>" << std::endl;
    exit(0);
  }
  
  // Assign input parameters:
  TString jobName = argv[1];
  TString DMSignal = argv[2];
  TString cateScheme = argv[3];
  options = argv[4];
  int seed = atoi(argv[5]);
  int nToysPerJob = atoi(argv[6]);
  int inputMuDM = atoi(argv[7]);
  
  // Copy the input workspace file locally:
  TString originFile = Form("%s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root",
			    DMAnalysis::masterOutput.Data(), jobName.Data(), 
			    DMSignal.Data());
  TString copiedFile = Form("workspaceDM_%s.root",DMSignal.Data());
  system(Form("cp %s %s", originFile.Data(), copiedFile.Data()));
  
  // Create output TTree:
  TString outputDir = Form("%s/%s/DMPseudoExp", DMAnalysis::masterOutput.Data(),
			   jobName.Data());
  TString tempOutputFileName = Form("%s/single_files/toy_mu%i_%i.root",
				    outputDir.Data(), inputMuDM, seed);
  
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
  double muDMVal, nllMu0, nllMu1, nllMuFree, llrL1L0, llrL0Lfree, llrL1Lfree;
  
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
  fOutputTree.Branch("muDMVal", &muDMVal, "muDMVal/D");
  fOutputTree.Branch("convergedMu0", &convergedMu0, "convergedMu0/O");
  fOutputTree.Branch("convergedMu1", &convergedMu1, "convergedMu1/O");
  fOutputTree.Branch("convergedMuFree", &convergedMuFree, "convergedMuFree/O");
  fOutputTree.Branch("nllMu0", &nllMu0, "nllMu0/D");
  fOutputTree.Branch("nllMu1", &nllMu1, "nllMu1/D");
  fOutputTree.Branch("nllMuFree", &nllMuFree, "nllMuFree/D");
  fOutputTree.Branch("llrL1L0", &llrL1L0, "llrL1L0/D");
  fOutputTree.Branch("llrL0Lfree", &llrL0Lfree, "llrL0Lfree/D");
  fOutputTree.Branch("llrL1Lfree", &llrL1Lfree, "llrL1Lfree/D");
  fOutputTree.Branch("numEventsPerCate", &numEventsPerCate);
  fOutputTree.Branch("namesNP", &namesNP);
  fOutputTree.Branch("valuesNPMu0", &valuesNPMu0);
  fOutputTree.Branch("valuesNPMu1", &valuesNPMu1);
  fOutputTree.Branch("valuesNPMuFree", &valuesNPMuFree);
  fOutputTree.Branch("namesGlobs", &namesGlobs);
  fOutputTree.Branch("valuesGlobsMu1", &valuesGlobsMu1);
  fOutputTree.Branch("valuesGlobsMu0", &valuesGlobsMu0);
  fOutputTree.Branch("valuesGlobsMuFree", &valuesGlobsMuFree);
  
  // Loop to generate pseudo experiments:
  std::cout << "DMPseudoExp: Generating " << nToysPerJob
	    << " toys with mu_DM = " << inputMuDM << endl;
  for (int i_t = 0; i_t < nToysPerJob; i_t++) {
    
    // Load model, data, etc. from workspace:
    TFile inputFile(copiedFile, "read");
    RooWorkspace *workspace = (RooWorkspace*)inputFile.Get("combinedWS");    
    ModelConfig *mc = (ModelConfig*)workspace->obj("modelConfig");
    TString dataName = (DMAnalysis::doBlind) ? "asimovDataMu0" : "obsData";
    if (options.Contains("FixBkgParam")) dataName = "asimovDataMu0";
    RooAbsData *obsData = workspace->data(dataName);
    
    // Get the POI (mu):
    RooRealVar* firstPOI = (RooRealVar*)mc->GetParametersOfInterest()->first();
    //statistics::setDefaultPrintLevel(0);
    //statistics::setDefaultMinimizer("Minuit");
    
    // Load snapshot of profiled data:
    workspace->loadSnapshot("paramsOrigin");
    RooArgSet* nuisAndPOI = new RooArgSet();
    nuisAndPOI->add(*(RooArgSet*)mc->GetNuisanceParameters());
    nuisAndPOI->add(*(RooArgSet*)mc->GetParametersOfInterest());
    double fixedMu = (double)inputMuDM;
    
    // Load the newly-created snapshot:
    workspace->loadSnapshot("paramsOrigin");
    if (options.Contains("FixMu")) {
      firstPOI->setVal(fixedMu);
      firstPOI->setConstant(true);
    }
    
    // Create a pseudo-dataset:
    createPseudoData(workspace, mc, seed);
    numEvents = workspace->data("toyData")->sumEntries();
    
    firstPOI->setConstant(false);
    
    DMTestStat *dmts = new DMTestStat(jobName, DMSignal, cateScheme, "new",
				      workspace);
    
    // Mu = 0 fits:
    nllMu0 = dmts->getFitNLL("toyData", 0, true, muDMVal);
    convergedMu0 = dmts->fitsAllConverged();
    namesNP = dmts->getNPNames();
    valuesNPMu0 = dmts->getNPValues();
    namesGlobs = dmts->getGlobsNames();
    valuesGlobsMu0 = dmts->getGlobsValues();
    
    // Mu = 1 fits:
    nllMu1 = dmts->getFitNLL("toyData", 1, true, muDMVal);
    convergedMu1 = dmts->fitsAllConverged();
    valuesNPMu1 = dmts->getNPValues();
    valuesGlobsMu1 = dmts->getGlobsValues();
    
    // Mu free fits:
    nllMuFree = dmts->getFitNLL("toyData", 1, false, muDMVal);
    convergedMuFree = dmts->fitsAllConverged();
    valuesNPMuFree = dmts->getNPValues();
    valuesGlobsMuFree = dmts->getGlobsValues();
    
    // Calculate profile likelihood ratios:
    llrL1L0 = nllMu1 - nllMu0;
    llrL1Lfree = muDMVal > 1.0 ? 0.0 : (nllMu1 - nllMuFree);
    llrL0Lfree = muDMVal < 0.0 ? 0.0 : (nllMu0 - nllMuFree);
    
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
