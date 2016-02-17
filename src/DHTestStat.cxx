////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHTestStat.cxx                                                      //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/11/2015                                                          //
//                                                                            //
//  This class allows the user to calculate p0, CL, and CLs based on an input //
//  workspace.                                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHTestStat.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the DHTestStat class. 
   @param newConfigFile - The analysis config file.
   @param newOptions - The job options ("New", "FromFile"), etc.
   @param newWorkspace - The workspace with the model for the test stats. 
*/
DHTestStat::DHTestStat(TString newConfigFile, TString newOptions,
		       RooWorkspace *newWorkspace) {
  std::cout << "DHTestStat: Initializing...\n\t" << newConfigFile << "\n\t"
	    << newOptions << "\n\t" << std::endl;
  
  // Assign input variables:
  m_options = newOptions;
  
  // Start with a clean class:
  clearData();
  
  // Set the analysis configuration:
  m_config = new Config(newConfigFile);
  m_jobName = m_config->getStr("JobName");
  m_anaType = m_config->getStr("AnalysisType");
  
  // Use Asimov data if the analysis is blind.
  m_dataForExpQ0 = "asimovDataMu1";
  m_dataForExpQMu = "asimovDataMu0";
  m_dataForObsQ0 = (m_config->getBool("DoBlind")) ? "asimovDataMu1" : "obsData";
  m_dataForObsQMu= (m_config->getBool("DoBlind")) ? "asimovDataMu0" : "obsData";
  
  // Load from file if the pointer passed is NULL:
  if (newWorkspace == NULL) {
    inputFile
      = new TFile(Form("%s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
		       (m_config->getStr("MasterOutput")).Data(),
		       m_jobName.Data(), m_anaType.Data()), "read");
    if (inputFile->IsOpen()) {
      printer("DHTestStat: Loading workspace.", false);
      m_workspace = (RooWorkspace*)inputFile->Get("combinedWS");
    }
    else {
      printer("DHTestStat: Error loading file, accessing with WS tool.", false);
      // Load the workspace from the nominal location.
      DHWorkspace *m_dhws = new DHWorkspace(newConfigFile, "FromFile");
      m_workspace = m_dhws->getCombinedWorkspace();
    }
  }
  // Use the workspace passed to the class constructor:
  else m_workspace = newWorkspace;
  
  // Then get the model configuration:
  m_mc = (ModelConfig*)m_workspace->obj("modelConfig");
  
  // Map storing all calculations:
  m_calculatedValues.clear();
  
  // Use linear plot y-axis by default:
  setPlotAxis(false, 0.0, 100.0, 1.0);

  // Create output directories:
  m_outputDir = Form("%s/%s/DHTestStat", 
		     (m_config->getStr("MasterOutput")).Data(),
		     m_jobName.Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  system(Form("mkdir -vp %s/CL/", m_outputDir.Data()));
  system(Form("mkdir -vp %s/p0/", m_outputDir.Data()));
  
  // Make new or load old values:
  if (m_options.Contains("FromFile")) loadStatsFromFile();
  
  // Finished instantiating the test statistic class. 
  std::cout << "DHTestStat: Initialized Successfully!" << std::endl;
  return;
}

/**
   -----------------------------------------------------------------------------
   Get the value of one of the test statistics.
   @param testStat - The test stat. name (p0, CL, CLs).
   @param observed - True iff observed, false if expected. 
   @param N - The standard deviation (-2, -1, 0, +1, +2). 
   @return - The test statistic value.
*/
double DHTestStat::accessValue(TString testStat, bool observed, int N) {
  TString currMapKey = getKey(testStat, observed, N);
  // Check that corresponding entry exists:
  if (mapValueExists(currMapKey)) return m_calculatedValues[currMapKey];
  else return 0;
}

/**
   -----------------------------------------------------------------------------
   Calculate the CL and CLs values using model fits.
*/
void DHTestStat::calculateNewCL() {
  printer("DHTestStat::calculateNewCL()", false);
  
  // Calculate observed qmu: 
  double muHatObs = 0.0;
  double nllMu1Obs = getFitNLL(m_dataForObsQMu, 1.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObsQMu, 1.0, false, muHatObs);
  double obsQMu = getQMuFromNLL(nllMu1Obs, nllMuHatObs, muHatObs, 1);
  
  // Calculate expected qmu:
  double muHatExp = 0.0;
  double nllMu1Exp = getFitNLL(m_dataForExpQMu, 1.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExpQMu, 0.0, false, muHatExp);
  double expQMu = getQMuFromNLL(nllMu1Exp, nllMuHatExp, muHatExp, 1);
  
  // Calculate CL:
  double expCLn2 = getCLFromQMu(expQMu, -2);
  double expCLn1 = getCLFromQMu(expQMu, -1);
  double expCLp1 = getCLFromQMu(expQMu, 1);
  double expCLp2 = getCLFromQMu(expQMu, 2);
  double expCL = getCLFromQMu(expQMu, 0);
  double obsCL = getCLFromQMu(obsQMu, 0);
  
  // Write CL values to file:
  ofstream textCL;
  textCL.open(Form("%s/CL/CL_values_%s.txt",
		   m_outputDir.Data(), m_anaType.Data()));
  textCL << m_anaType << " " << obsCL << " " << expCLn2 << " " << expCLn1
	 << " " << expCL << " " << expCLp1 << " " << expCLp2 << std::endl;
  textCL.close();
  
  // Print summary:
  std::cout << "  expected CL +2s = " << expCLp2 << std::endl;
  std::cout << "  expected CL +1s = " << expCLp1 << std::endl;
  std::cout << "  expected CL nom = " << expCL    << std::endl;
  std::cout << "  expected CL -1s = " << expCLn1 << std::endl;
  std::cout << "  expected CL -2s = " << expCLn2 << std::endl;
  std::cout << "  observed CL = " << obsCL << std::endl;
  std::cout << " " << std::endl;
  if (m_allGoodFits) std::cout << "All good fits? True" << std::endl;
  else std::cout << "All good fits? False" << std::endl;
  cout << " " << endl;
  if (obsQMu < 0) std::cout << "WARNING! obsQMu < 0 : " << obsQMu << std::endl;
  if (expQMu < 0) std::cout << "WARNING! expQMu < 0 : " << expQMu << std::endl;
  
  // Save CL and CLs for later access:
  m_calculatedValues[getKey("CL",0,-2)] = expCLn2;
  m_calculatedValues[getKey("CL",0,-1)] = expCLn1;
  m_calculatedValues[getKey("CL",0,0)] = expCL;
  m_calculatedValues[getKey("CL",0,1)] = expCLp1;
  m_calculatedValues[getKey("CL",0,2)] = expCLp2;
  m_calculatedValues[getKey("CL",1,0)] = obsCL;
  
  m_calculatedValues[getKey("CLs",0,-2)] = getCLsFromCL(expCLn2);
  m_calculatedValues[getKey("CLs",0,-1)] = getCLsFromCL(expCLn1);
  m_calculatedValues[getKey("CLs",0,0)] = getCLsFromCL(expCL);
  m_calculatedValues[getKey("CLs",0,1)] = getCLsFromCL(expCLp1);
  m_calculatedValues[getKey("CLs",0,2)] = getCLsFromCL(expCLp2);
  m_calculatedValues[getKey("CLs",1,0)] = getCLsFromCL(obsCL);
}

/**
   -----------------------------------------------------------------------------
   Calculate the p0 value using model fits.
*/
void DHTestStat::calculateNewP0() {
  printer("DHTestStat::calculateNewP0()", false);
  
  // Calculate observed q0: 
  double muHatObs = 0.0;
  double nllMu0Obs = getFitNLL(m_dataForObsQ0, 0.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObsQ0, 0.0, false, muHatObs);
  double obsQ0 = getQ0FromNLL(nllMu0Obs, nllMuHatObs, muHatObs);
  
  // Calculate expected q0:
  double muHatExp = 0.0;
  double nllMu0Exp = getFitNLL(m_dataForExpQ0, 0.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExpQ0, 0.0, false, muHatExp);
  double expQ0 = getQ0FromNLL(nllMu0Exp, nllMuHatExp, muHatExp);
  
  // Calculate p0 from q0:
  double expP0 = getP0FromQ0(expQ0);
  double obsP0 = getP0FromQ0(obsQ0);
  
  // Write p0 values to file:
  ofstream textP0;
  textP0.open(Form("%s/p0/p0_values_%s.txt", 
		   m_outputDir.Data(), m_anaType.Data()));
  textP0 << m_anaType << " " << expP0 << " " << obsP0 << std::endl;
  textP0.close();
  
  // Print summary:
  std::cout << "\n  Expected p0 = " << expP0 << std::endl;
  std::cout << "\n  Observed p0 = " << obsP0 << std::endl;
  if (fitsAllConverged()) std::cout << "All good fits? True\n" << std::endl;
  else std::cout << "All good fits? False\n" << std::endl;
  
  // Save p0 for later access:
  m_calculatedValues[getKey("p0", 1, 0)] = obsP0;
  m_calculatedValues[getKey("p0", 0, 0)] = expP0;
}

/**
   -----------------------------------------------------------------------------
   Clears all data stored by the class, but does not modify the workspace.
*/
void DHTestStat::clearData() {
  m_allGoodFits = true;
  m_calculatedValues.clear();
  m_mapGlobs.clear();
  m_mapNP.clear();
  m_mapPars.clear();
  m_doSaveSnapshot = false;
  m_doPlot = false;
  m_plotDir = "";
  clearFitParamSettings();
}

/**
   -----------------------------------------------------------------------------
   Clears all specifications for parameter values during fits.
*/
void DHTestStat::clearFitParamSettings() {
  m_paramValToSet.clear();
  m_paramConstToSet.clear();
}

/**
   -----------------------------------------------------------------------------
   Create a pseudo-dataset with a given value of DM and SM signal strength.
   @param seed - The random seed for dataset generation.
   @param valPoI - The value of the parameter of interest.
   @param fixPoI - True iff. the parameter of interest should be fixed. 
   @return - A pseudo-dataset.
*/
RooDataSet* DHTestStat::createPseudoData(int seed, int valPoI, bool fixPoI) {
  printer(Form("DHTestStatLLcreatePseudoData(seed=%d, PoI=%d, fixPoI=%d)",
	       seed, valPoI, (int)fixPoI), false);
  
  // Load the parameters from profiling data to create UNCONDITIONAL ENSEMBLE:
  TString snapshotName = Form("paramsProfilePoI%d", valPoI);
  if (m_workspace->getSnapshot(snapshotName)) {
    m_workspace->loadSnapshot(snapshotName);
    printer(Form("DHTestStat: Loaded snapshot %s",snapshotName.Data()), false);
  }
  else {
    // Load the original parameters from profiling:
    //m_workspace->loadSnapshot("paramsOrigin");
    printer(Form("DHTestStat: ERROR! No snapshot %s", snapshotName.Data()),
	    true);
  }
  
  //RooSimultaneous* combPdf = (RooSimultaneous*)m_mc->GetPdf();
  RooAbsPdf* combPdf = (RooAbsPdf*)m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  RooArgSet* observables = (RooArgSet*)m_mc->GetObservables();
  RooRealVar* firstPoI = (RooRealVar*)m_mc->GetParametersOfInterest()->first();
  
  //std::cout << "toy values BEFORE randomization etc." << std::endl;
  //printSet("nuisanceParameters", nuisanceParameters);
  //printSet("globalObservables", globalObservables);
  
  RooRandom::randomGenerator()->SetSeed(seed);
  statistics::constSet(nuisanceParameters, true);
  statistics::constSet(globalObservables, false);
      
  // Randomize the global observables and set them constant for now:
  statistics::randomizeSet(combPdf, globalObservables, seed); 
  statistics::constSet(globalObservables, true);
    
  // Set the parameter of interest value and status:
  firstPoI->setVal(valPoI);
  firstPoI->setConstant(fixPoI);
    
  // Check if other parameter settings have been specified for toys:
  // WARNING! This overrides the randomization settings above!
  for (std::map<TString,double>::iterator iterParam = m_paramValToSet.begin();
       iterParam != m_paramValToSet.end(); iterParam++) {
    // Check that workspace contains parameter:
    if (m_workspace->var(iterParam->first)) {
      m_workspace->var(iterParam->first)->setVal(iterParam->second);
      m_workspace->var(iterParam->first)
	->setConstant(m_paramConstToSet[iterParam->first]);
    }
    else {
      m_workspace->Print("v");
      printer(Form("DHTestStat: Error! Parameter %s not found in workspace! See printout above for clues...", (iterParam->first).Data()), true);
    }
  }
  
  //std::cout << "toy values AFTER randomization AND setting" << std::endl;
  //printSet("nuisanceParameters", nuisanceParameters);
  //printSet("globalObservables", globalObservables);
  
  // Store the toy dataset and number of events per dataset:
  map<string,RooDataSet*> toyDataMap; toyDataMap.clear();
  m_numEventsPerCate.clear();
  
  // Iterate over the categories:
  RooSimultaneous *simPdf = (RooSimultaneous*)m_workspace->pdf("combinedPdfSB");
  TIterator *cateIter = simPdf->indexCat().typeIterator();
  RooCatType *cateType = NULL;
  while ((cateType = (RooCatType*)cateIter->Next())) {
    //RooAbsPdf *currPdf = combPdf->getPdf(cateType->GetName());
    RooAbsPdf *currPdf = simPdf->getPdf(cateType->GetName());
    RooArgSet *currObs = currPdf->getObservables(observables);
        
    // If you want to bin the pseudo-data (speeds up calculation):
    if (m_options.Contains("Binned")) {
      currPdf->setAttribute("PleaseGenerateBinned");
      TIterator *iterObs = currObs->createIterator();
      RooRealVar *currObs = NULL;
      // Bin each of the observables:
      while ((currObs = (RooRealVar*)iterObs->Next())) currObs->setBins(120);
      toyDataMap[(std::string)cateType->GetName()]
	= (RooDataSet*)currPdf->generate(*currObs, AutoBinned(true),
					 Extended(currPdf->canBeExtended()),
					 GenBinned("PleaseGenerateBinned"));
    }
    // Construct unbinned pseudo-data by default:
    else {
      toyDataMap[(std::string)cateType->GetName()]
	= (RooDataSet*)currPdf->generate(*currObs,Extended(true));
    }
    double currEvt = toyDataMap[(std::string)cateType->GetName()]->sumEntries();
    m_numEventsPerCate.push_back(currEvt);
  }
  
  // Create the combined toy RooDataSet:
  RooCategory *categories = (RooCategory*)m_workspace->obj("categories");
  RooDataSet* pseudoData = new RooDataSet("toyData", "toyData", *observables, 
					  RooFit::Index(*categories),
					  RooFit::Import(toyDataMap));
  
  // Save the parameters used to generate toys:
  storeParams(nuisanceParameters, m_mapNP);
  storeParams(globalObservables, m_mapGlobs);
  storeParams((RooArgSet*)m_workspace->set("nonSysParameters"), m_mapPars);

  // release nuisance parameters (but don't change the values!):
  //m_workspace->loadSnapshot("paramsOrigin");
  statistics::constSet(nuisanceParameters, false);
  statistics::constSet(globalObservables, true);
  
  //std::cout << "toy values FINAL reloading." << std::endl;
  //printSet("nuisanceParameters", nuisanceParameters);
  //printSet("globalObservables", globalObservables);

  // Import into the workspace then return:
  m_workspace->import(*pseudoData);
  return pseudoData;
}

/*
   -----------------------------------------------------------------------------
   Create a pseudo-dataset with a given value of DM and SM signal strength.
   @param seed - The random seed for dataset generation.
   @param valPoI - The value of the parameter of interest.
   @param fixPoI - True iff. the parameter of interest should be fixed. 
   @return - A pseudo-dataset.

WARNING! THIS METHOD IS FLAWED. IT WILL BIAS YOUR TOYS AND DOES NOT CREATE
AN UNCONDITIONAL ENSEMBLE. IT WILL BIAS YOUR OBSERVATIONS...

RooDataSet* DHTestStat::createPseudoData(int seed, int valPoI, bool fixPoI) {
  printer(Form("DHTestStatLLcreatePseudoData(seed=%d, PoI=%d, fixPoI=%d)",
	       seed, valPoI, (int)fixPoI), false);
  
  //RooSimultaneous* combPdf = (RooSimultaneous*)m_mc->GetPdf();
  RooAbsPdf* combPdf = (RooAbsPdf*)m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  RooArgSet* observables = (RooArgSet*)m_mc->GetObservables();
  RooRealVar* firstPoI = (RooRealVar*)m_mc->GetParametersOfInterest()->first();
  
  // Load snapshot with all nuisance parameters and global observables at 0:
  m_workspace->loadSnapshot("paramsOrigin");

  //std::cout << "toy values FROM ORIGIN etc." << std::endl;
  //printSet("nuisanceParameters", nuisanceParameters);
  //printSet("globalObservables", globalObservables);
  
  // Randomize the global observables and set them constant for now:
  RooRandom::randomGenerator()->SetSeed(seed);
  statistics::constSet(globalObservables, false);
  statistics::randomizeSet(combPdf, globalObservables, seed); 
  statistics::constSet(globalObservables, true);
  
  //std::cout << "toy values AFTER randomization of globs" << std::endl;
  //printSet("nuisanceParameters", nuisanceParameters);
  //printSet("globalObservables", globalObservables);
  
  // Load the parameters from profiling data to create UNCONDITIONAL ENSEMBLE:
  RooArgSet *profiledSnapshot = NULL;
  TString snapshotName = Form("paramsProfilePoI%d", valPoI);
  if (m_workspace->getSnapshot(snapshotName)) {
    m_workspace->loadSnapshot(snapshotName);
    printer(Form("DHTestStat: Loaded snapshot %s",snapshotName.Data()), false);
  }
  else {
    // Load the original parameters from profiling:
    //m_workspace->loadSnapshot("paramsOrigin");
    printer(Form("DHTestStat: ERROR! No snapshot %s", snapshotName.Data()),
	    true);
  }
  
  statistics::constSet(nuisanceParameters, true, profiledSnapshot);
  
  //std::cout << "toy values AFTER randomization AND profiled NP" << std::endl;
  //printSet("nuisanceParameters", nuisanceParameters);
  //printSet("globalObservables", globalObservables);
  
  // Store the toy dataset and number of events per dataset:
  map<string,RooDataSet*> toyDataMap; toyDataMap.clear();
  m_numEventsPerCate.clear();
  
  // Set the parameter of interest value and status:
  firstPoI->setVal(valPoI);
  firstPoI->setConstant(fixPoI);
    
  // Check if other parameter settings have been specified for toys:
  // WARNING! This overrides the randomization settings above!
  for (std::map<TString,double>::iterator iterParam = m_paramValToSet.begin();
       iterParam != m_paramValToSet.end(); iterParam++) {
    // Check that workspace contains parameter:
    if (m_workspace->var(iterParam->first)) {
      m_workspace->var(iterParam->first)->setVal(iterParam->second);
      m_workspace->var(iterParam->first)
	->setConstant(m_paramConstToSet[iterParam->first]);
    }
    else {
      m_workspace->Print("v");
      printer(Form("DHTestStat: Error! Parameter %s not found in workspace! See printout above for clues...", (iterParam->first).Data()), true);
    }
  }
  
  std::cout << "toy values AFTER randomization AND setting" << std::endl;
  printSet("nuisanceParameters", nuisanceParameters);
  printSet("globalObservables", globalObservables);

  // Iterate over the categories:
  RooSimultaneous *simPdf = (RooSimultaneous*)m_workspace->pdf("combinedPdfSB");
  TIterator *cateIter = simPdf->indexCat().typeIterator();
  RooCatType *cateType = NULL;
  while ((cateType = (RooCatType*)cateIter->Next())) {
    //RooAbsPdf *currPdf = combPdf->getPdf(cateType->GetName());
    RooAbsPdf *currPdf = simPdf->getPdf(cateType->GetName());
    RooArgSet *currObs = currPdf->getObservables(observables);
        
    // If you want to bin the pseudo-data (speeds up calculation):
    if (m_options.Contains("Binned")) {
      currPdf->setAttribute("PleaseGenerateBinned");
      TIterator *iterObs = currObs->createIterator();
      RooRealVar *currObs = NULL;
      // Bin each of the observables:
      while ((currObs = (RooRealVar*)iterObs->Next())) currObs->setBins(120);
      toyDataMap[(std::string)cateType->GetName()]
	= (RooDataSet*)currPdf->generate(*currObs, AutoBinned(true),
					 Extended(currPdf->canBeExtended()),
					 GenBinned("PleaseGenerateBinned"));
    }
    // Construct unbinned pseudo-data by default:
    else {
      toyDataMap[(std::string)cateType->GetName()]
	= (RooDataSet*)currPdf->generate(*currObs,Extended(true));
    }
    double currEvt = toyDataMap[(std::string)cateType->GetName()]->sumEntries();
    m_numEventsPerCate.push_back(currEvt);
  }
  
  // Create the combined toy RooDataSet:
  RooCategory *categories = (RooCategory*)m_workspace->obj("categories");
  RooDataSet* pseudoData = new RooDataSet("toyData", "toyData", *observables, 
					  RooFit::Index(*categories),
					  RooFit::Import(toyDataMap));
  
  // Save the parameters used to generate toys:
  storeParams(nuisanceParameters, m_mapNP);
  storeParams(globalObservables, m_mapGlobs);
  storeParams((RooArgSet*)m_workspace->set("nonSysParameters"), m_mapPars);
  
  // release nuisance parameters (but don't change the values!):
  //m_workspace->loadSnapshot("paramsOrigin");
  statistics::constSet(nuisanceParameters, false);
  statistics::constSet(globalObservables, true);
  
  //std::cout << "toy values FINAL reloading." << std::endl;
  //printSet("nuisanceParameters", nuisanceParameters);
  //printSet("globalObservables", globalObservables);

  // Import into the workspace then return:
  m_workspace->import(*pseudoData);
  return pseudoData;
}
*/

/**
   -----------------------------------------------------------------------------
   Check if all of the fits done by this class have converged.
   @return - True iff. all of the fits have been successfully convergent.
*/
bool DHTestStat::fitsAllConverged() { 
  return m_allGoodFits;
}

/**
   -----------------------------------------------------------------------------
   Implements the functional form of qMu.
   @param x - The value of the test statistic.
   @return - The value of the asymptotic test statistic distribution.
*/
double DHTestStat::functionQMu(double x) {
  // This corresponds to the "special case" of mu=mu'
  double result = TMath::Exp(-1*x/2.0) / (2.0*sqrt(2.0*TMath::Pi()*x));
  return result;
}

/**
   -----------------------------------------------------------------------------
   Implements the functional form of qMuTilde.
   @param x - The value of the test statistic.
   @param asimovTestStat - The test stat value on Asimov data with mu=0 but
   fitting under mu=1 hypothesis.
   @return - The value of the asymptotic test statistic distribution.
*/
double DHTestStat::functionQMuTilde(double x, double asimovTestStat) {
  // This corresponds to the "special case" of mu=mu'
  double result = 0.0;
  double cutoff = asimovTestStat; // asimov test stat...
  if (x == 0) result = 0.0;
  else if (x <= cutoff) {
    result = (TMath::Exp(-1*x/2.0) / (2.0*sqrt(2.0*TMath::Pi()*x)));
  }
  else {
    result = (TMath::Exp(-1*TMath::Power((x+cutoff),2) / (8*cutoff))
	      / (2*sqrt(2*TMath::Pi()*cutoff)));
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Get the CL value from CLs.
   @param CLs - The CLs value to convert to CL.
   @return - The corresponding CL value.
*/
double DHTestStat::getCLFromCLs(double CLs) {
  return (1.0 - CLs);
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value from CL.
   @param CL - The CL value to convert to CLs.
   @return - The corresponding CLs value.
*/
double DHTestStat::getCLsFromCL(double CL) {
  return (1.0 - CL);
}

/**
   -----------------------------------------------------------------------------
   Get the CL value using qMu and the type.
   @param qMu - The value of the test statistic.
   @param N - The sigma value (-2,-1,0,1,2). Use 0 for median.
   @return - The CLs value.
*/
double DHTestStat::getCLFromQMu(double qMu, double N) {
  double CL = getCLFromCLs(getCLsFromQMu(qMu, N));
  return CL;
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value using qMu and the type.
   @param qMu - The value of the test statistic.
   @param N - The sigma value (-2,-1,0,1,2). Use 0 for median.
   @return - The CLs value.
*/
double DHTestStat::getCLsFromQMu(double qMu, double N) {
  // N = 0 for exp and obs
  double pMu = getPMuFromQMu(qMu);
  double pB = getPbFromN(N);
  double CLs = pMu / (1.0 - pB);
  return CLs;
}

/**
   -----------------------------------------------------------------------------
   Get the negative-log-likelihood for a fit of a specified type to a specified
   dataset.
   @param datasetName - The name of the dataset in the workspace.
   @param valPoI - The parameter of interest value to fix.
   @param fixPoI - True if PoI should be fixed to the specified value.
   @param &profiledValPoI - The profiled value of mu (passed by reference)
   @param resetParams - True iff original parameter values are used at start. 
   @return - The NLL value.
*/
double DHTestStat::getFitNLL(TString datasetName, double valPoI, bool fixPoI,
			     double &profiledValPoI, bool resetParams) {
  printer(Form("DHTestStat: getFitNLL(%s, PoI=%f, fixPoI=%d, resetPars=%d)",
	       datasetName.Data(),valPoI,(int)fixPoI,(int)resetParams),false);
  
  RooAbsPdf* combPdf = m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  if (resetParams) m_workspace->loadSnapshot("paramsOrigin");
  RooArgSet* origValNP = (RooArgSet*)m_workspace->getSnapshot("paramsOrigin");
  RooArgSet* poi = (RooArgSet*)m_mc->GetParametersOfInterest();
  RooRealVar* firstPoI = (RooRealVar*)poi->first();
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*nuisanceParameters);
  poiAndNuis->add(*poi);
  
  // Look for dataset. Create if non-existent & Asimov or requires binning.
  if (!m_workspace->data(datasetName)) {
    printer(Form("DHTestStat: Error! Requested data not available: %s",
		 datasetName.Data()), true);
  }
  
  // release nuisance parameters before fit and set to the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  // the global observables should be fixed to the nominal values...
  statistics::constSet(globalObservables, true);
  
  //std::cout << "Pre-fit parameter values" << std::endl;
  //printSet("nuisanceParameters", nuisanceParameters);
  //printSet("globalObservables", globalObservables);
  
  firstPoI->setVal(valPoI);
  firstPoI->setConstant(fixPoI);
  
  // NEW! If background-only fit, set spurious signal constant:
  //if (m_config->getBool("UseSystematics") && fixPoI && ((int)valPoI) == 0) {
  //m_workspace->var("MyyMODELING")->setVal(0);
  //m_workspace->var("MyyMODELING")->setConstant(true);
  //}
  
  // Check if other parameter settings have been specified for fit:
  for (std::map<TString,double>::iterator iterParam = m_paramValToSet.begin();
       iterParam != m_paramValToSet.end(); iterParam++) {
    // Check that workspace contains parameter:
    if (m_workspace->var(iterParam->first)) {
      m_workspace->var(iterParam->first)->setVal(iterParam->second);
      m_workspace->var(iterParam->first)
	->setConstant(m_paramConstToSet[iterParam->first]);
    }
    else {
      m_workspace->Print("v");
      printer(Form("DHTestStat: Error! Parameter %s not found in workspace! See printout above for clues...", ((TString)iterParam->first).Data()), true);
    }
  }
  
  // The actual fit command:
  RooNLLVar* varNLL
    = (RooNLLVar*)combPdf->createNLL(*m_workspace->data(datasetName),
				     Extended(combPdf->canBeExtended()));
  
  RooFitResult *fitResult = statistics::minimize(varNLL, "", NULL, true);
  if (!fitResult || fitResult->status() != 0) m_allGoodFits = false;
  
  /*
  std::cout << "PRINTING NORMALIZATIONS AND ERRORS" << std::endl;
  double bkgVal_jj = m_workspace->var("n_BkgNonHiggs_jj")->getVal();
  double bkgErr_jj =  m_workspace->var("n_BkgNonHiggs_jj")->getError();
  double SigSMVal_jj = m_workspace->function("n_SigSM_jj")->getVal();
  double SigSMErr_jj
    = m_workspace->function("n_SigSM_jj")->getPropagatedError(*fitResult);
  
  std::cout << "bkgVal_jj = " << bkgVal_jj << " +/- " << bkgErr_jj << std::endl;
  std::cout << "SigSMVal_jj = " << SigSMVal_jj << " +/- " << SigSMErr_jj
	    << std::endl;
  
  
  double bkgVal_bb = m_workspace->var("n_BkgNonHiggs_bb")->getVal();
  double bkgErr_bb =  m_workspace->var("n_BkgNonHiggs_bb")->getError();
  double SigSMVal_bb = m_workspace->function("n_SigSM_bb")->getVal();
  double SigSMErr_bb
    = m_workspace->function("n_SigSM_bb")->getPropagatedError(*fitResult);
  
  std::cout << "bkgVal_bb = " << bkgVal_bb << " +/- " << bkgErr_bb << std::endl;
  std::cout << "SigSMVal_bb = " << SigSMVal_bb << " +/- " << SigSMErr_bb
	    << std::endl;
  */

  
  //std::cout << "CONTINUUM: " 
  //	    << m_workspace->function("n_BkgNonHiggs_bbSideBand")->getVal() 
  //	    << " +/- "
  //	    << m_workspace->function("n_BkgNonHiggs_bbSideBand")->getError() 
  //	    << std::endl;
  
  





  //std::cout << "Post-fit parameter values" << std::endl;
  //printSet("nuisanceParameters", nuisanceParameters);
  //printSet("globalObservables", globalObservables);
  
  /*
  if ((m_config->getStr("AnalysisType")).EqualTo("NonResonant")) {
    
    double n_AllProcesses_bb
      = m_workspace->function("n_AllProcesses_bb")->getVal();
    double n_BkgNonHiggs_bb
      = m_workspace->function("n_BkgNonHiggs_bb")->getVal();
    double n_SigSM_bb = (m_workspace->function("n_SigSMggH_bb")->getVal() + 
			 m_workspace->function("n_SigSMVBF_bb")->getVal() +
			 m_workspace->function("n_SigSMWH_bb")->getVal() +
			 m_workspace->function("n_SigSMZH_bb")->getVal() +
			 m_workspace->function("n_SigSMttH_bb")->getVal() +
			 m_workspace->function("n_SigSMbbH_bb")->getVal());
    double n_SigBSM2H_bb = m_workspace->function("n_SigBSM2H_bb")->getVal();
    double n_bias_bb = 0.0;
    if (m_config->getBool("UseSystematics")) {
      n_bias_bb = m_workspace->function("n_Bias_bb")->getVal();
    }
    std::cout << "\nPrinting bb values for comparison" << std::endl;
    std::cout << "\tAllProcesses_bb = " << n_AllProcesses_bb << std::endl;
    std::cout << "\tBkgContinuum_bb = " << n_BkgNonHiggs_bb << std::endl;
    std::cout << "\tSigSM_bb        = " << n_SigSM_bb << std::endl;
    std::cout << "\tSigBSM2H_bb     = " << n_SigBSM2H_bb << std::endl;
    std::cout << "\tBias_bb         = " << n_bias_bb << std::endl;
    
    double n_AllProcesses_jj
      = m_workspace->function("n_AllProcesses_jj")->getVal();
    double n_BkgNonHiggs_jj
      = m_workspace->function("n_BkgNonHiggs_jj")->getVal();
    double n_SigSM_jj = (m_workspace->function("n_SigSMggH_jj")->getVal() + 
			 m_workspace->function("n_SigSMVBF_jj")->getVal() +
			 m_workspace->function("n_SigSMWH_jj")->getVal() +
			 m_workspace->function("n_SigSMZH_jj")->getVal() +
			 m_workspace->function("n_SigSMttH_jj")->getVal() +
			 m_workspace->function("n_SigSMbbH_jj")->getVal());
    double n_bias_jj = 0.0;
    if (m_config->getBool("UseSystematics")) {
      n_bias_jj = m_workspace->function("n_Bias_jj")->getVal();
    }
    std::cout << "\nPrinting jj values for comparison" << std::endl;
    std::cout << "\tAllProcesses_jj = " << n_AllProcesses_jj << std::endl;
    std::cout << "\tBkgContinuum_jj = " << n_BkgNonHiggs_jj << std::endl;
    std::cout << "\tSigSM_jj        = " << n_SigSM_jj << std::endl;
    std::cout << "\tBias_jj         = " << n_bias_jj << std::endl;
    std::cout << "\n" << std::endl;
  }
  else {
    double n_AllProcesses_bb
      = m_workspace->function("n_AllProcesses_bb")->getVal();
    double n_BkgNonHiggs_bb
      = m_workspace->function("n_BkgNonHiggs_bb")->getVal();
    double n_SigSM_bb = (m_workspace->function("n_SigSMggH_bb")->getVal() + 
			 m_workspace->function("n_SigSMVBF_bb")->getVal() +
			 m_workspace->function("n_SigSMWH_bb")->getVal() +
			 m_workspace->function("n_SigSMZH_bb")->getVal() +
			 m_workspace->function("n_SigSMttH_bb")->getVal() +
			 m_workspace->function("n_SigSMbbH_bb")->getVal() +
			 m_workspace->function("n_SigSMhh_bb")->getVal());
    
    double n_SigBSM2H_bb = m_workspace->function("n_SigBSM2H_bb")->getVal();
    std::cout << "\nPrinting bb values for comparison" << std::endl;
    std::cout << "\tAllProcesses_bb = " << n_AllProcesses_bb << std::endl;
    std::cout << "\tBkgContinuum_bb  = " << n_BkgNonHiggs_bb << std::endl;
    std::cout << "\tSigSM_bb        = " << n_SigSM_bb << std::endl;
    std::cout << "\tSigBSM2H_bb     = " << n_SigBSM2H_bb << std::endl;
    std::cout << "\n" << std::endl;
  }
  

  m_workspace->var("RNDM_MyyMODELING")->setVal(0);
  m_workspace->var("MyyMODELING")->setVal(0);
  std::cout << "CENTRAL BGM= " << varNLL->getVal() << std::endl;
  m_workspace->var("MyyMODELING")->setVal(1);
  std::cout << "PULLED BGM= " << varNLL->getVal() << std::endl;
  
  m_workspace->var("LUMI")->setVal(0);
  m_workspace->var("RNDM_LUMI")->setVal(0);
  std::cout << "CENTRAL LUMI= " << varNLL->getVal() << std::endl;
  m_workspace->var("LUMI")->setVal(1);
  std::cout << "PULLED LUMI= " << varNLL->getVal() << std::endl;
  */
  
  // Save a snapshot if requested:
  if (m_doSaveSnapshot) {
    TString textValPoI = fixPoI ? (Form("%d",(int)valPoI)) : "Free";
    m_workspace->saveSnapshot(Form("paramsProfilePoI%s", textValPoI.Data()),
			      *poiAndNuis);
  }
  
  // Plot the fit result if the user has set an output directory for plots:
  if (m_doPlot) {
    if (fixPoI && ((int)valPoI) == 1) plotFits("Mu1", datasetName);
    else if (fixPoI && ((int)valPoI) == 0) plotFits("Mu0", datasetName);
    else plotFits("MuFree", datasetName);
  }
  
  profiledValPoI = firstPoI->getVal();
  double nllValue = varNLL->getVal();
  delete varNLL;
  
  // Save names and values of nuisance parameters, globs, other parameters:
  storeParams(nuisanceParameters, m_mapNP);
  storeParams(globalObservables, m_mapGlobs);
  storeParams((RooArgSet*)m_workspace->set("nonSysParameters"), m_mapPars);
  
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  
  // Finish up, return NLL value.
  printer("DHTestStat: Fit has completed. Returning NLL.", false);
  return nllValue;
}

/**
   -----------------------------------------------------------------------------
   Get a map of global observable names to values from the most recent fit.
   @return - A map of global observable names and most recent fit values.
*/
std::map<std::string,double> DHTestStat::getGlobalObservables() {
  return m_mapGlobs;
}

/**
   -----------------------------------------------------------------------------
   Get the key for the value map.
   @param testStat - The test statistic.
   @param observed - True iff. observed.
   @param N - The sigma value (-2,-1,0,1,2).
*/
TString DHTestStat::getKey(TString testStat, bool observed, int N) {
  TString currKey = testStat;
  
  if (observed) currKey += "_obs";
  else currKey += "_exp";
  
  if (N > 0) currKey += Form("_p%d",N);
  else if (N < 0) currKey += Form("_n%d",N);
  
  return currKey;
}

/**
   -----------------------------------------------------------------------------
   Get the number of events in each category in the most recent pseudo-dataset.
   @return - A vector containing the weighted # events per category.
*/
std::vector<double> DHTestStat::getNEventsToys() {
  return m_numEventsPerCate;
}

/**
   -----------------------------------------------------------------------------
   Get a map of nuisance parameter names to values from the most recent fit.
   @return - A map of nuisance parameter names and most recent fit values.
*/
std::map<std::string,double> DHTestStat::getNuisanceParameters() {
  return m_mapNP;
}

/**
   -----------------------------------------------------------------------------
   Get a map of non-systematic parameter names to values from most recent fit.
   @return - A map of parameter names and most recent fit values.
*/
std::map<std::string,double> DHTestStat::getParameters() {
  return m_mapPars;
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of p0 based on the test statistic q0.
   @param q0 - The test statistic q0.
   @return - The value of p0.
*/
double DHTestStat::getP0FromQ0(double q0) {
  return (1 - ROOT::Math::gaussian_cdf(sqrt(fabs(q0))));
}

/**
   -----------------------------------------------------------------------------
   Calculate pB based on the standard deviation.
   @param N - The standard deviation.
   @return - The value of pB.
*/
double DHTestStat::getPbFromN(double N) {
  return (1 - ROOT::Math::gaussian_cdf(N));
}

/**
   -----------------------------------------------------------------------------
   Calculate the pB value based on qMu.
   @param qMu - The test statistic qMu.
   @param sigma - The sigma value...
   @param mu - The mu value... 
   @return - The value of pB.
*/
double DHTestStat::getPbFromQMu(double qMu, double sigma, double mu) {
  return (1 - ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qMu)));
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of pMu.
   @param qMu - The test statistic qMu.
   @return - The value of pMu.
*/
double DHTestStat::getPMuFromQMu(double qMu) {
  return (1 - ROOT::Math::gaussian_cdf(sqrt(fabs(qMu))));
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic q0 based on the nll.
   @param nllMu0 - NLL of a fit with signal strength 0;
   @param nllMuHat - NLL of a fit with profiled signal strength.
   @param muHat - Profiled signal strength.
   @return - The value of q0.
*/
double DHTestStat::getQ0FromNLL(double nllMu0, double nllMuHat, double muHat) {
  double q0 = (muHat < 0) ? 0 : (2 * (nllMu0 - nllMuHat));
  return q0;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic qMu based on the nll.
   @param nllMu - NLL of a fit with signal strength mu.
   @param nllMuHat - NLL of a fit with profiled signal strength.
   @param muHat - Profiled signal strength.
   @param muTest - Tested value of signal strength.
   @return - The value of qMu.
*/
double DHTestStat::getQMuFromNLL(double nllMu, double nllMuHat, double muHat,
				 double muTest) {
  double qMu = 0.0;
  if (muHat < muTest) qMu = 2 * (nllMu - nllMuHat);
  return qMu;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic qMuTilde based on the nll.
   @param nllMu - NLL of a fit with signal strength mu.
   @param nllMu0 - NLL of a fit with signal strength 0.
   @param nllMuHat - NLL of a fit with profiled signal strength.
   @param muHat - Profiled signal strength.
   @param muTest - Tested value of signal strength.
   @return - The value of qMuTilde.
*/
double DHTestStat::getQMuTildeFromNLL(double nllMu, double nllMu0,
				      double nllMuHat, double muHat,
				      double muTest) {
  double qMuTilde = 0;
  if (muHat <= 0) qMuTilde = 2 * (nllMu - nllMu0);
  else if (muHat > 0 && muHat <= muTest) qMuTilde = 2 * (nllMu - nllMuHat);
  else if (muHat > muTest) qMuTilde = 0;
  return qMuTilde;
}

/**
   -----------------------------------------------------------------------------
   Load the statistics files (p0 and CL) that were previously generated. If none
   are found, then create from scratch automatically.
*/
void DHTestStat::loadStatsFromFile() {
  
  // Load input p0 file:
  ifstream textP0;
  textP0.open(Form("%s/p0/p0_values_%s.txt", 
		   m_outputDir.Data(), m_anaType.Data()));
  
  // Load input CL file:
  ifstream textCL;
  textCL.open(Form("%s/CL/CL_values_%s.txt", 
		   m_outputDir.Data(), m_anaType.Data()));
  
  // If the input files don't exist, create from scratch:
  if (!textCL || !textP0) {
    calculateNewCL();
    calculateNewP0();
    return;
  }
  
  // Read p0 values:
  TString inName; double inExpP0; double inObsP0; 
  while (!textP0.eof()) {
    textP0 >> inName >> inExpP0 >> inObsP0;
    textP0.close();
  }
  m_calculatedValues[getKey("p0",1,0)] = inObsP0;
  m_calculatedValues[getKey("p0",0,0)] = inExpP0;
  
  // Read CL values:
  double inObsCL, inExpCLn2, inExpCLn1, inExpCL, inExpCLp1, inExpCLp2;
  while (!textCL.eof()) {
    textCL >> inName >> inObsCL >> inExpCLn2 >> inExpCLn1 >> inExpCL
	   >> inExpCLp1 >> inExpCLp2;
  }
  textCL.close();
  
  // save CL and CLs for later access:
  m_calculatedValues[getKey("CL", 0, -2)] = inExpCLn2;
  m_calculatedValues[getKey("CL", 0, -1)] = inExpCLn1;
  m_calculatedValues[getKey("CL", 0, 0)] = inExpCL;
  m_calculatedValues[getKey("CL", 0, 1)] = inExpCLp1;
  m_calculatedValues[getKey("CL", 0, 2)] = inExpCLp2;
  m_calculatedValues[getKey("CL", 1, 0)] = inObsCL;
  
  m_calculatedValues[getKey("CLs", 0, -2)] = getCLsFromCL(inExpCLn2);
  m_calculatedValues[getKey("CLs", 0, -1)] = getCLsFromCL(inExpCLn1);
  m_calculatedValues[getKey("CLs", 0, 0)] = getCLsFromCL(inExpCL);
  m_calculatedValues[getKey("CLs", 0, 1)] = getCLsFromCL(inExpCLp1);
  m_calculatedValues[getKey("CLs", 0, 2)] = getCLsFromCL(inExpCLp2);
  m_calculatedValues[getKey("CLs", 1, 0)] = getCLsFromCL(inObsCL);
}
/**
   -----------------------------------------------------------------------------
   Takes in a variable declaration such as "name[1,0,5]" and returns "name".
   @param varForm - The form of the variable declared.
   @return - The name of the variable without the rest of the expression.
*/
TString DHTestStat::nameOfVar(TString varForm) {
  TString name = varForm;
  name.Remove(name.First("["));
  return name;
}

/**
   -----------------------------------------------------------------------------
   Create a ratio plot:
   @param dataName - The name of the RooAbsData set in the workspace.
   @param pdfName - The name of the RooAbsPdf in the workspace.
   @param xMin - The minimum value of the observable range.
   @param xMax - The maximum value of the observable range.
   @param xBins - The number of bins for the observable.
   @return - A TGraphErrors to plot.
*/
TGraphErrors* DHTestStat::plotDivision(TString dataName, TString pdfName, 
				       TString obsName, double xMin, 
				       double xMax, double xBins){
  printer(Form("DHTestStat::plotDivision(%s, %s, %s, %f, %f, %f)",
	       dataName.Data(),pdfName.Data(),obsName.Data(),xMin,xMax,xBins),
	  false);
  
  RooRealVar *observable = m_workspace->var(obsName);
  RooAbsData *data = m_workspace->data(dataName);
  RooAbsPdf *pdf = m_workspace->pdf(pdfName); 
  double minOrigin = observable->getMin();
  double maxOrigin = observable->getMax();
  double nEvents = data->sumEntries();
    
  observable->setRange("fullRange", xMin, xMax);
  TH1F *originHist
    = (TH1F*)data->createHistogram("dataSub", *observable,
  				   RooFit::Binning(xBins, xMin, xMax));
  TGraphErrors *result = new TGraphErrors();
  double increment = (xMax - xMin) / ((double)xBins);
  
  RooAbsReal* intTot
    = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable),
				       RooFit::NormSet(*observable), 
				       RooFit::Range("fullRange"));
  double valTot = intTot->getVal();
  
  int pointIndex = 0;
  for (double i_m = xMin; i_m < xMax; i_m += increment) {
    observable->setRange(Form("range%2.2f",i_m), i_m, (i_m+increment));
    RooAbsReal* intCurr
      = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable), 
					 RooFit::NormSet(*observable), 
					 RooFit::Range(Form("range%2.2f",i_m)));
    double valCurr = intCurr->getVal();
    
    double currMass = i_m + (0.5*increment);
    double currPdfWeight = nEvents * (valCurr / valTot);
    TString varName = observable->GetName();
    double currDataWeight = data->sumEntries(Form("%s>%f&&%s<%f",varName.Data(),
						  i_m,varName.Data(),
						  (i_m+increment)));
    double currWeight = currDataWeight / currPdfWeight;
    if (currDataWeight == 0) currWeight = 1.0;
    result->SetPoint(pointIndex, currMass, currWeight);
    
    double currError = originHist->GetBinError(pointIndex+1) / currPdfWeight;
    result->SetPointError(pointIndex, 0.0, currError);
    pointIndex++;
  }
  observable->setMin(minOrigin);
  observable->setMax(maxOrigin);
  return result;
}

/**
   -----------------------------------------------------------------------------
   Plot the fits produced by the specified model.
   @param fitType - The type of fit.
   @param datasetName - The name of the dataset to be plotted. 
*/
void DHTestStat::plotFits(TString fitType, TString datasetName) {
  printer(Form("DHTestStat:plotFits(%s, %s)",fitType.Data(),datasetName.Data()),
	  false);
  
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  can->cd();
  TPad *pad1 = new TPad( "pad1", "pad1", 0.00, 0.33, 1.00, 1.00 );
  TPad *pad2 = new TPad( "pad2", "pad2", 0.00, 0.00, 1.00, 0.33 );
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.4);
  pad2->SetBorderMode(0);
  
  can->cd();
  pad1->Draw();
  pad2->Draw();
    
  // Loop over categories:
  std::vector<TString> cateNames = m_config->getStrV("CateNames");
  for (int i_c = 0; i_c < (int)cateNames.size(); i_c++) {
    pad1->cd();
    pad1->Clear();
    
    TString obsName = m_config->getStr(Form("OBS_%s",cateNames[i_c].Data()));
    obsName = nameOfVar(obsName);
    double obsMin = (*m_workspace->var(obsName)).getMin();
    double obsMax = (*m_workspace->var(obsName)).getMax();
    
    // Set the resonant analysis plot binning and axis scale to paper settings:
    int nBinsForPlot = (int)((obsMax - obsMin) / m_geVPerBin);
    
    // Plot everything on RooPlot:
    RooPlot* frame = (*m_workspace->var(obsName)).frame(nBinsForPlot);
    if ((m_workspace->pdf("model_"+cateNames[i_c]))) {
      (*m_workspace->data(Form("%s_%s",datasetName.Data(),cateNames[i_c].Data()))).plotOn(frame);
      if ((m_workspace->pdf("pdf_SigBSM2H_"+cateNames[i_c]))) {
	(*m_workspace->pdf("model_"+cateNames[i_c])).plotOn(frame, Components((*m_workspace->pdf("pdf_SigBSM2H_"+cateNames[i_c]))), LineColor(6));
      }
      if ((m_workspace->pdf("pdf_SigSM_"+cateNames[i_c]))) {
	(*m_workspace->pdf("model_"+cateNames[i_c])).plotOn(frame, Components((*m_workspace->pdf("pdf_SigSM_"+cateNames[i_c]))), LineColor(3));
      }
      if ((m_workspace->pdf("pdf_BkgNonHiggs_"+cateNames[i_c]))) {
	(*m_workspace->pdf("model_"+cateNames[i_c])).plotOn(frame, Components((*m_workspace->pdf("pdf_BkgNonHiggs_"+cateNames[i_c]))), LineColor(4));
      }
      (*m_workspace->pdf("model_"+cateNames[i_c])).plotOn(frame, LineColor(2));
    }
    TString xTitle =m_config->getStr(Form("OBSPrint_%s",cateNames[i_c].Data()));
    frame->SetXTitle(xTitle);
    frame->SetYTitle(Form("Events / %2.1f GeV", m_geVPerBin));
    frame->Draw();
    
    if (m_useLogScale) {
      gPad->SetLogy();
      frame->GetYaxis()->SetRangeUser(m_yMin, m_yMax);
    }
    
    //TLatex text; text.SetNDC(); text.SetTextColor(1);
    //text.DrawLatex(0.6, 0.89, cateNames[i_c]);
    TH1F *histDH = new TH1F("histDH", "histDH", 1, 0, 1);
    TH1F *histSH = new TH1F("histSH", "histSH", 1, 0, 1);
    TH1F *histBkg = new TH1F("histBkg", "histBkg", 1, 0, 1);
    TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
    histDH->SetLineColor(6);
    histSH->SetLineColor(3);
    histBkg->SetLineColor(4);
    histSig->SetLineColor(2);
    TLegend leg(0.60, 0.55, 0.92, 0.70);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.SetBorderSize(0);
    leg.AddEntry(histDH, "Di-Higgs", "l");
    leg.AddEntry(histSH, "Single Higgs", "l");
    leg.AddEntry(histBkg, "Non-resonant", "l");
    leg.AddEntry(histSig, "Sum", "l");
    leg.Draw("SAME");
    
    // Print ATLAS text on the plot:    
    TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
    t.SetTextFont(72); t.SetTextSize(0.05);
    t.DrawLatex(0.55, 0.88, "ATLAS");
    t.SetTextFont(42); t.SetTextSize(0.05);
    t.DrawLatex(0.67, 0.88, m_config->getStr("ATLASLabel"));
    //t.SetTextSize(0.04);
    t.DrawLatex(0.55, 0.82, 
		Form("#sqrt{s} = 13 TeV: #scale[0.7]{#int}Ldt = %2.1f fb^{-1}",
		     (m_config->getNum("AnalysisLuminosity")/1000.0)));
    TString printCateName 
      = m_config->getStr(Form("PrintCateName_%s", cateNames[i_c].Data()));
    t.DrawLatex(0.55, 0.76, printCateName);
    
    // Second pad on the canvas (Ratio Plot):
    pad2->cd();
    pad2->Clear();
    
    double ratioMin = 0.0;//-0.2
    double ratioMax = 2.0;//2.2
    TH1F *medianHist = new TH1F("median","median",nBinsForPlot,obsMin,obsMax);
    for (int i_b = 1; i_b <= nBinsForPlot; i_b++) {
      medianHist->SetBinContent(i_b, 1.0);
    }
    medianHist->SetLineColor(kRed);
    medianHist->SetLineWidth(2);
    medianHist->GetXaxis()->SetTitle(xTitle);
    medianHist->GetYaxis()->SetTitle("Data / Fit");
    medianHist->GetXaxis()->SetTitleOffset(0.95);
    medianHist->GetYaxis()->SetTitleOffset(0.7);
    medianHist->GetXaxis()->SetTitleSize(0.1);
    medianHist->GetYaxis()->SetTitleSize(0.1);
    medianHist->GetXaxis()->SetLabelSize(0.1);
    medianHist->GetYaxis()->SetLabelSize(0.1);
    medianHist->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    //medianHist->GetYaxis()->SetNdivisions(5);
    medianHist->GetYaxis()->SetNdivisions(4);
    medianHist->Draw();
    
    
    //medianHist->GetXaxis()->SetRangeUser(obsMin, obsMax);
    
    TGraphErrors* subData
      = plotDivision(Form("%s_%s",datasetName.Data(),cateNames[i_c].Data()),
		     Form("model_%s", cateNames[i_c].Data()), obsName, obsMin,
		     obsMax, nBinsForPlot);
    subData->Draw("EPSAME");
    
    TLine *line = new TLine();
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    //line->DrawLine(obsMin, 1.0, obsMax, 1.0); 
    line->SetLineWidth(1);
    line->SetLineStyle(2);
    line->DrawLine(obsMin, ((1.0+ratioMin)/2.0), obsMax, ((1.0+ratioMin)/2.0));
    line->DrawLine(obsMin, ((1.0+ratioMax)/2.0), obsMax, ((1.0+ratioMax)/2.0));
    subData->Draw("EPSAME");
        
    can->Print(Form("%s/fitPlot_%s_%s_%s.eps", m_plotDir.Data(),
		    m_anaType.Data(), fitType.Data(), cateNames[i_c].Data()));
    delete histDH;
    delete histSH;
    delete histBkg;
    delete histSig;
    delete frame;
  }
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Check whether the specified map entry exists.
    @param mapKey - The key for which we are finding a value:
    @return - True iff the categorization has been defined. 
*/
bool DHTestStat::mapValueExists(TString mapKey) {

  // Checks if there is a key corresponding to mapKey in the map: 
  bool nonExistent = (m_calculatedValues.find(mapKey) ==
		      m_calculatedValues.end());
  if (nonExistent) {
    std::cout << "DHTestStat: key " << mapKey << " not defined!" << std::endl;
  }
  return !nonExistent;
}

/**
   -----------------------------------------------------------------------------
   Prints a statement (if verbose) and exits (if fatal).
   @param statement - The statement to print.
   @param isFatal - True iff. this should trigger an exit command.
*/
void DHTestStat::printer(TString statement, bool isFatal) {
  if (m_config->getBool("Verbose") || isFatal) {
    std::cout << statement << std::endl;
  }
  if (isFatal) exit(0);
}

/**
   -----------------------------------------------------------------------------
   Print the names and values of parameters in a RooArgSet:
   @param setName - The name of the RooArgSet.
   @param set - The RooArgSet.
*/
void DHTestStat::printSet(TString setName, RooArgSet* set) {
  printer(Form("DHTestStat::printSet(%s)",setName.Data()), false);
  TIterator *iterSet = set->createIterator();
  RooRealVar *curr = NULL;
  // Bin each of the observables:
  while ((curr = (RooRealVar*)iterSet->Next())) {
    TString currLine = Form("\t %s = %f", ((TString)curr->GetName()).Data(),
			    curr->getVal());
    printer(currLine, false);
  }
}

/**
   -----------------------------------------------------------------------------
   Choose whether or not to save snapshots from profiling data.
   @param doSaveSnapshot - True iff you want to save snapshots in future fits.
*/
void DHTestStat::saveSnapshots(bool doSaveSnapshot) {
  m_doSaveSnapshot = doSaveSnapshot;
}

/**
   -----------------------------------------------------------------------------
   Set axis options for plots:
*/
void DHTestStat::setPlotAxis(bool useLogScale, double yMin, double yMax,
			     double GeVPerBin) {
  m_useLogScale = useLogScale;
  m_yMin = yMin;
  m_yMax = yMax;
  m_geVPerBin = GeVPerBin;
}

/**
   -----------------------------------------------------------------------------
   Set an output directory and enable plotting.
   @param directory - The output directory path.
*/
void DHTestStat::setPlotDirectory(TString directory) {
  m_plotDir = directory;
  m_doPlot = true;
}

/**
   -----------------------------------------------------------------------------
   Set the named parameter to a certain value and either fix or free it.
   @param paramName - The name of the fit parameter.
   @param paramVal - The new value of the fit parameter.
   @param doSetConstant - True iff the parameter should be set constant. 
*/
void DHTestStat::setParam(TString paramName, double paramVal,
			  bool doSetConstant) {
  m_paramValToSet[paramName] = paramVal;
  m_paramConstToSet[paramName] = doSetConstant;
  
  // Do this immediately?
  m_workspace->var(paramName)->setVal(paramVal);
  m_workspace->var(paramName)->setConstant(doSetConstant);
}

/**
   -----------------------------------------------------------------------------
   Store the parameter names and values.
   @param set - The RooArgSet to store.
   @param map - The map to hold the values.
*/
void DHTestStat::storeParams(RooArgSet *set, std::map<std::string,double>& map){
  map.clear();
  TIterator *iterSet = set->createIterator();
  RooRealVar *curr = NULL;
  while ((curr = (RooRealVar*)iterSet->Next())) {
    map[(std::string)curr->GetName()] = curr->getVal();
  }
}

/**
   -----------------------------------------------------------------------------
   Get the workspace.
   @return - A pointer to the class RooWorkspace object.
*/
RooWorkspace *DHTestStat::theWorkspace() {
  return m_workspace;
}

/**
   -----------------------------------------------------------------------------
   Get the ModelConfig.
   @return - A pointer to the class ModelConfig object.
*/

ModelConfig *DHTestStat::theModelConfig() {
  return m_mc;
}
