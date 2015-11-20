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
	    << "\n\t" << newOptions << "\n\t" << std::endl;
  
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
  m_dataForObsQ0 = (m_config->getBool("doBlind")) ? "asimovDataMu1" : "obsData";
  m_dataForObsQMu= (m_config->getBool("doBlind")) ? "asimovDataMu0" : "obsData";
  
  // Load from file if the pointer passed is NULL:
  if (newWorkspace == NULL) {
    inputFile
      = new TFile(Form("%s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
		       (m_config->getStr("masterOutput")).Data(),
		       m_jobName.Data(), m_anaType.Data()), "read");
    if (inputFile->IsOpen()) {
      std::cout << "DHTestStat: Loading workspace." << std::endl;
      m_workspace = (RooWorkspace*)inputFile->Get("combinedWS");
    }
    else {
      std::cout << "DHTestStat: Error loading file, accessing with WS tool."
		<< std::endl;
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
  setPlotAxis(false, 0.0, 100.0);

  // Create output directories:
  m_outputDir = Form("%s/%s/DHTestStat", 
		     (m_config->getStr("masterOutput")).Data(),
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
  std::cout << "DHTestStat: Calculating CLs" << std::endl;
  
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
  double expCLn2 = getCLFromQMu(expQMu, 0, -2);
  double expCLn1 = getCLFromQMu(expQMu, 0, -1);
  double expCLp1 = getCLFromQMu(expQMu, 0, 1);
  double expCLp2 = getCLFromQMu(expQMu, 0, 2);
  double expCL = getCLFromQMu(expQMu, 0, 0);
  double obsCL = getCLFromQMu(obsQMu, 1, 0);
  
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
  std::cout << "DHTestStat: calculating p0." << std::endl;
  
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
  std::cout << "DHTestStat: Create pseudodata with seed = " << seed 
	    << ", PoI = " << valPoI << std::endl;
  
  // Load the original parameters from profiling:
  m_workspace->loadSnapshot("paramsOrigin");
  
  RooSimultaneous* combPdf = (RooSimultaneous*)m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  RooArgSet* observables = (RooArgSet*)m_mc->GetObservables();
  RooArgSet* originValsNP
    = (RooArgSet*)m_mc->GetNuisanceParameters()->snapshot();
  RooRealVar* firstPoI = (RooRealVar*)m_mc->GetParametersOfInterest()->first();
  
  RooRandom::randomGenerator()->SetSeed(seed);
  statistics::constSet(nuisanceParameters, true);
  statistics::constSet(globalObservables, false);
  
  map<string,RooDataSet*> toyDataMap; 
  RooCategory *categories = (RooCategory*)m_workspace->obj("categories");
  TIterator *cateIter = combPdf->indexCat().typeIterator();
  RooCatType *cateType = NULL;
  RooDataSet *dataTemp[20];
  
  // Loop over all channels:
  int index = 0;
  // Previously this was commented and similar line below was uncommented
  statistics::randomizeSet(combPdf, globalObservables, seed); 
  statistics::constSet(globalObservables, true);
  
  //numEventsPerCate.clear();
  
  // Set the parameter of interest value and status:
  firstPoI->setVal(valPoI);
  firstPoI->setConstant(fixPoI);
  
  // Check if other parameter settings have been specified for toys:
  // WARNING! This overrides the randomization settings above!
  for (std::map<TString,double>::iterator iterParam = m_paramValToSet.begin();
       iterParam != m_paramValToSet.end(); iterParam++) {
    // CHeck that workspace contains parameter:
    if (m_workspace->var(iterParam->first)) {
      m_workspace->var(iterParam->first)->setVal(iterParam->second);
      m_workspace->var(iterParam->first)
	->setConstant(m_paramConstToSet[iterParam->first]);
    }
    else {
      m_workspace->Print("v");
      std::cout << "DHTestStat: Error! Parameter " << iterParam->first
		<< " not found in workspace! See printout above for clues..."
		<< std::endl;
      exit(0);
    }
  }
  
  // Iterate over the categories:
  while ((cateType = (RooCatType*)cateIter->Next())) {
    RooAbsPdf *currPDF = combPdf->getPdf(cateType->GetName());
    RooArgSet *currObs = currPDF->getObservables(observables);
    RooArgSet *currGlobs = currPDF->getObservables(globalObservables);
    RooRealVar *t = (RooRealVar*)currObs->first();
    
    //statistics::randomizeSet(currPDF, currGlobs, -1);
    //statistics::constSet(currGlobs, true);
    
    // If you want to bin the pseudo-data (speeds up calculation):
    if (m_options.Contains("Binned")) {
      currPDF->setAttribute("PleaseGenerateBinned");
      TIterator *iterObs = currObs->createIterator();
      RooRealVar *currObs = NULL;
      // Bin each of the observables:
      while ((currObs = (RooRealVar*)iterObs->Next())) {
	currObs->setBins(120);
      }
      dataTemp[index]
	= (RooDataSet*)currPDF->generate(*currObs, AutoBinned(true),
					 Extended(currPDF->canBeExtended()),
					 GenBinned("PleaseGenerateBinned"));
    }
    // Construct unbinned pseudo-data by default:
    else {
      dataTemp[index] = (RooDataSet*)currPDF->generate(*currObs,Extended(true));
    }
    
    toyDataMap[(std::string)cateType->GetName()] = dataTemp[index];
    //numEventsPerCate.push_back((double)dataTemp[index]->sumEntries());
    index++;
  }
  
  // Import the new data into the workspace:
  RooDataSet* pseudoData = new RooDataSet("toyData", "toyData", *observables, 
					  RooFit::Index(*categories),
					  RooFit::Import(toyDataMap));
  
  // release nuisance parameters:
  statistics::constSet(nuisanceParameters, false);
  
  // Import into the workspace:
  m_workspace->import(*pseudoData);
  
  return pseudoData;
}

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
   @param observed - True of observed stat., false if expected result.
   @param N - The sigma value (-2,-1,0,1,2). Use 0 for median.
   @return - The CLs value.
*/
double DHTestStat::getCLFromQMu(double qMu, bool observed, double N) {
  double CL = getCLFromCLs(getCLsFromQMu(qMu, observed, N));
  return CL;
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value using qMu and the type.
   @param qMu - The value of the test statistic.
   @param observed - True of observed stat., false if expected result.
   @param N - The sigma value (-2,-1,0,1,2). Use 0 for median.
   @return - The CLs value.
*/
double DHTestStat::getCLsFromQMu(double qMu, bool observed, double N) {
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
   @return - The NLL value.
*/
double DHTestStat::getFitNLL(TString datasetName, double valPoI, bool fixPoI,
			     double &profiledValPoI) { 
  std::cout << "DHTestStat: getFitNLL(" << datasetName << ", " << valPoI
	    << ", " << fixPoI << ")" << std::endl;
  
  RooAbsPdf* combPdf = m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  m_workspace->loadSnapshot("paramsOrigin");
  RooArgSet* origValNP = (RooArgSet*)m_workspace->getSnapshot("paramsOrigin");
  RooArgSet* poi = (RooArgSet*)m_mc->GetParametersOfInterest();
  RooRealVar* firstPoI = (RooRealVar*)poi->first();
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*nuisanceParameters);
  poiAndNuis->add(*poi);
  
  // Look for dataset. Create if non-existent & Asimov or requires binning.
  if (!m_workspace->data(datasetName)) {
    std::cout << "DHTestStat: Error! Requested data not available: " 
	      << datasetName << std::endl;
    exit(0);
  }
  
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  // the global observables should be fixed to the nominal values...
  statistics::constSet(globalObservables, true);
  
  firstPoI->setVal(valPoI);
  firstPoI->setConstant(fixPoI);
  
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
      std::cout << "DHTestStat: Error! Parameter " << iterParam->first
		<< " not found in workspace! See printout above for clues..."
		<< std::endl;
      exit(0);
    }
  }
  std::cout << "DHTestStat: Parameters have been set for fit." << std::endl;
  
  // The actual fit command:
  int status = 0; 
  RooNLLVar* varNLL = (RooNLLVar*)combPdf
    ->createNLL(*m_workspace->data(datasetName), Constrain(*nuisanceParameters),
		Extended(combPdf->canBeExtended()));
  RooFitResult *fitResult = statistics::minimize(varNLL, "", NULL, true);
  if (fitResult->status() != 0) m_allGoodFits = false;
  
  // Save a snapshot if requested:
  if (m_doSaveSnapshot) {
    TString muDHValue = fixPoI ? (Form("%d",(int)valPoI)) : "Free";
    m_workspace->saveSnapshot(Form("paramsProfileMu%s", muDHValue.Data()),
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
  
  // Save names and values of nuisance parameters:
  m_mapNP.clear();
  TIterator *iterNuis = nuisanceParameters->createIterator();
  RooRealVar *currNuis = NULL;
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    m_mapNP[(std::string)currNuis->GetName()] = currNuis->getVal();
  }
  
  // Save names and values of global observables:
  m_mapGlobs.clear();
  TIterator *iterGlobs = globalObservables->createIterator();
  RooRealVar *currGlob = NULL;
  while ((currGlob = (RooRealVar*)iterGlobs->Next())) {
    m_mapGlobs[(std::string)currGlob->GetName()] = currGlob->getVal();
  }
  
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  
  // Finish up, return NLL value.
  std::cout << "DHTestStat: Fit has competed. Returning NLL." << std::endl;
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
   Get a map of nuisance parameter names to values from the most recent fit.
   @return - A map of nuisance parameter names and most recent fit values.
*/
std::map<std::string,double> DHTestStat::getNuisanceParameters() {
  return m_mapNP;
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
   @param combWS - The combined workspace.
   @param fitType - The type of fit.
*/
void DHTestStat::plotFits(TString fitType, TString datasetName) {
  std::cout << "DHTestStat: Plot fit " << fitType << ", " << datasetName 
	    << std::endl;
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
  std::vector<TString> cateNames
    = m_config->getStrV(Form("cateNames%s", m_anaType.Data()));
  for (int i_c = 0; i_c < (int)cateNames.size(); i_c++) {
    //can->cd();
    //can->Clear();
    pad1->cd();
    pad1->Clear();
    
    TString obsName = m_anaType.EqualTo("NonRes") ? 
      Form("m_yy_%s", cateNames[i_c].Data()) : 
      Form("m_bbyy_%s", cateNames[i_c].Data());
    
    // Set the resonant analysis plot binning and axis scale to paper settings:
    int nGeVPerBin = 50; 
    if (cateNames[i_c].EqualTo("ResonantSR")) {
      nGeVPerBin = 5;
      m_yMin = 0.002;
      m_yMax = 40;
    }
    else if (cateNames[i_c].EqualTo("ResonantCR")) {
      nGeVPerBin = 20;
      m_yMin = 0.3;
      m_yMax = 40;
    }
    int nBinsForPlot = (int)(((*m_workspace->var(obsName)).getMax() - 
			      (*m_workspace->var(obsName)).getMin()) /
			     nGeVPerBin);
    
    RooPlot* frame = (*m_workspace->var(obsName)).frame(nBinsForPlot);
    (*m_workspace->data(Form("%s_%s",datasetName.Data(),cateNames[i_c].Data())))
      .plotOn(frame);
    
    (*m_workspace->pdf("model_"+cateNames[i_c]))
      .plotOn(frame,Components((*m_workspace->pdf("sigPdfDH_"+cateNames[i_c]))),
	      LineColor(6));
    (*m_workspace->pdf("model_"+cateNames[i_c]))
      .plotOn(frame,Components((*m_workspace->pdf("sigPdfSH_"+cateNames[i_c]))),
	      LineColor(3));
    (*m_workspace->pdf("model_"+cateNames[i_c]))
      .plotOn(frame, Components((*m_workspace->pdf("bkgPdf_"+cateNames[i_c]))), 
	      LineColor(4));
    (*m_workspace->pdf("model_"+cateNames[i_c])).plotOn(frame, LineColor(2));
    
    TString xTitle = m_anaType.EqualTo("NonRes") ? 
      "M_{#gamma#gamma} [GeV]" : "Constrained M_{#gamma#gammajj} [GeV]";
    frame->SetXTitle(xTitle);
    frame->SetYTitle(Form("Events / %d GeV", nGeVPerBin));
    frame->Draw();
    
    if (m_useLogScale) {
      gPad->SetLogy();
      frame->GetYaxis()->SetRangeUser(m_yMin, m_yMax);
    }
    
    TLatex text; text.SetNDC(); text.SetTextColor(1);
    text.DrawLatex(0.6, 0.89, cateNames[i_c]);
    TH1F *histDH = new TH1F("histDH", "histDH", 1, 0, 1);
    TH1F *histSH = new TH1F("histSH", "histSH", 1, 0, 1);
    TH1F *histBkg = new TH1F("histBkg", "histBkg", 1, 0, 1);
    TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
    histDH->SetLineColor(6);
    histSH->SetLineColor(3);
    histBkg->SetLineColor(4);
    histSig->SetLineColor(2);
    TLegend leg(0.6, 0.71, 0.92, 0.86);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.SetBorderSize(0);
    leg.AddEntry(histDH, "Di-Higgs", "l");
    leg.AddEntry(histSH, "Single Higgs", "l");
    leg.AddEntry(histBkg, "Non-resonant", "l");
    leg.AddEntry(histSig, "Sum", "l");
    leg.Draw("SAME");
    
    pad2->cd();
    pad2->Clear();
    
    double obsMin = (*m_workspace->var(obsName)).getMin();
    double obsMax = (*m_workspace->var(obsName)).getMax();
    
    double ratioMin = 0.0; double ratioMax = 2.0;
    TGraphErrors* subData
      = plotDivision(Form("%s_%s",datasetName.Data(),cateNames[i_c].Data()),
		     Form("model_%s", cateNames[i_c].Data()), obsName, obsMin,
		     obsMax, nBinsForPlot);
    subData->GetYaxis()->SetTitle("Data / Fit");
    subData->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    subData->GetXaxis()->SetTitleOffset(0.95);
    subData->GetYaxis()->SetTitleOffset(0.7);
    subData->GetXaxis()->SetTitleSize(0.1);
    subData->GetYaxis()->SetTitleSize(0.1);
    subData->GetXaxis()->SetLabelSize(0.1);
    subData->GetYaxis()->SetLabelSize(0.1);
    subData->GetYaxis()->SetNdivisions(4);
    subData->GetXaxis()->SetRangeUser(obsMin, obsMax);
    subData->Draw("AEP");
    
    TLine *line = new TLine();
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->DrawLine(obsMin, 1.0, obsMax, 1.0); 
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
void DHTestStat::setPlotAxis(bool useLogScale, double yMin, double yMax) {
  m_useLogScale = useLogScale;
  m_yMin = yMin;
  m_yMax = yMax;
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
}
