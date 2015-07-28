////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHTestStat.cxx                                                      //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 06/07/2015                                                          //
//                                                                            //
//  This class allows the user to calculate p0, CL, and CLs based on an input //
//  workspace.                                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHTestStat.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the DHTestStat class. 
   @param newJobName - The name of the job
   @param newDHSignal - The Di-Higgs signal to incorporate in the model.
   @param newOptions - The job options ("New", "FromFile"), etc.
   @param newWorkspace - The workspace with the model for the test stats. 
*/
DHTestStat::DHTestStat(TString newJobName, TString newDHSignal,
		       TString newCateScheme, TString newOptions,
		       RooWorkspace *newWorkspace) {
  std::cout << "DHTestStat: Initializing...\n\t" << newJobName << "\n\t"
	    << newDHSignal << "\n\t" << newCateScheme << "\n\t" 
	    << newOptions << "\n\t" << std::endl;
  
  // Assign input variables: 
  m_cateScheme = newCateScheme;
  m_jobName = newJobName;
  m_DHSignal = newDHSignal;
  m_options = newOptions;
  
  // Start with a clean class:
  clearData();
  
  m_anaType = DHAnalysis::getAnalysisType(m_DHSignal);
    
  // Use Asimov data if the analysis is blind.
  m_dataForObs = (DHAnalysis::doBlind) ?
    Form("asimovDataMu0_%s", m_anaType.Data()) :
    Form("obsData_%s", m_anaType.Data());
  m_dataForExp = Form("asimovDataMu0_%s", m_anaType.Data());
  
  // Load from file if the pointer passed is NULL:
  if (newWorkspace == NULL) {
    inputFile
      = new TFile(Form("%s/%s/DHWorkspace/rootfiles/workspaceDH.root",
		       DHAnalysis::masterOutput.Data(), m_jobName.Data()),
		  "read");
    if (inputFile->IsOpen()) {
      std::cout << "DHTestStat: Loading workspace." << std::endl;
      m_workspace = (RooWorkspace*)inputFile->Get("combinedWS");
    }
    else {
      std::cout << "DHTestStat: Error loading file, accessing with WS tool."
		<< std::endl;
      // Load the workspace from the nominal location.
      DHWorkspace *m_dhws = new DHWorkspace(m_jobName,m_cateScheme,"FromFile");
      m_workspace = m_dhws->getCombinedWorkspace();
    }
  }
  // Use the workspace passed to the class constructor:
  else {
    m_workspace = newWorkspace;
  }
  
  // Based on resonant or nonresonant:
  m_mc
    = (ModelConfig*)m_workspace->obj(Form("modelConfig_%s", m_anaType.Data()));
  
  // Map storing all calculations:
  m_calculatedValues.clear();
  
  // Create output directories:
  m_outputDir = Form("%s/%s/TestStat/", DHAnalysis::masterOutput.Data(), 
		     m_jobName.Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  system(Form("mkdir -vp %s/CL/", m_outputDir.Data()));
  system(Form("mkdir -vp %s/p0/", m_outputDir.Data()));

  // Make new or load old values:
  if (m_options.Contains("FromFile")) {
    loadStatsFromFile();
  }
  
  std::cout << "DHTestStat: Initialized Successfully!" << std::endl;
  return;
}

/**
   -----------------------------------------------------------------------------
   Get the value of one of the test statistics.
   @param testStat - the test stat. name (p0, CL, CLs).
   @param observed - true iff observed, false if expected. 
   @param N - the standard deviation (-2, -1, 0, +1, +2). 
   @returns - the test statistic value.
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
  double nllMu1Obs = getFitNLL(m_dataForObs, 1.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObs, 1.0, false, muHatObs);
  double obsQMu = getQMuFromNLL(nllMu1Obs, nllMuHatObs, muHatObs, 1);
  
  // Calculate expected qmu:
  double muHatExp = 0.0;
  double nllMu1Exp = getFitNLL(m_dataForExp, 1.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExp, 0.0, false, muHatExp);
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
		   m_outputDir.Data(), m_DHSignal.Data()));
  textCL << m_DHSignal << " " << obsCL << " " << expCLn2 << " " << expCLn1
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
  
  // save CL and CLs for later access:
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
  double nllMu0Obs = getFitNLL(m_dataForObs, 0.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObs, 0.0, false, muHatObs);
  double obsQ0 = getQ0FromNLL(nllMu0Obs, nllMuHatObs, muHatObs);
  
  // Calculate expected q0:
  double muHatExp = 0.0;
  double nllMu0Exp = getFitNLL(m_dataForExp, 0.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExp, 0.0, false, muHatExp);
  double expQ0 = getQ0FromNLL(nllMu0Exp, nllMuHatExp, muHatExp);
  
  // Calculate p0 from q0:
  double expP0 = getP0FromQ0(expQ0);
  double obsP0 = getP0FromQ0(obsQ0);
  
  // Write p0 values to file:
  ofstream textP0;
  textP0.open(Form("%s/p0/p0_values_%s.txt", 
		   m_outputDir.Data(), m_DHSignal.Data()));
  textP0 << m_DHSignal << " " << expP0 << " " << obsP0 << std::endl;
  textP0.close();
  
  // Print summary:
  std::cout << "\n  Expected p0 = " << expP0 << std::endl;
  std::cout << "  Observed p0 = " << obsP0 << std::endl;
  if (fitsAllConverged()) {
    std::cout << "All good fits? True\n" << std::endl;
  }
  else {
    std::cout << "All good fits? False\n" << std::endl;
  }
  
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
  m_namesGlobs.clear();
  m_namesNP.clear();
  m_valuesGlobs.clear();
  m_valuesNP.clear();
  
  //m_nBins = 240;

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
  m_setParamConsts.clear();
  m_setParamNames.clear();
  m_setParamVals.clear();
}

/**
   -----------------------------------------------------------------------------
   Create Asimov data for the statistical model, using a fit to observed data
   for the shape and normalizaiton of the background.
   @param cateWS - the current category workspace.
   @param valMuDH - the value of the di-Higgs signal strength to use.
   @param valMuSH - the value of the single Higgs signal strength to use.
   @returns - A RooDataSet with Asimov data, also imports in workspace. 
*/
RooDataSet* DHTestStat::createAsimovData(int valMuDH, int valMuSH) {
  std::cout << "DHWorkspace: Creating Asimov data, mu_DH = " << valMuDH
	    << " and mu_SH = " << valMuSH << std::endl;
  
  // Set mu_DH and mu_SH to the specified values:
  RooRealVar *poi = m_workspace->var("mu_DH");
  double initialMuDH = poi->getVal();
  double initialMuSH = m_workspace->var("mu_DH")->getVal();
  poi->setVal(valMuDH);
  poi->setConstant(true);
  m_workspace->var("mu_DH")->setVal(valMuSH);
  m_workspace->var("mu_DH")->setConstant(true);
  
  // Also use m_DHSignal to set the resonance mass.
  if (m_DHSignal.EqualTo("Res")) {
    int mediatorMass = DHAnalysis::getMediatorMass(m_DHSignal);
    if (mediatorMass > 0) {
      m_workspace->var("resMassDH")->setVal(mediatorMass);
      m_workspace->var("resMassDH")->setConstant(true);
    }
  }
  
  RooDataSet *asimovData = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*m_workspace->pdf(Form("combinedPdf_%s", m_anaType.Data())), *m_workspace->set(Form("observables_%s", m_anaType.Data())));
  TString asimovName = Form("asimovDataMu%d_%s", valMuDH, m_DHSignal.Data());
  asimovData->SetNameTitle(asimovName, asimovName);
  
  m_workspace->import(*asimovData);
  m_workspace->var("mu_DH")->setVal(initialMuDH);
  m_workspace->var("mu_SH")->setVal(initialMuSH);
  m_workspace->var("mu_DH")->setConstant(false);
  m_workspace->var("mu_SH")->setConstant(true);
  std::cout << "DHWorkspace: Asimov data has " << asimovData->sumEntries() 
	    << " entries" << std::endl;
  return asimovData;
}

/**
   -----------------------------------------------------------------------------
   Create an Asimov dataset using the name.
   @param datasetName - the name of the Asimov data set.
   @returns - A RooDataSet with Asimov data, also imports in workspace.
*/
RooDataSet* DHTestStat::createAsimovData(TString datasetName) {
  if (datasetName.Contains("Mu1")) return createAsimovData(1.0, 1.0);
  else if (datasetName.Contains("Mu0")) return createAsimovData(0.0, 1.0);
  else {
    std::cout << "DHTestStat: Asimov data creation error for " << datasetName
	      << std::endl;
    return NULL;
  }
}


/*
   -----------------------------------------------------------------------------
   Creates a binned dataset from a specified unbinned dataset.
   @param unbinnedName - The name of the dataset to be binned.
   @param nBins - The number of bins for the binning procedure. 
   @returns - A binned dataset.

RooDataSet* DHTestStat::createBinnedData(TString unbinnedName, int nBins) {
  
  RooDataSet* dataUnbinned = NULL;
  if (m_workspace->data(datasetName)) {
    dataUnbinned = m_workspace->data(unbinnedName);
  }
  // Create Asimov data if it is missing:
  else if (datasetName.Contains("asimovData")) {
    dataUnbinned = createAsimovData(datasetName);
  }
  else {
    std::cout << "DHTestStat: Error! Requested data unavailable for binning: " 
	      << datasetName << std::endl;
    exit(0);
  }
  
  // Load the RooCategory object from the workspace:
  RooCategory *categories
    = m_workspace->var(Form("categories_%s", currAna.Data()));
  for (int currCateIndex = 0; currCateIndex < categories->numTypes();
       currCateIndex++){
    
    // Create a binned observed data set:
    RooArgSet* obsPlusWt = new RooArgSet();
    RooRealVar wt("wt","wt",1);
    obsPlusWt->add(wt);
    TString obsName = m_anaType.EqualTo("NonRes") ? 
      Form("m_yy_%s",currCateName.Data()) : 
      Form("m_bbyy_%s",currCateName.Data());
    obsPlusWt->add(*m_observable);
    
    // Create a histogram to store binned data:
    TH1F* dataHist = new TH1F("dataHist", "dataHist", nBins, 
			      m_observable->getMin(), m_observable->getMax());
    for (int i_e = 0; i_e < dataUnbinned->numEntries(); i_e++) {
      dataUnbinned->get(i_e);
      dataHist->Fill(m_observable->getVal());
    }
    
    // Fill obsdatabinned dataset with binned data:
    TString unbinnedName = dataUnbinned->GetName();
    TString binnedName = Form("%s_binned", unbinnedName.Data());
    RooDataSet *obsDataBinned = new RooDataSet(binnedName,binnedName,*obsPlusWt,
					       RooFit::WeightVar(wt));
    for (int i_b = 1; i_b <= dataHist->GetNbinsX(); i_b++) {
      double massVal = dataHist->GetBinCenter(i_b);
      double weightVal = dataHist->GetBinContent(i_b);
      m_observable->setVal(massVal);
      wt.setVal(weightVal);
      obsDataBinned->add(RooArgSet(*m_observable, wt), weightVal);
    }
  }
  
  RooDataSet* binnedData = new RooDataSet(Form("%s_binned",unbinnedName.Data()),
					  Form("%s_binned",unbinnedName.Data()),
					  *args, Index(*categories), 
					  Import(dataMap),
					  WeightVar(wt));
  return binnedData
}
*/

/**
   -----------------------------------------------------------------------------
   Check if all of the fits done by this class have converged.
   @returns - true iff. all of the fits have been successfully convergent.
*/
bool DHTestStat::fitsAllConverged() { 
  return m_allGoodFits;
}

/**
   -----------------------------------------------------------------------------
   Implements the functional form of qMu.
   @param x - the value of the test statistic.
   @returns - the value of the asymptotic test statistic distribution.
*/
double DHTestStat::functionQMu(double x) {
  // This corresponds to the "special case" of mu=mu'
  double result = TMath::Exp(-1*x/2.0) / (2.0*sqrt(2.0*TMath::Pi()*x));
  return result;
}

/**
   -----------------------------------------------------------------------------
   Implements the functional form of qMuTilde.
   @param x - the value of the test statistic.
   @param asimovTestStat - the test stat value on Asimov data with mu=0 but
                           fitting under mu=1 hypothesis.
   @returns - the value of the asymptotic test statistic distribution.
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
   @param CLs - the CLs value to convert to CL.
   @returns - the corresponding CL value.
*/
double DHTestStat::getCLFromCLs(double CLs) {
  return (1.0 - CLs);
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value from CL.
   @param CL - the CL value to convert to CLs.
   @returns - the corresponding CLs value.
*/
double DHTestStat::getCLsFromCL(double CL) {
  return (1.0 - CL);
}

/**
   -----------------------------------------------------------------------------
   Get the CL value using qMu and the type.
   @param qMu - the value of the test statistic.
   @param observed - true of observed stat., false if expected result.
   @param N - the sigma value (-2,-1,0,1,2). Use 0 for median.
   @returns - the CLs value.
*/
double DHTestStat::getCLFromQMu(double qMu, bool observed, double N) {
  double CL = getCLFromCLs(getCLsFromQMu(qMu, observed, N));
  return CL;
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value using qMu and the type.
   @param qMu - the value of the test statistic.
   @param observed - true of observed stat., false if expected result.
   @param N - the sigma value (-2,-1,0,1,2). Use 0 for median.
   @returns - the CLs value.
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
   @param datasetName - the name of the dataset in the workspace.
   @param muVal - the mu value to fix.
   @param fixMu - true if mu should be fixed to the specified value.
   @param &profiledMu - the profiled value of mu (passed by reference)
   @returns - the nll value.
*/
double DHTestStat::getFitNLL(TString datasetName, double muVal, bool fixMu,
			     double &profiledMu) { 
  std::cout << "DHTestStat: getFitNLL(" << datasetName << ", " << muVal
	    << ", " << fixMu << ")" << std::endl;
  
  RooAbsPdf* combPdf = m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  m_workspace->loadSnapshot(Form("paramsOrigin_%s", m_anaType.Data()));
  RooArgSet* origValNP 
    = (RooArgSet*)m_workspace->getSnapshot(Form("paramsOrigin_%s",
						m_anaType.Data()));
  RooArgSet* poi = (RooArgSet*)m_mc->GetParametersOfInterest();
  RooRealVar* firstPoI = (RooRealVar*)poi->first();
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*nuisanceParameters);
  poiAndNuis->add(*poi);
  
  // Look for dataset. Create if non-existent & Asimov or requires binning.
  RooAbsData *data = NULL;
  if (m_workspace->data(datasetName)) data = m_workspace->data(datasetName);
  // Create binned dataset:
  //else if (m_options.Contains("binned")) {
  //data = createBinnedData(datasetName, m_nBins);
  //}
  // Create Asimov data:
  else if (datasetName.Contains("asimovData")) {
    data = createAsimovData(datasetName);
  }
  else {
    std::cout << "DHTestStat: Error! Requested data not available: " 
	      << datasetName << std::endl;
    exit(0);
  }
  
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  // the global observables should be fixed to the nominal values...
  statistics::constSet(globalObservables, true);
  
  firstPoI->setVal(muVal);
  firstPoI->setConstant(fixMu);
  
  // Check if parameters have been specified during fit:
  for (int i_p = 0; i_p < (int)m_setParamNames.size(); i_p++) {
    std::cout << "i_p=" << i_p << ", " << m_setParamNames[i_p] << std::endl;
    m_workspace->var(m_setParamNames[i_p])->setVal(m_setParamVals[i_p]);
    m_workspace->var(m_setParamNames[i_p])->setConstant(m_setParamConsts[i_p]);
  }
  
  // Iterate over SM mu values and fix all to 1:
  RooArgSet* muConstants = (RooArgSet*)m_workspace->set(Form("muSHConstants_%s",m_anaType.Data()));
  TIterator *iterMuConst = muConstants->createIterator();
  RooRealVar *currMuConst = NULL;
  
  while ((currMuConst = (RooRealVar*)iterMuConst->Next())) {
    std::cout << "DHTestStat: Setting " << currMuConst->GetName()
	      << " constant." << std::endl;
    currMuConst->setVal(1.0);
    currMuConst->setConstant(true);
  }
  
  int status = 0; 
  RooNLLVar* varNLL = (RooNLLVar*)combPdf->createNLL(*data, Constrain(*nuisanceParameters), Extended(combPdf->canBeExtended()));
  statistics::minimize(status, varNLL, "", NULL, false);
  if (status != 0) {
    m_allGoodFits = false;
  }
  
  // Save a snapshot if requested:
  if (m_doSaveSnapshot) {
    TString muDHValue = fixMu ? (Form("%d",(int)muVal)) : "Free";
    m_workspace->saveSnapshot(Form("paramsProfileMu%s", muDHValue.Data()),
			      *poiAndNuis);
  }
  
  // Plot the fit result if the user has set an output directory for plots:
  if (m_doPlot) {
    if (fixMu && ((int)muVal) == 1 ) plotFits("Mu1", datasetName);
    else if (fixMu && ((int)muVal) == 0 ) plotFits("Mu0", datasetName);
    else plotFits("MuFree", datasetName);
  }
  
  profiledMu = firstPoI->getVal();
  double nllValue = varNLL->getVal();
  delete varNLL;
  
  // Save names and values of nuisance parameters:
  m_namesNP.clear();
  m_valuesNP.clear();
  TIterator *iterNuis = nuisanceParameters->createIterator();
  RooRealVar *currNuis = NULL;
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    m_namesNP.push_back((std::string)currNuis->GetName());
    m_valuesNP.push_back(currNuis->getVal());
  }
  
  // Save names and values of global observables:
  m_namesGlobs.clear();
  m_valuesGlobs.clear();
  TIterator *iterGlobs = globalObservables->createIterator();
  RooRealVar *currGlob = NULL;
  while ((currGlob = (RooRealVar*)iterGlobs->Next())) {
    m_namesGlobs.push_back((std::string)currGlob->GetName());
    m_valuesGlobs.push_back(currGlob->getVal());
  }
    
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  return nllValue;
}

/**
   -----------------------------------------------------------------------------
   Get a vector of global observable names from the most recent fit.
*/
std::vector<std::string> DHTestStat::getGlobsNames() {
  return m_namesGlobs;
}

/**
   -----------------------------------------------------------------------------
   Get a vector of global observable values from the most recent fit.
*/
std::vector<double> DHTestStat::getGlobsValues() {
  return m_valuesGlobs;
}

/**
   -----------------------------------------------------------------------------
   Get the key for the value map.
   @param testStat - the test statistic.
   @param observed - true iff. observed.
   @param N - the sigma value (-2,-1,0,1,2).
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
   Get a vector of nuisance parameter names from the most recent fit.
*/
std::vector<std::string> DHTestStat::getNPNames() {
  return m_namesNP;
}

/**
   -----------------------------------------------------------------------------
   Get a vector of nuisance parameter values from the most recent fit.
*/
std::vector<double> DHTestStat::getNPValues() {
  return m_valuesNP;
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of p0 based on the test statistic q0.
   @param q0 - the test statistic q0.
   @returns - the value of p0.
*/
double DHTestStat::getP0FromQ0(double q0) {
  double p0 = 1 - ROOT::Math::gaussian_cdf(sqrt(fabs(q0)));
  return p0;
}

/**
   -----------------------------------------------------------------------------
   Calculate pB based on the standard deviation.
   @param N - the standard deviation.
   @returns - the value of pB.
*/
double DHTestStat::getPbFromN(double N) {
  double pB = 1 - ROOT::Math::gaussian_cdf(N);
  return pB;
}

/**
   -----------------------------------------------------------------------------
   Calculate the pB value based on qMu.
   @param qMu - the test statistic qMu.
   @param sigma - the sigma value...
   @param mu - the mu value... 
   @returns - the value of pB.
*/
double DHTestStat::getPbFromQMu(double qMu, double sigma, double mu) {
  double pB = 1 - ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qMu));
  return pB;
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of pMu.
   @param qMu - the test statistic qMu.
   @returns - the value of pMu.
*/
double DHTestStat::getPMuFromQMu(double qMu) {
  double pMu = 1 - ROOT::Math::gaussian_cdf(sqrt(fabs(qMu)));
  return pMu;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic q0 based on the nll.
   @param nllMu0 - nll of a fit with signal strength 0;
   @param nllMuHat - nll of a fit with profiled signal strength.
   @param muHat - profiled signal strength.
   @returns - the value of q0.
*/
double DHTestStat::getQ0FromNLL(double nllMu0, double nllMuHat, double muHat) {
  double q0 = (muHat < 0) ? 0 : (2 * (nllMu0 - nllMuHat));
  return q0;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic qMu based on the nll.
   @param nllMu - nll of a fit with signal strength mu.
   @param nllMuHat - nll of a fit with profiled signal strength.
   @param muHat - profiled signal strength.
   @param muTest - tested value of signal strength.
   @returns - the value of qMu.
*/
double DHTestStat::getQMuFromNLL(double nllMu, double nllMuHat, double muHat,
				 double muTest) {
  double qMu = 0;
  if (muHat < muTest) qMu = 2 * (nllMu - nllMuHat);
  else qMu = 0.0;
  return qMu;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic qMuTilde based on the nll.
   @param nllMu - nll of a fit with signal strength mu.
   @param nllMu0 - nll of a fit with signal strength 0.
   @param nllMuHat - nll of a fit with profiled signal strength.
   @param muHat - profiled signal strength.
   @param muTest - tested value of signal strength.
   @returns - the value of qMuTilde.
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
		   m_outputDir.Data(), m_DHSignal.Data()));
  
  // Load input CL file:
  ifstream textCL;
  textCL.open(Form("%s/CL/CL_values_%s.txt", 
		   m_outputDir.Data(), m_DHSignal.Data()));
  
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
  m_calculatedValues[getKey("CL",0,-2)] = inExpCLn2;
  m_calculatedValues[getKey("CL",0,-1)] = inExpCLn1;
  m_calculatedValues[getKey("CL",0,0)] = inExpCL;
  m_calculatedValues[getKey("CL",0,1)] = inExpCLp1;
  m_calculatedValues[getKey("CL",0,2)] = inExpCLp2;
  m_calculatedValues[getKey("CL",1,0)] = inObsCL;
  
  m_calculatedValues[getKey("CLs",0,-2)] = getCLsFromCL(inExpCLn2);
  m_calculatedValues[getKey("CLs",0,-1)] = getCLsFromCL(inExpCLn1);
  m_calculatedValues[getKey("CLs",0,0)] = getCLsFromCL(inExpCL);
  m_calculatedValues[getKey("CLs",0,1)] = getCLsFromCL(inExpCLp1);
  m_calculatedValues[getKey("CLs",0,2)] = getCLsFromCL(inExpCLp2);
  m_calculatedValues[getKey("CLs",1,0)] = getCLsFromCL(inObsCL);
}

/**
   -----------------------------------------------------------------------------
   Plot the fits produced by the specified model.
   @param combWS - the combined workspace.
   @param fitType - the type of fit.
*/
void DHTestStat::plotFits(TString fitType, TString datasetName) {
  std::cout << "DHTestStat: Plot final fits for " << fitType << std::endl;
  std::cout << "DEBUG1" << std::endl;
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  // loop over categories:
  int nCategories = DHAnalysis::getNumCategories(m_cateScheme, m_anaType);
  for (int i_c = 0; i_c < nCategories; i_c++) {
    can->Clear();
    std::cout << "DEBUG2" << std::endl;
    TString currCateName
      = DHAnalysis::cateIndexToName(m_cateScheme, m_anaType, i_c);
    TString obsName = m_anaType.EqualTo("NonRes") ? 
      Form("m_yy_%s", currCateName.Data()) : 
      Form("m_bbyy_%s", currCateName.Data());
    RooPlot* frame =  (*m_workspace->var(obsName)).frame(50);
    TString cutName = Form("categories_%s==categories_%s::%s", m_anaType.Data(),
			   m_anaType.Data(), currCateName.Data());
    std::cout << "DEBUG3" << std::endl;
    RooDataSet *currData
      =(RooDataSet*)(m_workspace->data(Form("%s", datasetName.Data())));
    std::cout << "DEBUG3.1" << std::endl;
    //currData->plotOn(frame, RooFit::Cut(cutName));
    /*
    //RooCategory *categories
    //=(RooCategory*)(m_workspace->var(Form("categories_%s",m_anaType.Data())));
    m_workspace->Print("v");
    std::cout << "DEBUG3.2" << std::endl;
    RooArgSet tempSet = m_workspace->allCats();
    tempSet.Print("v");
    std::cout << "DEBUG3.3" << std::endl;
    RooCategory *categories
      =(RooCategory*)(m_workspace->cat(Form("categories_%s",m_anaType.Data())));

    std::cout << "DEBUG4" << std::endl;
    (*m_workspace->pdf(Form("model_%s",m_anaType.Data())))
      .plotOn(frame, Components((*m_workspace->pdf("sigPdfDH_"+currCateName))),
	      LineColor(6));
    (*m_workspace->pdf(Form("model_%s",m_anaType.Data())))
      .plotOn(frame, Components((*m_workspace->pdf("sigPdfSH_"+currCateName))),
	      LineColor(3));
    (*m_workspace->pdf(Form("model_%s",m_anaType.Data())))
      .plotOn(frame, Components((*m_workspace->pdf("bkgPdf_"+currCateName))), 
	      LineColor(4));
    (*m_workspace->pdf(Form("model_%s",m_anaType.Data())))
      .plotOn(frame, LineColor(2));
    
    std::cout << "DEBUG5" << std::endl;
    //double chi2 = frame->chiSquare();
    frame->SetYTitle("Events / GeV");
    frame->SetXTitle("M_{#gamma#gamma} [GeV]");
    frame->Draw();
    
    TLatex text; text.SetNDC(); text.SetTextColor(1);
    text.DrawLatex(0.2, 0.81, currCateName);
    TH1F *histDH = new TH1F("histDH", "histDH", 1, 0, 1);
    TH1F *histSH = new TH1F("histSH", "histSH", 1, 0, 1);
    TH1F *histBkg = new TH1F("histBkg", "histBkg", 1, 0, 1);
    TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
    histDH->SetLineColor(6);
    histSH->SetLineColor(3);
    histBkg->SetLineColor(4);
    histSig->SetLineColor(2);
    TLegend leg(0.61, 0.63, 0.89, 0.77);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.SetBorderSize(0);
    leg.AddEntry(histDH, "Di-Higgs", "l");
    leg.AddEntry(histSH, "Single Higgs", "l");
    leg.AddEntry(histBkg, "Non-resonant", "l");
    leg.AddEntry(histSig, "Sum", "l");
    leg.Draw("SAME");
    can->Print(Form("%s/fitPlot_%s_%s_%s.eps", m_plotDir.Data(),
		    m_DHSignal.Data(), fitType.Data(), currCateName.Data()));
    */  
  }
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Check whether the specified map entry exists.
    @param mapKey - the key for which we are finding a value:
    @returns - true iff the categorization has been defined. 
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
   @param doSaveSnapshot - true iff you want to save snapshots in future fits.
*/
void DHTestStat::saveSnapshots(bool doSaveSnapshot) {
  m_doSaveSnapshot = doSaveSnapshot;
}

/**
   -----------------------------------------------------------------------------
   Set an output directory and enable plotting.
   @param directory - the output directory path.
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
void DHTestStat::setParams(TString paramName, double paramVal,
			   bool doSetConstant) {
  m_setParamConsts.push_back(doSetConstant);
  m_setParamNames.push_back(paramName);
  m_setParamVals.push_back(paramVal);
}
