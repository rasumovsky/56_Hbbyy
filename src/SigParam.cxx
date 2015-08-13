////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParam.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch  <-- Please use for reporting issues!                //
//  Date: 25/06/2015                                                          //
//                                                                            //
//  This class implements the resonance modeling for the ATLAS Hgamma group.  //
//                                                                            //
//  General notes:                                                            //
//                                                                            //
//    - Function names can be "DoubleCB" for double-sided Crystal Ball or     //
//      CBGA for Crystal Ball + Gaussian.                                     //
//                                                                            //
//    - Category indices should start at zero.                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "HGamTools/SigParam.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the SigParam class.
   @param signalType - The type of signal (ggH, VBF, WH, ZH, ttH, bbH, etc).
   @param directory - The directory for input and output files.
*/
SigParam::SigParam(TString signalType, TString directory) {
  std::cout << "\nSigParam::Initializing..." << std::endl;
  
  // Assign output directory:
  setSignalType(signalType);
  setDirectory(directory);
  
  // Create RooWorkspace and RooCategory to store all signal models:
  m_ws = new RooWorkspace("signalWS");
  m_cat = new RooCategory("signalCates","signalCates");
  
  // Import mass and weight variables, then make pointers:
  m_ws->factory(Form("m_yy[%f,%f]", 10.0, 5000.0));
  m_ws->factory(Form("wt[%f]", 1.0));
  m_ws->factory(Form("mResonance[%f,%f]", 10.0, 5000.0));
  m_yy = m_ws->var("m_yy");
  m_wt = m_ws->var("wt");
  m_mResonance = m_ws->var("mResonance");
  m_ws->factory(Form("expr::mRegularized('(@0-100.0)/100.0',{mResonance})"));
  
  // Define the data sets at each mass in each category for resonance fit:
  m_massCatePairs.clear();

  // Lists of mass resolution and mass scale systematics:
  m_listMRS = "";
  m_listMSS = "";
    
  // Set the default initial values and ranges for fit parameters:
  setParamState("a_muCBNom", "[-0.0,-2.0,2.0]");
  setParamState("b_muCBNom", "[-0.1,-0.5,0.5]");
  setParamState("c_muCBNom", "[-0.02,-0.5,0.5]");
  setParamState("a_sigmaCBNom", "[0.0,-10.0,10.0]");
  setParamState("b_sigmaCBNom", "[3.90,0.01,10.0]");
  setParamState("a_alphaCB", "[2.2,0.0,4.0]");
  setParamState("b_alphaCB", "[0.0,-0.1,0.1]");
  setParamState("nCB", "[5.0,0.1,10.0]");
  setParamState("a_muGANom", "[-0.0,-2.0,2.0]");
  setParamState("b_muGANom", "[-0.1,-0.5,0.5]");
  setParamState("c_muGANom", "[-0.02,-0.5,0.5]");
  setParamState("a_sigmaGANom", "[5.0,-1.0,20.0]");
  setParamState("b_sigmaGANom", "[1.8,0.1,4.0]");
  setParamState("fracCB", "[0.9,0.0,1.0]");
  setParamState("a_alphaCBLo", "[2.42,1.0,4.0]");
  setParamState("b_alphaCBLo", "[-483,-1000,0]");
  setParamState("c_alphaCBLo", "[380,100,500]");
  setParamState("nCBLo", "[9.0,0.1,20.0]");
  setParamState("a_alphaCBHi", "[2.2,0.0,4.0]");
  setParamState("b_alphaCBHi", "[0.0,-0.5,0.5]");
  setParamState("nCBHi", "[5.0,0.1,10.0]");
  
  // Some basic fit options:
  useCommonCBGAMean(false);
  
  // Set the plot format:
  doBinnedFit(false, 1);
  setLogYAxis(false);
  setPlotFormat(".eps");
  setRatioPlot(true, 0.0, 2.0);
  m_currFunction = "DoubleCB";
  
  std::cout << "SigParam: Successfully initialized!" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Adds the data from a RooDataSet object to the dataset for parameterization.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The event category index (starting at 0).
   @param dataSet - The RooDataSet object to import.
   @param observableName - The name of the observable.
*/
void SigParam::addDataSet(double resonanceMass, int cateIndex,
			  RooDataSet* dataSet, TString observableName) {
  
  // Loop over events stored in the input dataset:
  for (int i_d = 0; i_d < dataSet->numEntries(); i_d++) {
    // Get the RooRealVars and values in the current event:
    const RooArgSet *currArgs = (RooArgSet*)dataSet->get(i_d);
    
    double massValue = -999.0; double weightValue = -999.0;
    
    // Iterate through the RooArgSet, find mass and weight variables:
    TIterator *iterArgs = currArgs->createIterator();
    RooRealVar* currArg;
    while ((currArg = (RooRealVar*)iterArgs->Next())) {
      if (TString(currArg->GetName()).EqualTo(observableName)) {
	massValue = currArg->getVal();
      }
    }
    weightValue = dataSet->weight();
    // Then add the mass and weight values to the dataset:
    if (massValue > -990 && weightValue > -990) {
      addMassPoint(resonanceMass, cateIndex, massValue, weightValue);
    }
    // Exit if the variable names were incorrect:
    else {
      std::cout << "SigParam: No match for observable or weight." << std::endl;
      exit(0);
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Adds the data from a Root TTree object or an MxAOD to the dataset for
   signal parameterization. Note: the input 
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The event category index (starting at 0).
   @param dataTree - The TTree or MxAOD object to import.
   @param massBranchName - The name of the branch storing mass values.
   @param weightBranchName - The name of the branch storing event weights.
*/
void SigParam::addDataTree(double resonanceMass, int cateIndex,
			   TTree* dataTree, TString massBranchName,
			   TString weightBranchName) {
  // Get the mass and weight branches.
  double massValue;
  double weightValue;
  dataTree->SetBranchAddress(massBranchName, &massValue);
  dataTree->SetBranchAddress(weightBranchName, &weightValue);
  
  // Loop over the TTree.
  for (Long64_t i_t = 0; i_t < dataTree->GetEntries(); i_t++) {
    
    // Add current event values to the dataset used for parameterization.
    addMassPoint(resonanceMass, cateIndex, massValue, weightValue);
  }
}

/**
   -----------------------------------------------------------------------------
   Adds a parameter to the workspace for the fit.
   @param paramName - The name of the shape parameter of interest.
   @param cateIndex - The event category index (starting at 0).
*/
void SigParam::addFitParameter(TString paramName, int cateIndex) {
  TString valAndRange = getParamState(paramName);
  m_ws->factory(Form("%s_%sc%d%s", paramName.Data(), m_signalType.Data(),
		     cateIndex, valAndRange.Data()));
}

/**
   -----------------------------------------------------------------------------
   Adds a mass point to the dataset for fitting the signal diphoton resonance.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The event category index (starting at 0).
   @param diphotonMass - The reconstructed diphoton invariant mass.
   @param eventWeight - The weight of the event.
*/
void SigParam::addMassPoint(double resonanceMass, int cateIndex, 
			    double diphotonMass, double eventWeight) {

  TString currKey = getKey(resonanceMass, cateIndex);
  
  // Create new dataset if the corresponding one doesn't yet exist.
  if (!dataExists(resonanceMass, cateIndex)) {

    RooDataSet* newData 
      = new RooDataSet(Form("data_%s",currKey.Data()),
		       Form("data_%s",currKey.Data()),
		       RooArgSet(*m_yy,*m_wt), RooFit::WeightVar(*m_wt));
    m_ws->import(*newData);
    
    // Each dataset also corresponds to a unique category in the fit:
    m_cat->defineType(currKey);
    
    // Keep track of all resonance mass - category pairs for fitting.
    std::pair<double,int> currPair;
    currPair.first = resonanceMass;
    currPair.second = cateIndex;
    m_massCatePairs.push_back(currPair);
  }
  
  // Set the observable and weight values and then fill dataset:
  m_yy->setVal(diphotonMass);
  m_wt->setVal(eventWeight);
  RooDataSet* currData =(RooDataSet*)m_ws->data(Form("data_%s",currKey.Data()));
  currData->add(RooArgSet(*m_yy,*m_wt), eventWeight);
}

/**
   -----------------------------------------------------------------------------
   Add a single mass resolution systematic uncertainty to the signal shape.
   Note: the constraint terms must be defined separately. 
   @param nameMResSys - The name of the systematic nuisance parameter.
*/
void SigParam::addMResSystematic(TString nameMResSys) {

  // Name of nuisance parameter:
  TString atlasExpMRS = Form("atlas_expected_%s",nameMResSys.Data());
  
  // Check to see if the nuisance parameter is already in the workspace:
  if (!(bool)m_ws->obj(atlasExpMRS)) {
    m_ws->factory(Form("%s[1]",atlasExpMRS.Data()));
  }
  
  // Add to the list storing signal resolution nuisance parameters:
  if (!m_listMRS.Contains(nameMResSys)) {
    m_listMRS.Append(Form(",%s",atlasExpMRS.Data()));
  }
}

/**
   -----------------------------------------------------------------------------
   Add several mass resolution systematic uncertainties to the signal shape. 
   Note: the constraint terms must be defined separately. 
   @param namesMResSys - A list of the names of systematic nuisance parameters.
*/
void SigParam::addMResSystematics(std::vector<TString> namesMResSys) {
  for (std::vector<TString>::iterator sysIter = namesMResSys.begin();
       sysIter != namesMResSys.end(); sysIter++) {
    SigParam::addMResSystematic(*sysIter);
  }
}

/**
   -----------------------------------------------------------------------------
   Add a single mass scale systematic uncertainty to the signal shape.
   Note: the constraint terms must be defined separately. 
   @param nameMScaleSys - The name of the systematic nuisance parameter.
*/
void SigParam::addMScaleSystematic(TString nameMScaleSys) {
  
  // Name of nuisance parameter:
  TString atlasExpMSS = Form("atlas_expected_%s", nameMScaleSys.Data());
  
  // Check to see if the nuisance parameter is already in the workspace:
  if (!(bool)m_ws->obj(atlasExpMSS)) {
    m_ws->factory(Form("%s[1]",atlasExpMSS.Data()));
  }
  
  // Add to the list storing signal mass scale nuisance parameters:
  if (!m_listMSS.Contains(nameMScaleSys)) {
    m_listMSS.Append(Form(",%s",atlasExpMSS.Data()));
  }
}

/**
   -----------------------------------------------------------------------------
   Add several mass scale systematic uncertainties to the signal shape. 
   Note: the constraint terms must be defined separately. 
   @param namesMScaleSys - List of the names of systematic nuisance parameters.
*/
void SigParam::addMScaleSystematics(std::vector<TString> namesMScaleSys) {
  for (std::vector<TString>::iterator sysIter = namesMScaleSys.begin();
       sysIter != namesMScaleSys.end(); sysIter++) {
    SigParam::addMScaleSystematic(*sysIter);
  }
}

/**
   -----------------------------------------------------------------------------
   Adds a signal PDF from this class to a pre-existing workspace. 
   @param workspace - The pre-existing workspace.
   @param cateIndex - The index of the category of the desired PDF.
   @returns - True iff the PDF and yield parameter were imported.
*/
bool SigParam::addSigToWS(RooWorkspace *&workspace, int cateIndex) {
  std::cout << "SigParam: Adding parameterized " << m_signalType 
	    << " signal in category " << cateIndex
	    << " to pre-existing workspace." << std::endl;
  bool goodStatus = true;
  
  // Add the signal model to the workspace:
  RooAbsPdf* currSignal
    = m_ws->pdf(Form("sigPdf_%sc%d",m_signalType.Data(),cateIndex));
  if (currSignal && !workspace->pdf(currSignal->GetName())) {
    workspace->import(*currSignal);
  }
  else { 
    std::cout << "SigParam: Signal doesn't exist or is duplicate!" << std::endl;
    return false;
  }
  // Then explicitly add the parameters to the workspace:
  RooArgSet *currSet = currSignal->getVariables();
  TIterator *iterParam = currSet->createIterator();
  RooRealVar* currParam = NULL;
  while ((currParam = (RooRealVar*)iterParam->Next())) {
    workspace->import(*currParam);
  }
  
  // Import the yield formula and associated parameters:
  workspace->import(*m_ws->var(Form("yieldVar_a_%sc%d",
				    m_signalType.Data(), cateIndex)));
  workspace->import(*m_ws->var(Form("yieldVar_b_%sc%d",
				    m_signalType.Data(), cateIndex)));
  workspace->import(*m_ws->var(Form("yieldVar_c_%sc%d",
				    m_signalType.Data(), cateIndex)));
  workspace->import(*m_ws->var(Form("yieldVar_d_%sc%d",
				    m_signalType.Data(), cateIndex)));
  TString yieldName = Form("sigYield_%sc%d", m_signalType.Data(), cateIndex);
  if (m_ws->function(yieldName)) {
    if (!workspace->function(yieldName)) {
      workspace->import((*m_ws->function(yieldName)));
    }
    else {
      std::cout << "SigParam: yield parameter " << yieldName 
		<< " already found in workspace." << std::endl;
    }
  }
  else {
    std::cout << "SigParam: Error! Yield param was not created by this program."
	      << std::endl;
    goodStatus = false;
  }
  
  std::cout << "SigParam:: Finished adding parameterized signal." << std::endl;
  return goodStatus;
}

/**
   -----------------------------------------------------------------------------
   Adds a signal PDF from this class to a pre-existing workspace. 
   @param workspace - The pre-existing workspace.
   @param resonanceMass - The mass of the resonance.
   @param cateIndex - The index of the category of the desired PDF.
   @returns - True iff the PDF and yield parameter were imported.
*/
bool SigParam::addSigToWS(RooWorkspace *&workspace, double resonanceMass,
			  int cateIndex) {
  std::cout << "SigParam: Adding individual " << m_signalType
	    << " signal in category " << cateIndex << " with mass " 
	    << resonanceMass << " to pre-existing workspace." << std::endl;
  bool goodStatus = true;

  TString currKey = getKey(resonanceMass, cateIndex);
  // Add the signal model to the workspace, if not already contained:
  RooAbsPdf* currSignal
    = m_ws->pdf(Form("sigPdf_%s%s",m_signalType.Data(), currKey.Data()));
  if (currSignal && !workspace->pdf(currSignal->GetName())) {
    workspace->import(*currSignal);
  }
  else { 
    std::cout << "SigParam: Signal doesn't exist or is duplicate!" << std::endl;
    return false;
  }
  
  // Then explicitly add the parameters to the workspace:
  RooArgSet *currSet = currSignal->getVariables();
  TIterator *iterParam = currSet->createIterator();
  RooRealVar* currParam = NULL;
  while ((currParam = (RooRealVar*)iterParam->Next())) {
    if (!workspace->var(currParam->GetName())) {
      workspace->import(*currParam);
    }
    else {
      std::cout << "SigParam: parameter " << currParam->GetName()
		<< " already found in workspace." << std::endl;
    }
  }
  
  // add yield factor to workspace:
  TString yieldName = Form("sigYield_%s%s",m_signalType.Data(), currKey.Data());
  if (m_ws->var(yieldName)) {
    if (!workspace->var(yieldName)) {
      workspace->import((*m_ws->var(yieldName)));
    }
    else {
      std::cout << "SigParam: yield parameter " << yieldName 
		<< " already found in workspace." << std::endl;
    }
  }
  else {
    std::cout << "SigParam: Error! Yield param was not created by this program."
	      << std::endl;
    goodStatus = false;
  }
  std::cout << "SigParam:: Finished adding individual signal." << std::endl;
  return goodStatus;
}

/**
   Creates a single or combined binned dataset. 
   @param unbinnedDataName - The name of the unbinned data in the workspace.
   @param cateIndex - The index of the category of the desired PDF.
*/
void SigParam::binTheData(TString unbinnedDataName, int cateIndex) {
  if (!unbinnedDataName.Contains(Form("data_c%d",cateIndex))) {
    TString currBinnedName = unbinnedDataName; 
    currBinnedName.ReplaceAll("data","dataBinned");
    binSingleDataSet(unbinnedDataName, currBinnedName);
  }
  else {
    // Create a dataset map to store the RooDataSet objects:
    std::map<std::string,RooDataSet*> currDataMap; currDataMap.clear();
    
    // Create a list of the corresponding unbinned datasets, and loop over it:
    std::vector<double> currMassPoints = massPointsForCategory(cateIndex); 
    for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
      TString currKey = getKey(currMassPoints[i_m], cateIndex);

      // bin each one, add to map
      TString currUnbinnedName = Form("data_%s", currKey.Data());
      TString currBinnedName = Form("dataBinned_%s", currKey.Data());
      binSingleDataSet(currUnbinnedName, currBinnedName);
      currDataMap[((std::string)currKey)]
	= (RooDataSet*)m_ws->data(currBinnedName);
    }
    
    // Create the combined binned dataset:
    RooArgSet *args = new RooArgSet();
    args->add(*(m_ws->var("m_yy")));
    args->add(*(m_ws->var("wt")));
    RooDataSet *dataBinned = new RooDataSet(Form("dataBinned_c%d", cateIndex),
					    Form("dataBinned_c%d", cateIndex),
					    *args, RooFit::Index(*m_cat),
					    RooFit::Import(currDataMap),
					    RooFit::WeightVar(*m_wt));
    m_ws->import(*dataBinned);
  }
}

/**
   Create a binned dataset from an unbinned dataset.
   @param unbinnedName - The name of the unbinned data in the workspace.
   @param unbinnedName - The name of the binned data to create.
*/
void SigParam::binSingleDataSet(TString unbinnedName, TString binnedName) {
  std::cout << "SigParam: Rebin sample " << unbinnedName << std::endl;
  
  RooDataSet *unbinnedData = (RooDataSet*)m_ws->data(unbinnedName);
  
  // Create a histogram to automatically bin the points:
  int nBins = (int)(m_nBinsPerGeV*(m_yy->getMax() - m_yy->getMin()));
  TH1 *hist = unbinnedData
    ->createHistogram(Form("hist_%s",unbinnedName.Data()), *m_yy, 
		      RooFit::Binning(nBins, m_yy->getMin(), m_yy->getMax()));
  
  // Create a dataset to fill with binned entries:
  RooDataSet* binnedData = new RooDataSet(binnedName, binnedName,
					  RooArgSet(*m_yy,*m_wt), 
					  RooFit::WeightVar(*m_wt));
  for (int i_b = 1; i_b <= hist->GetNbinsX(); i_b++) {
    m_yy->setVal(hist->GetXaxis()->GetBinCenter(i_b));
    m_wt->setVal(hist->GetBinContent(i_b));
    binnedData->add(RooArgSet(*m_yy,*m_wt), hist->GetBinContent(i_b));
  }
  
  // Import binned dataset to the workspace:
  m_ws->import(*binnedData);
}

/**
   -----------------------------------------------------------------------------
   Get a list of categories corresponding to a single mass point.
   @param resonanceMass - The mass value.
   @returns - A vector of category indices.
*/
std::vector<int> SigParam::categoriesForMass(double resonanceMass) {
  // Create a list of mass points in this category.
  std::vector<int> currCategories; currCategories.clear();
  for (int i_p = 0; i_p < (int)m_massCatePairs.size(); i_p++) {
    if (equalMasses((m_massCatePairs[i_p]).first, resonanceMass)) {
      currCategories.push_back((m_massCatePairs[i_p]).second);
    }
  }
  std::sort(currCategories.begin(), currCategories.end());
  return currCategories;
}

/**
   -----------------------------------------------------------------------------
   Check if the dataset being requested has been instantiated (exists in map).
   @param massIndex - The index of the signal mass.
   @param cateIndex - The index of the category.
   @returns - True iff the dataset has been defined.
*/
bool SigParam::dataExists(double resonanceMass, int cateIndex) {
  for (int i_p = 0; i_p < (int)m_massCatePairs.size(); i_p++) {
    if (equalMasses((m_massCatePairs[i_p]).first, resonanceMass) &&
	m_massCatePairs[i_p].second == cateIndex) {
      return true;
    }
  }
  return false;
}

/**
   Option to do a binned fit. Sets the private variable m_binned.
   @param binned - True iff the fit should be binned.
   @param nBinsPerGeV - The number of bins per GeV for the binned data.
*/
void SigParam::doBinnedFit(bool binned, double nBinsPerGeV) {
  m_binned = binned;
  m_nBinsPerGeV = nBinsPerGeV;
  std::cout << "SigParam: Binned bool = " << m_binned << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Check if two doubles are equal.
   @param massValue1 - The first mass value to compare.
   @param massValue2 - The second mass value to compare.
   @returns - True iff the masses are equal within 0.001 GeV.
*/
bool SigParam::equalMasses(double massValue1, double massValue2) {
  return (fabs(massValue1 - massValue2) <= 0.001);// mass precision in GeV
}

/**
   -----------------------------------------------------------------------------
   Perform a single or simultaneous fit.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category.
   @returns - The RooFitResult, which gives fit status.
*/
RooFitResult* SigParam::fitResult(double resonanceMass, int cateIndex) {
  std::cout << "SigParam: Preparing to fit resonance" << std::endl;
  
  // Clock the fit:
  clock_t time;
  time = clock();
  
  // Define the PDF and dataset names:
  TString sigName = (resonanceMass < 0.0) ?
    Form("sigPdfTmp_%sc%d", m_signalType.Data(), cateIndex) :
    Form("sigPdf_%s%s", m_signalType.Data(), 
	 (getKey(resonanceMass,cateIndex)).Data());
  TString dataName = (resonanceMass < 0.0) ? Form("data_c%d",cateIndex) :
    Form("data_%s",(getKey(resonanceMass,cateIndex)).Data());
  
  // Change the dataset name if the fit is binned:
  if (m_binned) {
    TString binnedDataName = dataName;
    binnedDataName.ReplaceAll("data", "dataBinned");
    
    // Also check that binned data have been created, otherwise create:
    if (m_ws->data(dataName) && !m_ws->data(binnedDataName)) {
      binTheData(dataName, cateIndex);
    }
    
    dataName = binnedDataName;
  }

  RooAbsPdf *currSignal = m_ws->pdf(sigName);
  RooAbsData *currData = m_ws->data(dataName);
  
  // Free, fit, fix, then save the nuisance parameter values:
  SigParam::setParamsConstant(currSignal, false);
  RooFitResult *result = currSignal->fitTo(*currData, RooFit::PrintLevel(0),
					   RooFit::SumW2Error(kTRUE),
					   RooFit::Save(true));
  /*
  RooFitResult *result = currSignal->fitTo(*currData, RooFit::PrintLevel(0),
					   RooFit::SumW2Error(kFALSE),
					   RooFit::Save(true));
  */
  SigParam::setParamsConstant(currSignal, true);
  std::cout << "SigParam: Fit procedure concluded." << std::endl;
  
  // Clock the fit:
  time = clock() - time;
  printf("\tFit required %d clock cycles (%f seconds).\n\n",
	 (int)time, ((float)time/CLOCKS_PER_SEC));
  return result;
}

/**
   -----------------------------------------------------------------------------
   Perform a simultaneous fit across multiple masses.
   @param cateIndex - The index of the category.
   @returns - The RooFitResult, which gives fit status.
*/
RooFitResult* SigParam::fitResult(int cateIndex) {
  return SigParam::fitResult(-999.9, cateIndex);
}

/**
   -----------------------------------------------------------------------------
   Retrieve a key string for the given mass and category index.
   @param resonanceMass - The floating value of the signal mass.
   @param cateIndex - The index of the category.
   @returns - A key string specific to the mass and category.
*/
TString SigParam::getKey(double resonanceMass, int cateIndex) {
  TString key = Form("m%d_c%d", massDoubleToInt(resonanceMass), cateIndex);
  return key;
}

/**
   -----------------------------------------------------------------------------
   Get the number of categories contained in the datasets for fitting. Note:
   it is possible that there are different numbers of categories defined for 
   different mass points. This is up to the user to sort out. 
   @returns - The total number of categories for the parameterization.
*/
int SigParam::getNCategories() {
  m_nCategories = 0;
  // Loop over mass-category pairs, find highest category index. 
  for (int i_p = 0; i_p < (int)m_massCatePairs.size(); i_p++) {
    // The number of categories is equal to the highest index + 1:
    if ((m_massCatePairs[i_p]).second > m_nCategories + 1) {
      m_nCategories = (m_massCatePairs[i_p]).second + 1;
    }
  }
  return m_nCategories;
}

/**
   -----------------------------------------------------------------------------
   Get the value of the fit error for a particular parameter of the signal PDF. 
   @param paramName - The name of the shape parameter of interest.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category.
   @returns - The value of the specified signal parameter. 
*/
double SigParam::getParameterError(TString paramName, double resonanceMass,
				   int cateIndex) {
  if (paramName.Contains("nCB") || paramName.Contains("fracCB")) {
    return SigParam::getParameterError(paramName, cateIndex);
  }
  else {
    RooRealVar *var
      = m_ws->var(Form("%s_%s%s", paramName.Data(), m_signalType.Data(), 
		       (getKey(resonanceMass,cateIndex)).Data()));
    if (!var) {
      std::cout << "SigParam: requested parameter not found: " 
		<< paramName << std::endl;
      return 0.0;
    }
    else {
      return var->getError();
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Get the value of the fit error for a particular parameter of the signal PDF. 
   @param paramName - The name of the shape parameter of interest.
   @param cateIndex - The index of the category.
   @returns - The value of the specified signal parameter. 
*/
double SigParam::getParameterError(TString paramName, int cateIndex) {
  RooRealVar *var = m_ws->var(Form("%s_%sc%d", paramName.Data(),
				   m_signalType.Data(), cateIndex));
  if (!var) {
    std::cout << "SigParam: requested parameter not found: "
	      << paramName << std::endl;
    return 0.0;
  }
  else {
    return var->getError();
  }
}

/**
   -----------------------------------------------------------------------------
   Get the value of a particular parameter of the signal PDF. 
   @param paramName - The name of the shape parameter of interest.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category.
   @returns - The value of the specified signal parameter. 
*/
double SigParam::getParameterValue(TString paramName, double resonanceMass, 
				   int cateIndex) {
  if (paramName.Contains("nCB") || paramName.Contains("fracCB")) {
    return SigParam::getParameterValue(paramName, cateIndex);
  }
  else {
    RooRealVar *var = m_ws->var(Form("%s_%s%s", paramName.Data(), 
				     m_signalType.Data(),
				     (getKey(resonanceMass,cateIndex)).Data()));
    if (!var) {
      std::cout << "SigParam: requested parameter not found: "
		<< paramName << std::endl;
      return 0.0;
    }
    else {
      return var->getVal();
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Get the value of a particular parameter of the signal PDF. 
   @param paramName - The name of the shape parameter of interest.
   @param cateIndex - The index of the category.
   @returns - The value of the specified signal parameter. 
*/
double SigParam::getParameterValue(TString paramName, int cateIndex) {
  RooRealVar *var = m_ws->var(Form("%s_%sc%d", paramName.Data(),
				   m_signalType.Data(), cateIndex));
  if (!var) {
    std::cout << "SigParam: requested parameter not found: "
	      << paramName << std::endl;
    return 0.0;
  }
  else {
    return var->getVal();
  }
}

/**
   -----------------------------------------------------------------------------
   Get the initial value [a] and range [b,c] of a parameter: "[a,b,c]"
   @param paramName - The name of the parameter of interest.
   @returns - The initial value and range of the parameter. 
*/
TString SigParam::getParamState(TString paramName) {
  return m_paramState[paramName];
}

/**
   -----------------------------------------------------------------------------
   Retrieves the resonance parameterized as a function of mResonance.
   @param cateIndex - The index of the category for the PDF.
   @returns - A pointer to the signal PDF.
*/
RooAbsPdf* SigParam::getResonance(int cateIndex) {
  std::cout << "SigParam: Get parameterized shape in category = " 
	    << cateIndex << std::endl;
  TString pdfName = Form("sigPdf_%sc%d",m_signalType.Data(), cateIndex);
  RooAbsPdf* pdf = m_ws->pdf(pdfName);
  std::cout << "SigParam: Returning parameterized pdf " << pdfName << std::endl;
  return pdf;
}

/**
   -----------------------------------------------------------------------------
   Get the resonance shape for a single category and mass.
   @param resonanceMass - The truth mass of the resonance
   @param cateIndex - The index of the category for which we want the PDF.
   @returns - A pointer to the signal PDF.
*/
RooAbsPdf* SigParam::getSingleResonance(double resonanceMass, int cateIndex) {
  std::cout << "SigParam: Get parameterized shape in category = " 
	    << cateIndex << " and mass = " << resonanceMass << std::endl;
  TString pdfName = Form("sigPdf_%s%s",m_signalType.Data(),
			 (getKey(resonanceMass,cateIndex)).Data());
  RooAbsPdf* pdf = m_ws->pdf(pdfName);
  std::cout << "SigParam: Returning parameterized pdf " << pdfName << std::endl;
  return pdf;
}

/**
   -----------------------------------------------------------------------------
   Get a list of parameters associated with the PDF in the given category.
   @param resonanceMass - The truth mass of the resonance
   @param cateIndex - The index of the category for which we want the PDF.
   @returns - A vector of parameter names.
*/
std::vector<TString> SigParam::getVariableNames(double resonanceMass, 
						int cateIndex) {
  std::cout << "SigParam: Get PDF variables in category = " << cateIndex
	    << " and mass = " << resonanceMass << std::endl;
  std::vector<TString> result; result.clear();
  TString pdfName = Form("sigPdf_%s%s",m_signalType.Data(),
			 (getKey(resonanceMass,cateIndex)).Data());
  RooArgSet* currSet = (*m_ws->pdf(pdfName)).getVariables();
  TIterator *iterArg = currSet->createIterator();
  RooRealVar *currArg = NULL;
  while ((currArg = (RooRealVar*)iterArg->Next())) {
    TString currArgName = currArg->GetName();
    currArgName = currArgName.ReplaceAll(getKey(resonanceMass,cateIndex), "");
    currArgName = currArgName.ReplaceAll(Form("_%s",m_signalType.Data()), "");
    currArgName = currArgName.ReplaceAll(Form("c%d",cateIndex), "");
    if (!currArgName.EqualTo(m_yy->GetName())) result.push_back(currArgName);
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Get the workspace that stores the SigParam class data.
   @returns - A pointer to the workspace storing all datasets and PDFs.
*/
RooWorkspace* SigParam::getWorkspace() {
  return m_ws;
}

/**
   -----------------------------------------------------------------------------
   Get the signal yield for a particular mass in a particular category.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category for which we want the PDF.
   @returns - The signal yield for the specified mass in the given category.
*/
double SigParam::getYieldInCategory(double resonanceMass, int cateIndex) {
  std::cout << "SigParam: Get yield in category = " 
	    << cateIndex << " at mass = " << resonanceMass << std::endl;
  
  TString currKey = getKey(resonanceMass,cateIndex);
  // First check if parameterized yield is available:
  if (m_ws->function(Form("sigYield_%sc%d", m_signalType.Data(), cateIndex))) {
    (*m_ws->var("mResonance")).setVal(resonanceMass);
    return (*m_ws->function(Form("sigYield_%sc%d", m_signalType.Data(),
				 cateIndex))).getVal();
  }
  // Then try to get individual yield:
  else if(m_ws->var(Form("sigYield_%s%s",m_signalType.Data(),currKey.Data()))) {
    return (*m_ws->var(Form("sigYield_%s%s",m_signalType.Data(),
			    currKey.Data()))).getVal();
  }
  // Then try to get directly from dataset normalization:
  else if ((m_ws->data(Form("data_%s",currKey.Data())))) {
    return (*m_ws->data(Form("data_%s",currKey.Data()))).sumEntries();
  }
  // Or return error message:
  else {
    std::cout << "SigParam: requested yield not found." << std::endl;
    return 0.0;
  }
}

/**
   -----------------------------------------------------------------------------
   Get the signal yield for a particular resonance mass in all categories.
   @param resonanceMass - The truth mass of the resonance.
   @returns - The signal yield in all categories for the specified mass.
*/
double SigParam::getYieldTotal(double resonanceMass) {
  std::cout << "SigParam: Get total yield at mass = " 
	    << resonanceMass << std::endl;
  
  // Create a list of categories for this mass point:
  std::vector<int> currCategories = SigParam::categoriesForMass(resonanceMass);
  // Loop through names of datasets, add components:
  double sum = 0.0;
  for (int i_c = 0; i_c < (int)currCategories.size(); i_c++) {
    sum += getYieldInCategory(resonanceMass, currCategories[i_c]);
  }
  return sum;
}

/**
   -----------------------------------------------------------------------------
   Load the signal parameterization from file.
   @param directory - Name of the directory housing the input workspace.
   @returns - True iff the file is successfully loaded.
*/
bool SigParam::loadParameterization(TString directory, TString signalType){
  std::cout << "SigParam: Load parameterization from" << directory
	    << std::endl;
  setDirectory(directory);
  bool parameterizationExists = false;
  TFile inputFile(Form("%s/res_%sworkspace.root", m_directory.Data(),m_signalType.Data()));
  if (inputFile.IsOpen()) {
    m_ws = (RooWorkspace*)inputFile.Get("signalWS");
    if (m_ws) parameterizationExists = true;
    m_yy = m_ws->var("m_yy");
    m_wt = m_ws->var("wt");
    m_mResonance = m_ws->var("mResonance");
    setSignalType(signalType);
    std::cout << "SigParam: Successfully loaded from file!" << std::endl;
  }
  return parameterizationExists;
}

/**
   -----------------------------------------------------------------------------
   Parameterize the resonance shape in all categories.
   @param function - The functional form of the resonance.
   @returns - True iff. all fits converge.
*/
bool SigParam::makeAllParameterizations(TString function) {
  std::cout << "SigParam: Engage full signal parameterization!" << std::endl;
  if (!function.EqualTo("DoubleCB") && !function.EqualTo("CBGA")) {
    std::cout << "SigParam: Error! Improper function " << function << std::endl;
    return false;
  }
    
  bool result = true;
  // Define models in each category independently:
  for (int i_c = 0; i_c < getNCategories(); i_c++) {
    if (!makeCategoryParameterization(i_c, function)) result = false;
  }
  
  std::cout << "SigParam: Completed full signal parameterization" << std::endl;
  return result;
}

/**
   -----------------------------------------------------------------------------
   Parameterize the resonance shape as a function of mass for a single category.
   @param cateIndex - The index of the category to fit.
   @param function - The functional form of the resonance.
   @returns - True iff. all fits converge.
*/
bool SigParam::makeCategoryParameterization(int cateIndex, TString function) {
  std::cout << "SigParam: parameterizing category " << cateIndex << std::endl;
  if (!function.EqualTo("DoubleCB") && !function.EqualTo("CBGA")) {
    std::cout << "SigParam: Error! Improper function " << function << std::endl;
    return false;
  }
  
  // Create a list of mass points in this category.
  std::vector<double> currMassPoints = massPointsForCategory(cateIndex);
    
  // Create a simultaneous PDF just for this category:
  RooSimultaneous *currSim
    = new RooSimultaneous(Form("sigPdfTmp_%sc%d",m_signalType.Data(),cateIndex),
			  Form("sigPdfTmp_%sc%d",m_signalType.Data(),cateIndex),
			  *m_cat);
  
  // If no mass points, cannot fit signal, silly!
  if (currMassPoints.size() == 0) {
    std::cout << "SigParam: No masspoints for cate. " << cateIndex << std::endl;
    return false;
  }
  // If only one mass point, no need for parameterization:
  else if (currMassPoints.size() == 1) {
    std::cout << "SigParam: 1 mass point -> no parameterization." << std::endl;
    return makeSingleResonance(currMassPoints[0], cateIndex, function);
  }
  // If more than 1 mass points, parameterize variables:
  else {
    addFitParameter("a_muCBNom", cateIndex);
    addFitParameter("b_muCBNom", cateIndex);
    addFitParameter("c_muCBNom", cateIndex);
    addFitParameter("a_sigmaCBNom", cateIndex);
    addFitParameter("b_sigmaCBNom", cateIndex);
    if (function.EqualTo("CBGA")) {
      addFitParameter("a_alphaCB", cateIndex);
      addFitParameter("b_alphaCB", cateIndex);
      addFitParameter("nCB", cateIndex);
      if (!m_sameCBGAMean) {
	addFitParameter("a_muGANom", cateIndex);
	addFitParameter("b_muGANom", cateIndex);
	addFitParameter("c_muGANom", cateIndex);
      }
      addFitParameter("a_sigmaGANom", cateIndex);
      addFitParameter("b_sigmaGANom", cateIndex);
      addFitParameter("fracCB", cateIndex);
    }
    else if (function.EqualTo("DoubleCB")) {
      addFitParameter("a_alphaCBLo", cateIndex);
      addFitParameter("b_alphaCBLo", cateIndex);
      addFitParameter("c_alphaCBLo", cateIndex);
      addFitParameter("nCBLo", cateIndex);
      addFitParameter("a_alphaCBHi", cateIndex);
      addFitParameter("b_alphaCBHi", cateIndex);
      addFitParameter("nCBHi", cateIndex);
    }
    
    // Loop over mass points, define resonance model in each:
    for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
      TString currKey = getKey(currMassPoints[i_m],cateIndex);
      double mRegVal = regularizedMass(currMassPoints[i_m]);
      double mResVal = currMassPoints[i_m];
      // Define the RooFormulaVars which control mH parameterization:
      m_ws->factory(Form("expr::muCBNom_%s%s('@0+@1*%f+@2*%f*%f+%f',{a_muCBNom_%sc%d,b_muCBNom_%sc%d,c_muCBNom_%sc%d})", m_signalType.Data(), currKey.Data(), mRegVal, mRegVal, mRegVal, mResVal, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
      m_ws->factory(Form("expr::sigmaCBNom_%s%s('@0+@1*%f',{a_sigmaCBNom_%sc%d,b_sigmaCBNom_%sc%d})", m_signalType.Data(), currKey.Data(), mRegVal, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
      // Crystal Ball + Gaussian-specific parameters:
      if (function.EqualTo("CBGA")) {
	m_ws->factory(Form("expr::alphaCB_%s%s('@0+@1*%f',{a_alphaCB_%sc%d,b_alphaCB_%sc%d})", m_signalType.Data(), currKey.Data(), mRegVal, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
	if (!m_sameCBGAMean) {
	  m_ws->factory(Form("expr::muGANom_%s%s('@0+@1*%f+@2*%f*%f+%f',{a_muGANom_%sc%d,b_muGANom_%sc%d,c_muGANom_%sc%d})", m_signalType.Data(), currKey.Data(), mRegVal, mRegVal, mRegVal, mResVal, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
	}
	m_ws->factory(Form("expr::sigmaGANom_%s%s('@0+@1*%f',{a_sigmaGANom_%sc%d,b_sigmaGANom_%sc%d})", m_signalType.Data(), currKey.Data(), mRegVal, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
      }
      // Double Crystal Ball-specific parameters:
      else if (function.EqualTo("DoubleCB")) {
	m_ws->factory(Form("expr::alphaCBLo_%s%s('@0+@1/(%f+@2)',{a_alphaCBLo_%sc%d,b_alphaCBLo_%sc%d,c_alphaCBLo_%sc%d})", m_signalType.Data(), currKey.Data(), mRegVal, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
	m_ws->factory(Form("expr::alphaCBHi_%s%s('@0+@1*%f',{a_alphaCBHi_%sc%d,b_alphaCBHi_%sc%d})", m_signalType.Data(), currKey.Data(), mRegVal, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
      }
      
      // Create, individual resonance shapes, add to simultaneous PDF:
      resonanceCreator(currMassPoints[i_m], cateIndex, function);
      currSim->addPdf(*m_ws->pdf(Form("sigPdf_%s%s",m_signalType.Data(),
				      currKey.Data())), currKey);
    }
    
    // Import simultaneous PDF into workspace:
    m_ws->import(*currSim);
    
    // Get the RooDataSets:
    std::map<std::string,RooDataSet*> currDataMap; currDataMap.clear();
    for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
      // maybe add if statement here...
      TString currKey = getKey(currMassPoints[i_m], cateIndex);
      currDataMap[((std::string)currKey)]
	= (RooDataSet*)m_ws->data(Form("data_%s",currKey.Data()));
    }
    
    RooArgSet *args = new RooArgSet();
    args->add(*(m_ws->var("m_yy")));
    args->add(*(m_ws->var("wt")));
    RooDataSet *obsData = new RooDataSet(Form("data_c%d",cateIndex),
					 Form("data_c%d",cateIndex), *args,
					 RooFit::Index(*m_cat),
					 RooFit::Import(currDataMap),
					 RooFit::WeightVar(*m_wt));
    m_ws->import(*obsData);
    
    // Then fit simultaneous PDF to combined dataset:
    RooFitResult *result = fitResult(cateIndex);
    
    // Then construct parametric resonance (function of mResonance):
    m_ws->factory(Form("expr::muCBNom_%sc%d('@0+@1*@3+@2*@3*@3+@4',{a_muCBNom_%sc%d,b_muCBNom_%sc%d,c_muCBNom_%sc%d,mRegularized,mResonance})", m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
    m_ws->factory(Form("expr::sigmaCBNom_%sc%d('@0+@1*@2',{a_sigmaCBNom_%sc%d,b_sigmaCBNom_%sc%d,mRegularized})", m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
    // Crystal Ball + Gaussian-specific parameters:
    if (function.EqualTo("CBGA")) {
      m_ws->factory(Form("expr::alphaCB_%sc%d('@0+@1*@2',{a_alphaCB_%sc%d,b_alphaCB_%sc%d,mRegularized})", m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
      if (!m_sameCBGAMean) {
	m_ws->factory(Form("expr::muGANom_%sc%d('@0+@1*@3+@2*@3*@3+@4',{a_muGANom_%sc%d,b_muGANom_%sc%d,c_muGANom_%sc%d,mRegularized,mResonance})", m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
      }
      m_ws->factory(Form("expr::sigmaGANom_%sc%d('@0+@1*@2',{a_sigmaGANom_%sc%d,b_sigmaGANom_%sc%d,mRegularized})", m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
    }
    // Double Crystal Ball-specific parameters:
    else if (function.EqualTo("DoubleCB")) {
      m_ws->factory(Form("expr::alphaCBLo_%sc%d('@0+@1/(@3+@2)',{a_alphaCBLo_%sc%d,b_alphaCBLo_%sc%d,c_alphaCBLo_%sc%d,mRegularized})", m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
      m_ws->factory(Form("expr::alphaCBHi_%sc%d('@0+@1*@2',{a_alphaCBHi_%sc%d,b_alphaCBHi_%sc%d,mRegularized})", m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
    }
    // Create the parameterized resonance shape (function of mResonance):
    resonanceCreator(-999, cateIndex, Form("%s_Parameterized",function.Data()));
    
    // Get the yield parameterization:
    makeYieldParameterization(cateIndex);

    // Return the fit status:
    return (result->status() == 0);
  }
}

/**
   -----------------------------------------------------------------------------
   Create the resonance for a single mass point and category. 
   @param resonanceMass - The truth mass of the resonance
   @param cateIndex - The index of the category to fit.
   @param function - The functional form of the resonance.
   @returns - True iff. all fits converge.
*/
bool SigParam::makeSingleResonance(double resonanceMass, int cateIndex,
				   TString function) {
  std::cout << "SigParam: parameterizing single point in category " << cateIndex
	    << " and at mass " << resonanceMass << std::endl;
  if (!function.EqualTo("DoubleCB") && !function.EqualTo("CBGA")) {
    std::cout << "SigParam: Error! Improper function " << function << std::endl;
    return false;
  }
  
  // Before calling resonanceCreator, need to define dependent variables.
  TString currKey = getKey(resonanceMass, cateIndex);
  m_ws->factory(Form("muCBNom_%s%s[%f,%f,%f]", m_signalType.Data(), 
		     currKey.Data(), resonanceMass, 0.9*resonanceMass, 
		     1.1*resonanceMass));
  m_ws->factory(Form("sigmaCBNom_%s%s[2.0,0.01,200.0]", m_signalType.Data(), 
		     currKey.Data()));
  
  if (function.EqualTo("CBGA")) {
    m_ws->factory(Form("alphaCB_%s%s[1.5,1.0,2.5]", m_signalType.Data(), 
		       currKey.Data()));
    if (!m_ws->var(Form("nCB_%sc%d",m_signalType.Data(), cateIndex))) {
      //m_ws->factory(Form("nCB_%sc%d[9.0,0.01,20.0]", m_signalType.Data(), 
      //		 cateIndex));
      m_ws->factory(Form("nCB_%sc%d[9.0]", m_signalType.Data(), cateIndex));
    }
    //m_ws->factory(Form("sigmaGANom_%s%s[2.0,0.01,40.0]", m_signalType.Data(),
    //		       currKey.Data()));
    if (!m_sameCBGAMean) {
      m_ws->factory(Form("muGANom_%s%s[%f,%f,%f]", m_signalType.Data(), 
			 currKey.Data(), resonanceMass, 0.8*resonanceMass, 
			 1.2*resonanceMass));
    }
    m_ws->factory(Form("sigmaGANom_%s%s[10.0,0.01,80.0]", m_signalType.Data(),
		       currKey.Data()));
    //m_ws->factory(Form("fracCB_%sc%d[0.9,0.001,0.999]", m_signalType.Data(), 
    //		       cateIndex));
    m_ws->factory(Form("fracCB_%sc%d[0.5,0.001,0.999]", m_signalType.Data(), 
		       cateIndex));
  }
  else if (function.EqualTo("DoubleCB")) {
    m_ws->factory(Form("alphaCBLo_%s%s[1.5,1.0,2.5]", m_signalType.Data(),
		       currKey.Data()));
    if (!m_ws->var(Form("nCBLo_%sc%d", m_signalType.Data(), cateIndex))) {
      m_ws->factory(Form("nCBLo_%sc%d[17.0,0.01,30.0]", m_signalType.Data(),
			 cateIndex));
    }
    m_ws->factory(Form("alphaCBHi_%s%s[2.2,1.0,3.0]", m_signalType.Data(),
		       currKey.Data()));
    if (!m_ws->var(Form("nCBHi_%sc%d", m_signalType.Data(), cateIndex))) {
      m_ws->factory(Form("nCBHi_%sc%d[5.2,0.01,10.0]", m_signalType.Data(), 
			 cateIndex));
    }
  }
  
  // Create the single resonance PDF:
  resonanceCreator(resonanceMass, cateIndex, function);
  
  // Then fit single PDF to single dataset:
  RooFitResult *result = fitResult(resonanceMass,cateIndex);
  // Return the fit status:
  return (result->status() == 0);
}

/**
   -----------------------------------------------------------------------------
   Parameterizes the signal yields in a category as a function of mH.
   @param cateIndex - The index of the category to fit.
*/
void SigParam::makeYieldParameterization(int cateIndex) {
  std::cout << "SigParam: Parameterizing the signal yield." << std::endl;
  
  // Create arrays to store fit data:
  int nResPoints = 0;
  double mResValues[100] = {0.0};
  double yieldValues[100] = {0.0};
  
  // Retrieve dataset yield for each mass point that was imported:
  std::vector<double> currMassPoints = massPointsForCategory(cateIndex);
  for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
    TString dataName
      = Form("data_%s",(getKey(currMassPoints[i_m],cateIndex)).Data());
    if ((m_ws->data(dataName))) {
      mResValues[nResPoints] = regularizedMass(currMassPoints[i_m]);
      yieldValues[nResPoints] = (*m_ws->data(dataName)).sumEntries();
      nResPoints++;
    }
  }
  
  // Use TF1 and TGraph to fit the yield:
  m_yieldFunc[cateIndex] = new TF1("yieldFunc", "pol3", 0.0, 0.5);
  m_yieldGraph[cateIndex] = new TGraph(nResPoints, mResValues, yieldValues);
  m_yieldGraph[cateIndex]->Fit(m_yieldFunc[cateIndex]);
  
  // Create the yield parameters:
  m_ws->factory(Form("yieldVar_a_%sc%d[%f]", m_signalType.Data(), cateIndex,
		     m_yieldFunc[cateIndex]->GetParameter(0)));
  m_ws->factory(Form("yieldVar_b_%sc%d[%f]", m_signalType.Data(), cateIndex,
		     m_yieldFunc[cateIndex]->GetParameter(1)));
  m_ws->factory(Form("yieldVar_c_%sc%d[%f]", m_signalType.Data(), cateIndex, 
		     m_yieldFunc[cateIndex]->GetParameter(2)));
  m_ws->factory(Form("yieldVar_d_%sc%d[%f]", m_signalType.Data(), cateIndex,
		     m_yieldFunc[cateIndex]->GetParameter(3)));
  
  // Then create a yield RooFormulaVar.
  m_ws->factory(Form("expr::sigYield_%sc%d('@0+@1*@4+@2*@4*@4+@3*@4*@4*@4',{yieldVar_a_%sc%d,yieldVar_b_%sc%d,yieldVar_c_%sc%d,yieldVar_d_%sc%d,mRegularized})", m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex, m_signalType.Data(), cateIndex));
}

/**
   -----------------------------------------------------------------------------
   Convert the resonance mass integer to a floating value.
   @param massInteger - An integer value representing the mass.
   @returns - The floating value of the mass in GeV.
*/
double SigParam::massIntToDouble(int massInteger) {
  return ((double)massInteger) / 1000.0;
}

/**
   -----------------------------------------------------------------------------
   Convert the resonance mass value to an integer representation. 
   @param resonanceMass - The value of the mass.
   @returns - The integer representation of the mass.
*/
int SigParam::massDoubleToInt(double resonanceMass) {
  return (int)(resonanceMass * 1000.0);
}

/**
   -----------------------------------------------------------------------------
   Get a list of mass points corresponding to a single category.
   @param cateIndex - The index of the category.
   @returns - A vector of mass values.
*/
std::vector<double> SigParam::massPointsForCategory(int cateIndex) {
  // Create a list of mass points in this category.
  std::vector<double> currMassPoints; currMassPoints.clear();
  for (int i_p = 0; i_p < (int)m_massCatePairs.size(); i_p++) {
    if ((m_massCatePairs[i_p]).second == cateIndex) {
      currMassPoints.push_back((m_massCatePairs[i_p]).first);
    }
  }
  std::sort(currMassPoints.begin(), currMassPoints.end());
  return currMassPoints;
}

/**
   -----------------------------------------------------------------------------
   Create a ratio plot (or subtraction plot, for the moment...)
   @param data - The RooAbsData set for comparison.
   @param pdf - The PDF for comparison.
   @param xMin - The minimum value of the observable range.
   @param xMax - The maximum value of the observable range.
   @param xBins - The number of bins for the observable.
   @returns - A TGraphErrors to plot.
*/
TGraphErrors* SigParam::plotSubtraction(RooAbsData *data, RooAbsPdf *pdf, 
					double xMin, double xMax, double xBins){
  double minOrigin = m_yy->getMin();
  double maxOrigin = m_yy->getMax();
  double nEvents = data->sumEntries();
  
  m_yy->setRange("fullRange", xMin, xMax);
  TH1F *originHist
    = (TH1F*)data->createHistogram("dataSub", *m_yy, 				
  				   RooFit::Binning(xBins, xMin, xMax));
  TGraphErrors *result = new TGraphErrors();
  double increment = (xMax - xMin) / ((double)xBins);
  
  RooAbsReal* intTot
    = (RooAbsReal*)pdf->createIntegral(RooArgSet(*m_yy),
				       RooFit::NormSet(*m_yy), 
				       RooFit::Range("fullRange"));
  double valTot = intTot->getVal();
  
  int pointIndex = 0;
  for (double i_m = xMin; i_m < xMax; i_m += increment) {
    m_yy->setRange(Form("range%2.2f",i_m), i_m, (i_m+increment));
    RooAbsReal* intCurr
      = (RooAbsReal*)pdf->createIntegral(RooArgSet(*m_yy), 
					 RooFit::NormSet(*m_yy), 
					 RooFit::Range(Form("range%2.2f",i_m)));
    double valCurr = intCurr->getVal();
    
    double currMass = i_m + (0.5*increment);
    double currPdfWeight = nEvents * (valCurr / valTot);
    TString varName = m_yy->GetName();
    double currDataWeight = data->sumEntries(Form("%s>%f&&%s<%f",varName.Data(),
						  i_m,varName.Data(),
						  (i_m+increment)));
    double currWeight = currDataWeight - currPdfWeight;
    result->SetPoint(pointIndex, currMass, currWeight);
  
    double currError = originHist->GetBinError(pointIndex+1);
    result->SetPointError(pointIndex, 0.0, currError);
    pointIndex++;
  }
  m_yy->setMin(minOrigin);
  m_yy->setMax(maxOrigin);
  return result;
}

/**
   -----------------------------------------------------------------------------
   Create a ratio plot (or subtraction plot, for the moment...)
   @param data - The RooAbsData set for comparison.
   @param pdf - The PDF for comparison.
   @param xMin - The minimum value of the observable range.
   @param xMax - The maximum value of the observable range.
   @param xBins - The number of bins for the observable.
   @returns - A TGraphErrors to plot.
*/
TGraphErrors* SigParam::plotDivision(RooAbsData *data, RooAbsPdf *pdf, 
				     double xMin, double xMax, double xBins){
  double minOrigin = m_yy->getMin();
  double maxOrigin = m_yy->getMax();
  double nEvents = data->sumEntries();
  
  m_yy->setRange("fullRange", xMin, xMax);
  TH1F *originHist
    = (TH1F*)data->createHistogram("dataSub", *m_yy, 				
  				   RooFit::Binning(xBins, xMin, xMax));
  TGraphErrors *result = new TGraphErrors();
  double increment = (xMax - xMin) / ((double)xBins);
  
  RooAbsReal* intTot
    = (RooAbsReal*)pdf->createIntegral(RooArgSet(*m_yy),
				       RooFit::NormSet(*m_yy), 
				       RooFit::Range("fullRange"));
  double valTot = intTot->getVal();
  
  int pointIndex = 0;
  for (double i_m = xMin; i_m < xMax; i_m += increment) {
    m_yy->setRange(Form("range%2.2f",i_m), i_m, (i_m+increment));
    RooAbsReal* intCurr
      = (RooAbsReal*)pdf->createIntegral(RooArgSet(*m_yy), 
					 RooFit::NormSet(*m_yy), 
					 RooFit::Range(Form("range%2.2f",i_m)));
    double valCurr = intCurr->getVal();
    
    double currMass = i_m + (0.5*increment);
    double currPdfWeight = nEvents * (valCurr / valTot);
    TString varName = m_yy->GetName();
    double currDataWeight = data->sumEntries(Form("%s>%f&&%s<%f",varName.Data(),
						  i_m,varName.Data(),
						  (i_m+increment)));
    double currWeight = currDataWeight / currPdfWeight;
    result->SetPoint(pointIndex, currMass, currWeight);
  
    double currError = originHist->GetBinError(pointIndex+1) / currPdfWeight;
    result->SetPointError(pointIndex, 0.0, currError);
    pointIndex++;
  }
  m_yy->setMin(minOrigin);
  m_yy->setMax(maxOrigin);
  return result;
}

/**
   -----------------------------------------------------------------------------
   Plot a resonance PDF for all masses defined for one category.
   @param cateIndex - The index of the category.
*/
void SigParam::plotCategoryResonances(int cateIndex) {
  std::cout << "SigParam: Plot resonances in cate. " << cateIndex << std::endl;
  
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
    
  // Get a list of the resonance masses:
  std::vector<double> currMassPoints = massPointsForCategory(cateIndex);
  int xMin = currMassPoints[0] - 10;
  int xMax = currMassPoints[currMassPoints.size()-1] + 10;
  int xBins = (int)(m_nBinsPerGeV) * (xMax - xMin);
    
  RooPlot* frame = m_yy->frame(RooFit::Bins(xBins), RooFit::Range(xMin, xMax));
  frame->SetYTitle("Events/0.5 GeV");
  frame->SetXTitle("Mass [GeV]");
    
  // Loop over mass points, drawing data and PDF for each.
  for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
    
    pad1->cd();
    
    RooAbsData *currData = NULL;
    RooAbsPdf *currPdf = NULL;
    
    TString currKey = getKey(currMassPoints[i_m], cateIndex);
    if ((m_ws->data(Form("data_%s",currKey.Data())))) {
      currData = (m_ws->data(Form("data_%s",currKey.Data())));
      currData->plotOn(frame);
    }
    else {
      std::cout << "SigParam: data for plotting undefined." << std::endl;
      return;
    }
    
    if ((m_ws->pdf(Form("sigPdf_%sc%d",m_signalType.Data(),cateIndex)))) {
      (*m_ws->var("mResonance")).setVal(currMassPoints[i_m]);
      currPdf = (m_ws->pdf(Form("sigPdf_%sc%d",m_signalType.Data(),cateIndex)));
      currPdf->plotOn(frame, RooFit::LineColor(4));
    }
    else if ((m_ws->pdf(Form("sigPdf_%s%s",
			      m_signalType.Data(),currKey.Data())))) {
      currPdf
	= (m_ws->pdf(Form("sigPdf_%s%s",m_signalType.Data(),currKey.Data())));
      currPdf->plotOn(frame, RooFit::LineColor(4));
    }
    else {
      std::cout << "SigParam: resonance for plotting undefined." << std::endl;
      return;
    }
    
    if (i_m == 0) frame->Draw();
    else frame->Draw("SAME");
    
    TLatex text; text.SetNDC(); text.SetTextColor(1);
    text.DrawLatex(0.72, 0.88, Form("category %d", cateIndex));
    
    // Switch to sub-plot:
    pad2->cd();
    TGraphErrors* subData = (m_doRatioPlot) ? 
      plotDivision(currData,currPdf,xMin,xMax,xBins) :
      plotSubtraction(currData,currPdf,xMin,xMax,xBins);
    
    if (i_m == 0) {
      subData->GetXaxis()->SetTitle("Mass [GeV]");
      if (m_doRatioPlot) {
	subData->GetYaxis()->SetTitle("Data / Fit");
	subData->GetYaxis()->SetRangeUser(m_ratioMin, m_ratioMax);
      }
      else subData->GetYaxis()->SetTitle("Data - Fit");
      subData->GetXaxis()->SetTitleOffset(0.95);
      subData->GetYaxis()->SetTitleOffset(0.7);
      subData->GetXaxis()->SetTitleSize(0.1);
      subData->GetYaxis()->SetTitleSize(0.1);
      subData->GetXaxis()->SetLabelSize(0.1);
      subData->GetYaxis()->SetLabelSize(0.1);
      subData->GetYaxis()->SetNdivisions(4);
      subData->SetMarkerColor(1);
      subData->GetXaxis()->SetRangeUser(xMin, xMax);
      subData->Draw("AEP");
      
      TLine *line = new TLine();
      line->SetLineStyle(1);
      line->SetLineWidth(2);
      line->SetLineColor(kRed);
      if (m_doRatioPlot) {
	line->DrawLine(xMin, 1.0, xMax, 1.0);
	line->SetLineWidth(1);
	line->SetLineStyle(2);
	line->DrawLine(xMin,((1.0+m_ratioMin)/2.0),xMax,((1.0+m_ratioMin)/2.0));
	line->DrawLine(xMin,((1.0+m_ratioMax)/2.0),xMax,((1.0+m_ratioMax)/2.0));
      }
      else line->DrawLine(xMin, 0.0, xMax, 0.0);
      subData->Draw("EPSAME");
    }
    subData->Draw("EPSAME");
  }
  can->Print(Form("%s/plot_paramResonance_c%d%s", 
		  m_directory.Data(), cateIndex, m_fileFormat.Data()));
}

/**
   -----------------------------------------------------------------------------
   Plot a resonance PDF for one value of the resonance mass in one category.
   @param resonanceMass - The mass value in GeV.
   @param cateIndex - The index of the category.
*/
void SigParam::plotSingleResonance(double resonanceMass, int cateIndex) {
  std::cout << "SigParam: Plotting resonance at mass " << resonanceMass 
	    << " in category " << cateIndex << std::endl;
  
  TCanvas *can = new TCanvas("can","can",800,800);
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
  pad1->cd();
  
  double rMin = 0.75*resonanceMass;
  double rMax = 1.25*resonanceMass;
  int rBins = (int)(m_nBinsPerGeV*(rMax - rMin));
  
  RooPlot* frame = m_yy->frame(RooFit::Bins(rBins), RooFit::Range(rMin,rMax));
  frame->SetYTitle("Events/0.5 GeV");
  frame->SetXTitle("Mass [GeV]");
  
  RooAbsData *currData = NULL;
  RooAbsPdf *currPdf = NULL;
  
  TString currKey = getKey(resonanceMass,cateIndex);
  if ((m_ws->data(Form("data_%s",currKey.Data())))) {
    currData = (m_ws->data(Form("data_%s",currKey.Data())));
    (*m_ws->data(Form("data_%s",currKey.Data()))).plotOn(frame);
  }
  else {
    std::cout << "SigParam: data for plotting undefined." << std::endl;
    return;
  }
  
  // First check to see if parameterized shape exists:
  if ((m_ws->pdf(Form("sigPdf_%sc%d",m_signalType.Data(),cateIndex)))) {
    (*m_ws->var("mResonance")).setVal(resonanceMass);    
    currPdf = (m_ws->pdf(Form("sigPdf_%sc%d",m_signalType.Data(),cateIndex)));
    if (m_currFunction.EqualTo("CBGA")) {
      currPdf->plotOn(frame, RooFit::Components((*m_ws->pdf(Form("pdfGA_%sc%d",m_signalType.Data(), cateIndex)))), RooFit::LineColor(3), RooFit::LineStyle(2));
      currPdf->plotOn(frame, RooFit::Components((*m_ws->pdf(Form("pdfCB_%sc%d",m_signalType.Data(), cateIndex)))), RooFit::LineColor(4), RooFit::LineStyle(2));
    }
    currPdf->plotOn(frame, RooFit::LineColor(2));
  }
  else if((m_ws->pdf(Form("sigPdf_%s%s",m_signalType.Data(),currKey.Data())))){
    currPdf
      = (m_ws->pdf(Form("sigPdf_%s%s",m_signalType.Data(),currKey.Data())));
    if (m_currFunction.EqualTo("CBGA")) {
      currPdf->plotOn(frame, RooFit::Components((*m_ws->pdf(Form("pdfGA_%s%s",m_signalType.Data(), currKey.Data())))), RooFit::LineColor(3), RooFit::LineStyle(2));
      currPdf->plotOn(frame, RooFit::Components((*m_ws->pdf(Form("pdfCB_%s%s",m_signalType.Data(), currKey.Data())))), RooFit::LineColor(4), RooFit::LineStyle(2));
    }
    currPdf->plotOn(frame, RooFit::LineColor(2));
  }
  else {
    std::cout << "SigParam: resonance for plotting undefined." << std::endl;
    return;
  }
  frame->Draw();
  if (m_useLogYAxis) {
    gPad->SetLogy();
    frame->GetYaxis()->SetRangeUser(0.00001 *(*m_ws->data(Form("data_%s",currKey.Data()))).sumEntries(), (*m_ws->data(Form("data_%s",currKey.Data()))).sumEntries());
  }
  TLatex text; text.SetNDC(); text.SetTextColor(1);
  text.DrawLatex(0.2, 0.88, Form("category %d", cateIndex));
  pad2->cd();
  
  TGraphErrors* subData = (m_doRatioPlot) ?
    plotDivision(currData, currPdf, rMin, rMax, rBins) : 
    plotSubtraction(currData, currPdf, rMin, rMax, rBins);
  
  subData->GetXaxis()->SetTitle("Mass [GeV]");
  if (m_doRatioPlot) {
    subData->GetYaxis()->SetTitle("Data / Fit");
    subData->GetYaxis()->SetRangeUser(m_ratioMin, m_ratioMax);
  }
  else subData->GetYaxis()->SetTitle("Data - Fit");
  subData->GetXaxis()->SetTitleOffset(0.95);
  subData->GetYaxis()->SetTitleOffset(0.7);
  subData->GetXaxis()->SetTitleSize(0.1);
  subData->GetYaxis()->SetTitleSize(0.1);
  subData->GetXaxis()->SetLabelSize(0.1);
  subData->GetYaxis()->SetLabelSize(0.1);
  subData->GetYaxis()->SetNdivisions(4);
  subData->GetXaxis()->SetRangeUser(rMin, rMax);
  subData->Draw("AEP");
  
  TLine *line = new TLine();
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  if (m_doRatioPlot) {
    line->DrawLine(rMin, 1.0, rMax, 1.0); 
    line->SetLineWidth(1);
    line->SetLineStyle(2);
    line->DrawLine(rMin,((1.0+m_ratioMin)/2.0),rMax,((1.0+m_ratioMin)/2.0));
    line->DrawLine(rMin,((1.0+m_ratioMax)/2.0),rMax,((1.0+m_ratioMax)/2.0));
  }
  else line->DrawLine(rMin, 0.0, rMax, 0.0);
  subData->Draw("EPSAME");
  can->Print(Form("%s/plot_singleRes_m%2.2f_c%d%s", m_directory.Data(),
		  resonanceMass, cateIndex, m_fileFormat.Data()));
}

/**
   -----------------------------------------------------------------------------
   Graph the signal yields as function of resonance mass in the given category.
   @param cateIndex - The index of the category for yield plotting. 
*/
void SigParam::plotYields(int cateIndex) {
  std::cout << "SigParam: Plotting yields in category " << cateIndex
	    << std::endl;
  
  // Get a list of the mass points in the category:
  //TCanvas *canY = new TCanvas("can", "can", 800, 600);
  //canY->cd();
  
  TCanvas *can = new TCanvas("can","can",800,800);
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
  pad1->cd();
  
  m_yieldFunc[cateIndex]->SetLineColor(kBlue);
  m_yieldGraph[cateIndex]->Draw("AP");
  m_yieldFunc[cateIndex]->Draw("LSAME");
  TLatex text; text.SetNDC(); text.SetTextColor(1);
  text.DrawLatex(0.4, 0.88, Form("%s signal, category %d",
				 m_signalType.Data(), cateIndex));
  pad2->cd();
  
  TGraphErrors *gRatio = new TGraphErrors(Form("ratio_c%d",cateIndex),
					  Form("ratio_c%d",cateIndex));
  for (int i_p = 0; i_p < m_yieldGraph[cateIndex]->GetN(); i_p++) {
    double currX = 0.0; double currY = 0.0;
    m_yieldGraph[cateIndex]->GetPoint(i_p, currX, currY);
    gRatio->SetPoint(i_p, currX, (currY/m_yieldFunc[cateIndex]->Eval(currX)));
  }
  gRatio->GetXaxis()->SetTitle("Mass [GeV]");
  gRatio->GetYaxis()->SetTitle("Data - Fit");
  gRatio->GetXaxis()->SetTitleOffset(0.95);
  gRatio->GetYaxis()->SetTitleOffset(0.7);
  gRatio->GetXaxis()->SetTitleSize(0.1);
  gRatio->GetYaxis()->SetTitleSize(0.1);
  gRatio->GetXaxis()->SetLabelSize(0.1);
  gRatio->GetYaxis()->SetLabelSize(0.1);
  gRatio->GetYaxis()->SetNdivisions(4);
  gRatio->GetXaxis()
    ->SetRangeUser(m_yieldGraph[cateIndex]->GetXaxis()->GetXmin(),
		   m_yieldGraph[cateIndex]->GetXaxis()->GetXmin());
  gRatio->Draw("AEP");
  gRatio->Draw("AP");
  
  TLine *line = new TLine();
  line->SetLineStyle(1);
  line->SetLineWidth(1);
  line->SetLineColor(kBlue);
  line->DrawLine(gRatio->GetXaxis()->GetXmin(), 1,
		 gRatio->GetXaxis()->GetXmax(), 1);
  
  can->Print(Form("%s/plot_paramYield_%sc%d%s", m_directory.Data(),
		  m_signalType.Data(), cateIndex, m_fileFormat.Data()));
}

/**
   -----------------------------------------------------------------------------
   Create a regularized mass variable that helps minimizers converge.
   @param resonanceMass - The mass value in GeV.
   @returns - The regularized mass mR = (m-100)/100;
*/
double SigParam::regularizedMass(double resonanceMass) {
  return ((resonanceMass - 100.0) / 100.0);
}

/**
   -----------------------------------------------------------------------------
   Create the resonance shape corresponding to a single mass point in a single
   analysis category. The signal will be stored in the class workspace. 
   @param resonanceMass - the truth mass of the resonance.
   @param cateIndex - the index of the category to fit.
   @param function - the functional form of the resonance.
*/
void SigParam::resonanceCreator(double resonanceMass, int cateIndex, 
				TString function) {
  std::cout << "SigParam: Adding resonance to workspace" << "\tresonanceMass: "
	    << resonanceMass << "\tcategory: " << cateIndex << std::endl;
  m_currFunction = function;
  
  // Check that the dataset has been defined and is not empty.
  if (!dataExists(resonanceMass, cateIndex) &&
      !function.Contains("Parameterized")) {
    std::cout << "SigParam: Cannot fit resonance, no dataset." << std::endl;
    exit(0);
  }
  
  TString currKey = getKey(resonanceMass, cateIndex);
  TString suffix = Form("%s%s", m_signalType.Data(), currKey.Data());
  if (function.Contains("Parameterized")) {
    suffix = Form("%sc%d", m_signalType.Data(), cateIndex);
  }
  
  // Define the Crystal Ball + Gaussian shape:
  if (function.Contains("CBGA")) {
    // Cystal Ball component:
    m_ws->factory(Form("RooCBShape::pdfCB_%s(m_yy, prod::muCB_%s(muCBNom_%s%s), prod::sigmaCB_%s(sigmaCBNom_%s%s), alphaCB_%s, nCB_%sc%d)", suffix.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data(), suffix.Data(), m_signalType.Data(), cateIndex));
    // Gaussian component:
    if (m_sameCBGAMean) {
      m_ws->factory(Form("RooGaussian::pdfGA_%s(m_yy, prod::muGA_%s(muCBNom_%s%s), prod::sigmaGA_%s(sigmaGANom_%s%s))", suffix.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
    }
    else {
      m_ws->factory(Form("RooGaussian::pdfGA_%s(m_yy, prod::muGA_%s(muGANom_%s%s), prod::sigmaGA_%s(sigmaGANom_%s%s))", suffix.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
    }
    
    // Crystal Ball + Gaussian:
    m_ws->factory(Form("SUM::sigPdf_%s(fracCB_%sc%d*pdfCB_%s,pdfGA_%s)", suffix.Data(), m_signalType.Data(), cateIndex, suffix.Data(), suffix.Data()));
  }
  
  // Define double-sided Crystal Ball shape:
  else if (function.Contains("DoubleCB")) {
    m_ws->factory(Form("HggTwoSidedCBPdf::sigPdf_%s(m_yy, prod::muCB_%s(muCBNom_%s%s), prod::sigmaCB_%s(sigmaCBNom_%s%s), alphaCBLo_%s, nCBLo_%sc%d, alphaCBHi_%s, nCBHi_%sc%d)", suffix.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data(), suffix.Data(), m_signalType.Data(), cateIndex, suffix.Data(), m_signalType.Data(), cateIndex));
  }
  
  // Define the yield:
  
  if (!function.Contains("Parameterized")) {
    double yieldValue
      = (*m_ws->data(Form("data_%s",currKey.Data()))).sumEntries();
    m_ws->factory(Form("sigYield_%s%s[%f]",m_signalType.Data(),currKey.Data(),yieldValue));
  }
  
  std::cout << "SigParam: Resonance successfully added." << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Save parameterization workspace, list of parameter values, and signal yields.
*/
void SigParam::saveAll() {
  saveParameterization();
  saveParameterList();
  saveYieldList();
}

/**
   -----------------------------------------------------------------------------
   Save the workspace containing the parameterization data to file.
*/
void SigParam::saveParameterization() {
  m_ws->writeToFile(Form("%s/res_%sworkspace.root", m_directory.Data(),
			 m_signalType.Data()));
}

/**
   -----------------------------------------------------------------------------
   Save a list of parameter names and values.
*/
void SigParam::saveParameterList() {
  std::ofstream outFile(Form("%s/resonance paramList.txt",
			m_directory.Data()));
  RooArgSet args = m_ws->allVars();
  TIterator *iterArgs = args.createIterator();
  RooRealVar *currIter = NULL;
  while ((currIter = (RooRealVar*)iterArgs->Next())) {
    outFile << currIter->GetName() << " " << currIter->getVal() << std::endl;
  }
  outFile.close();
}

/**
   -----------------------------------------------------------------------------
   Save a list of signal yields in all categories and at all mass points. 
*/
void SigParam::saveYieldList() {
  std::ofstream outFile(Form("%s/resonance_yieldList.txt", m_directory.Data()));
  // Loop over datasets:
  for (int i_d = 0; i_d < (int)m_massCatePairs.size(); i_d++) {
    double currMass = (m_massCatePairs[i_d]).first;
    int currCate = (m_massCatePairs[i_d]).second;
    outFile << currMass << " " << currCate << " "
	    << getYieldInCategory(currMass, currCate) << std::endl;
  }
  outFile.close();
}

/**
   -----------------------------------------------------------------------------
   Set the output directory for files.
   @param directory - The new input/output directory for files.
*/
void SigParam::setDirectory(TString directory) {
  TString simpSigName = m_signalType;
  simpSigName.Remove(simpSigName.Length()-1);
  if (directory.EqualTo("")) {
    if (simpSigName.EqualTo("")) m_directory = "";
    else m_directory = simpSigName;
  }
  else {
    if (simpSigName.EqualTo("")) m_directory = directory;
    else m_directory = Form("%s/%s", directory.Data(), simpSigName.Data());
  }
  
  // Create the output directory if it doesn't exist:
  system(Form("mkdir -vp %s", m_directory.Data()));
  std::cout << "SigParam: I/O directory set to " << directory << std::endl;
}

/**
   Use a logarithmic Y axis for the plots.
   @param useLogYAxis - True iff plots should have log axis.
*/
void SigParam::setLogYAxis(bool useLogYAxis) {
  m_useLogYAxis = useLogYAxis;
}

/**
   -----------------------------------------------------------------------------
   Set the initial value and range for a fit parameter.
   @param paramName - The name of the shape parameter of interest.
   @param valueAndRange - A string with the initial value and range: "[a,b,c]".
*/
void SigParam::setParamState(TString paramName, TString valueAndRange) {
  m_paramState[paramName] = valueAndRange;
}

/**
   -----------------------------------------------------------------------------
   Set the format for plots. Examples are ".eps", ".pdf", ".png", etc. Standard
   ROOT file extensions.
*/
void SigParam::setPlotFormat(TString fileFormat) {
  m_fileFormat = fileFormat;
}

/**
   Set the sub-plot option to a ratio plot instead of a subtraction plot.
   @param doRatioPlot - True iff you want to do a ratio plot.
   @param ratioMin - The minimum ratio plot y-axis range.
   @param ratioMax - The maximum ratio plot y-axis range.
*/
void SigParam::setRatioPlot(bool doRatioPlot, double ratioMin, double ratioMax){
  m_doRatioPlot = doRatioPlot;
  m_ratioMin = ratioMin;
  m_ratioMax = ratioMax;
}

/**
   -----------------------------------------------------------------------------
   Make the parameters of a PDF free or fixed.
   @param pdf - The PDF containing the parameters to be freed/fixed.
   @param isConstant - True iff setting the parameters constant.
*/
void SigParam::setParamsConstant(RooAbsPdf* pdf, bool isConstant) {
  RooArgSet *currArgs = pdf->getVariables();
  TIterator *iterArgs = currArgs->createIterator();
  RooRealVar* currIter = NULL;
  while ((currIter = (RooRealVar*)iterArgs->Next())) {
    currIter->setConstant(isConstant);
  }
}

/**
   -----------------------------------------------------------------------------
   Set the signal type to avoid collisions when using many production modes.
   @param signalType - The type of signal. 
*/
void SigParam::setSignalType(TString signalType) {
  if (signalType == "") {
    m_signalType = "";
  }
  else {
    m_signalType = Form("%s_", signalType.Data());
  }
  std::cout << "SigParam: Signal type set to " << signalType << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Tell the program whether or not to use the same value for the Crystal Ball
   and Gaussian means.
   @param sameCBGAMean - True iff the same mean value should be used. 
*/
void SigParam::useCommonCBGAMean(bool sameCBGAMean) {
  m_sameCBGAMean = sameCBGAMean;
}
