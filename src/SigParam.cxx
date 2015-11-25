////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParam.cxx                                                        //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch  <-- Please use for reporting issues!                //
//  Date: 08/11/2015                                                          //
//                                                                            //
//  This class implements the resonance modeling for the ATLAS Hgamma group.  //
//                                                                            //
//  General notes:                                                            //
//                                                                            //
//    - Function names can be:                                                //
//       * "DoubleCB" for double-sided Crystal Ball                           //
//       * "CBGA" for Crystal Ball + Gaussian                                 //
//       * "GAx3" for 3 Gaussians                                             //
//       * "BifurGA" for bifurcated Gaussian                                  //
//       * "Landau"                                                           //
//       * "CBPlusVoigt"                                                      //
//       * "Voigt"                                                            //
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
  m_ws->factory("wt[1.0]");
  m_ws->factory("mResonance[10,10000]");
  m_ws->factory("expr::mRegularized('(@0-100.0)/100.0',{mResonance})");
  m_verbose = true;
  m_nCategories = 0;
  
  // Define the data sets at each mass in each category for resonance fit:
  m_massCatePairs.clear();
  m_cateNames.clear();
  
  // Lists of mass resolution and mass scale systematics:
  m_listMRS = "";
  m_listMSS = "";
  
  // Set the default initial values and ranges for fit parameters:
  m_paramState.clear();
  m_varParameterization.clear();
  
  // Some basic fit options:
  useCommonCBGAMean(false);
  
  // Fit result information:
  m_currChi2 = 0.0;
  m_currNLL = 0.0;
  m_currExtendVal = 0.0;
  m_generatedDataNorm = 0.0;
  
  // Set the plot format:
  doBinnedFit(false, 1);
  setLogYAxis(false);
  setPlotFormat(".eps");
  setRatioPlot(true, 0.0, 2.0);
  m_currFunction = "DoubleCB";
  RooRealVar nameVar("functionName", m_currFunction, 0);
  m_ws->import(nameVar);
  
  // Default values for the parameterization functions and parameters are set
  // below. These will be overwritten by any calls to setVarParameterizat() or
  // setParamState() by the user after initializing the SigParam class.
  
  //--------------------------------------------------------------------------//
  // Default parameterization functions for variables:
  // Common parameters:
  setVarParameterization("muCBNom", "@0+@1*obs+@2*obs*obs+mRes");
  setVarParameterization("sigmaCBNom", "@0+@1*obs");
  // Crystal Ball + Gaussian:
  setVarParameterization("alphaCB", "@0+@1*obs");
  setVarParameterization("muGANom", "@0+@1*obs+@2*obs*obs+mRes");
  setVarParameterization("sigmaGANom", "@0+@1*obs");
  // Double-sided Crystal Ball:
  setVarParameterization("alphaCBLo", "@0+@1/(obs+@2)");
  //setVarParameterization("alphaCBHi", "@0+@1*obs");
  setVarParameterization("alphaCBHi", "@0+@1/(obs+@2)");
  // Triple Gaussian:
  setVarParameterization("muGA1Nom", "@0+@1*obs+mRes");
  setVarParameterization("muGA2Nom", "@0+@1*obs+mRes");
  setVarParameterization("muGA3Nom", "@0+@1*obs+mRes");
  setVarParameterization("sigmaGA1Nom", "@0+@1*obs");
  setVarParameterization("sigmaGA2Nom", "@0+@1*obs");
  setVarParameterization("sigmaGA3Nom", "@0+@1*obs");
  // Bifurcated Gaussian:
  setVarParameterization("sigmaGALowNom", "@0+@1*obs");
  setVarParameterization("sigmaGAHiNom", "@0+@1*obs");
  // Landau:
  setVarParameterization("muLANom", "@0+@1*obs+@2*obs*obs+mRes");
  setVarParameterization("sigmaLANom", "@0+@1*obs");
  // Voigt:
  setVarParameterization("muVoigtNom", "@0+@1*obs+@2*obs*obs+mRes");
  setVarParameterization("widthVoigtNom", "@0+@1*obs");
  setVarParameterization("sigmaVoigtNom", "@0+@1*obs"); 
  
  //--------------------------------------------------------------------------//
  // Default parameter values below are for parameterized fit over mResonance.
  // Common parameters:
  setParamState("a_muCBNom", "[-0.0,-2.0,2.0]");
  setParamState("b_muCBNom", "[-0.1,-0.5,0.5]");
  setParamState("c_muCBNom", "[-0.02,-0.5,0.5]");
  setParamState("a_sigmaCBNom", "[0.0,-10.0,10.0]");
  setParamState("b_sigmaCBNom", "[3.90,0.01,10.0]");
  // Crystal Ball + Gaussian:
  setParamState("a_alphaCB", "[2.2,0.0,4.0]");
  setParamState("b_alphaCB", "[0.0,-0.1,0.1]");
  setParamState("nCB", "[5.0,0.1,10.0]");
  setParamState("a_muGANom", "[-0.0,-2.0,2.0]");
  setParamState("b_muGANom", "[-0.1,-0.5,0.5]");
  setParamState("c_muGANom", "[-0.02,-0.5,0.5]");
  setParamState("a_sigmaGANom", "[5.0,-1.0,20.0]");
  setParamState("b_sigmaGANom", "[1.8,0.1,4.0]");
  setParamState("fracCB", "[0.9,0.01,0.99]");
  // Double-sided Crystal Ball:
  setParamState("a_alphaCBLo", "[2.42,1.0,4.0]");
  setParamState("b_alphaCBLo", "[-483,-1000,0]");
  setParamState("c_alphaCBLo", "[380,100,500]");
  setParamState("nCBLo", "[9.0,0.1,20.0]");
  setParamState("a_alphaCBHi", "[2.2,0.0,4.0]");
  setParamState("b_alphaCBHi", "[0.0,-0.5,0.5]");
  setParamState("c_alphaCBHi", "[0.0,-2.0,1.0]");
  setParamState("nCBHi", "[5.0,0.1,10.0]");
  // Triple Gaussian:
  setParamState("a_muGA1Nom", "[-0.01,-2.0,2.0]");
  setParamState("b_muGA1Nom", "[-0.1,-0.5,0.5]");
  setParamState("a_muGA2Nom", "[0.01,-2.0,2.0]");
  setParamState("b_muGA2Nom", "[-0.1,-0.5,0.5]");
  setParamState("a_muGA3Nom", "[0.0,-2.0,2.0]");
  setParamState("b_muGA3Nom", "[-0.1,-0.5,0.5]");
  setParamState("a_sigmaGA1Nom", "[0.0,-10.0,10.0]");
  setParamState("b_sigmaGA1Nom", "[3.90,0.01,10.0]");
  setParamState("a_sigmaGA2Nom", "[0.0,-10.0,10.0]");
  setParamState("b_sigmaGA2Nom", "[3.90,0.01,10.0]");
  setParamState("a_sigmaGA3Nom", "[0.0,-10.0,10.0]");
  setParamState("b_sigmaGA3Nom", "[3.90,0.01,10.0]");
  setParamState("fracGA1", "[0.32,0.01,0.99]");
  setParamState("fracGA2", "[0.33,0.01,0.99]");
  //Bifurcated Gaussian:
  setParamState("a_sigmaGALowNom", "[0.0,-10.0,10.0]");
  setParamState("b_sigmaGALowNom", "[3.90,0.01,10.0]");
  setParamState("a_sigmaGAHiNom", "[0.0,-10.0,10.0]");
  setParamState("b_sigmaGAHiNom", "[3.90,0.01,10.0]");
  // Landau:
  setParamState("a_muLANom", "[-0.0,-2.0,2.0]");
  setParamState("b_muLANom", "[-0.1,-0.5,0.5]");
  setParamState("c_muLANom", "[-0.02,-0.5,0.5]");
  setParamState("a_sigmaLANom", "[0.0,-10.0,10.0]");
  setParamState("b_sigmaLANom", "[3.90,0.01,10.0]");
  // Voigt:
   setParamState("a_muVoigtNom", "[-0.0,-2.0,2.0]");
  setParamState("b_muVoigtNom", "[-0.1,-0.5,0.5]");
  setParamState("c_muVoigtNom", "[-0.02,-0.5,0.5]");
  setParamState("a_widthVoigtNom", "[5.0,-1.0,20.0]");
  setParamState("b_widthVoigtNom", "[1.8,0.1,4.0]");
  setParamState("a_sigmaVoigtNom", "[5.0,-1.0,20.0]");
  setParamState("b_sigmaVoigtNom", "[1.8,0.1,4.0]");
  
  //--------------------------------------------------------------------------//
  // Parameter values below are for fits to individual mass points:
  // Common parameters:
  setParamState("muCBNom", "[125,0.1,10000]");
  setParamState("sigmaCBNom", "[1.7,0.01,200.0]");
  //setParamState("alphaCB", "[1.5,1.0,3.0]");
  setParamState("alphaCB", "[1.5,0.1,3.0]");
  //setParamState("nCB", "[9.0,0.01,30.0]");
  setParamState("nCB", "[10.0]");
  setParamState("muGANom", "[125,0.1,10000]");
  setParamState("sigmaGANom", "[10.0,0.01,80.0]");
  setParamState("fracCB", "[0.9,0.001,0.999]");
  
  // Double-sided Crystal Ball PDF:
  setParamState("alphaCBLo", "[1.5,0.1,2.5]");
  setParamState("nCBLo", "[17.0,0.01,30.0]");
  setParamState("alphaCBHi", "[2.2,0.1,3.0]");
  setParamState("nCBHi", "[5.2,0.01,10.0]");
  
  // Triple Gaussian PDF:
  setParamState("muGA1Nom", "[125,0.0,10000]");
  setParamState("muGA2Nom", "[125,0.0,10000]");
  setParamState("muGA3Nom", "[125,0.0,10000]");
  setParamState("sigmaGA1Nom", "[2.0,0.01,200.0]");
  setParamState("sigmaGA2Nom", "[4.0,0.01,200.0]");
  setParamState("sigmaGA3Nom", "[6.0,0.01,200.0]");
  setParamState("fracGA1", "[0.5,0.001,0.999]");
  setParamState("fracGA2", "[0.4,0.001,0.999]");
  
  // Bifurcated Gaussian PDF:
  setParamState("sigmaGALowNom", "[3.0,0.01,200.0]");
  setParamState("sigmaGAHiNom", "[1.0,0.01,200.0]");
  
  // Landau parameterization:
  setParamState("muLANom", "[125,0.0,10000]");
  setParamState("sigmaLANom", "[2.0,0.01,200.0]");
  
  // Voigt-specific parameters:
  setParamState("muVoigtNom", "[125,0.0,10000]");
  setParamState("widthVoigtNom", "[10.0,0.01,80.0]");
  setParamState("sigmaVoigtNom", "[10.0,0.01,80.0]");
  
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
    
    // Create mass observable for this category:
    if (!m_ws->var(Form("m_yy_%s",currKey.Data()))) {
      m_ws->factory(Form("m_yy_%s[10,10000]",currKey.Data()));
    }
    
    RooDataSet* newData 
      = new RooDataSet(Form("data_%s",currKey.Data()),
		       Form("data_%s",currKey.Data()),
		       RooArgSet(*m_ws->var(Form("m_yy_%s",currKey.Data())),
				 *m_ws->var("wt")),
		       RooFit::WeightVar(*m_ws->var("wt")));
    m_ws->import(*newData);
    
    // Each dataset also corresponds to a unique category in the fit:
    m_cat->defineType(currKey);
    
    // Keep track of all resonance mass - category pairs for fitting.
    std::pair<double,int> currPair;
    //currPair.first = resonanceMass
    currPair.first = massIntToDouble(massDoubleToInt(resonanceMass));
    currPair.second = cateIndex;
    m_massCatePairs.push_back(currPair);
    
  }
  
  // Set the observable and weight values and then fill dataset:
  m_ws->var(Form("m_yy_%s",currKey.Data()))->setVal(diphotonMass);
  m_ws->var("wt")->setVal(eventWeight);
  RooDataSet* currData =(RooDataSet*)m_ws->data(Form("data_%s",currKey.Data()));
  currData->add(RooArgSet(*m_ws->var(Form("m_yy_%s",currKey.Data())),
			  *m_ws->var("wt")),
		eventWeight);
  
  // Also update the number of categories:
  if (cateIndex >= m_nCategories) {
    m_nCategories = cateIndex + 1;
  }
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
   @return - True iff the PDF and yield parameter were imported.
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
   @return - True iff the PDF and yield parameter were imported.
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
   -----------------------------------------------------------------------------
   Adds a parameter to the workspace for the fit.
   @param varName - The name of the shape variable of interest.
   @param cateIndex - The event category index (starting at 0).
*/
void SigParam::addVariable(TString varName, int cateIndex) {
  if (m_verbose) {
    std::cout << "SigParam: addVariable(" << varName << ", " << cateIndex << ")"
	      << std::endl;
  }
  std::vector<TString> currParams = listParamsForVar(varName);
  // Loop over the parameterization variables corresponding to varName.
  for (int i_p = 0; i_p < (int)currParams.size(); i_p++) {
    TString valAndRange = getParamState(currParams[i_p]);
    TString fullName = Form("%s_%sc%d", (currParams[i_p]).Data(),
			    m_signalType.Data(), cateIndex);
    // Check to see whether the variable has already been added before creating:
    if (!m_ws->var(fullName)) {
      m_ws->factory(Form("%s%s", fullName.Data(), valAndRange.Data()));
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Creates a single or combined binned dataset. 
   @param unbinnedName - The name of the unbinned data in the workspace.
   @param resonanceMass - The mass value.
   @param cateIndex - The index of the category.
*/
void SigParam::binTheData(TString unbinnedName, double resonanceMass, 
			  int cateIndex) {
  // Create a single binned dataset:
  if (!unbinnedName.Contains(Form("data_c%d",cateIndex))) {
    TString currBinnedName = unbinnedName; 
    currBinnedName.ReplaceAll("data","dataBinned");
    binSingleDataSet(unbinnedName,currBinnedName, resonanceMass, cateIndex);
  }
  // Create the combined binned dataset:
  else {
    RooArgSet *args = new RooArgSet();
    
    // Create a dataset map to store the RooDataSet objects:
    std::map<std::string,RooDataSet*> currDataMap; currDataMap.clear();
    
    // Create a list of the corresponding unbinned datasets, and loop over it:
    std::vector<double> currMassPoints = massPointsForCategory(cateIndex); 
    for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
      TString currKey = getKey(currMassPoints[i_m], cateIndex);

      // Bin each one, add to map
      TString currUnbinnedName = Form("data_%s", currKey.Data());
      TString currBinnedName = Form("dataBinned_%s", currKey.Data());
      binSingleDataSet(currUnbinnedName, currBinnedName, currMassPoints[i_m],
		       cateIndex);
      currDataMap[((std::string)currKey)]
	= (RooDataSet*)m_ws->data(currBinnedName);
      args->add(*m_ws->var(Form("m_yy_%s",currKey.Data())));
    }
        
    args->add(*(m_ws->var("wt")));
    RooDataSet *dataBinned =new RooDataSet(Form("dataBinned_c%d", cateIndex),
					   Form("dataBinned_c%d", cateIndex),
					   *args, RooFit::Index(*m_cat),
					   RooFit::Import(currDataMap),
					   RooFit::WeightVar(*m_ws->var("wt")));
    m_ws->import(*dataBinned);
  }
}

/**
   -----------------------------------------------------------------------------
   Create a binned dataset from an unbinned dataset.
   @param unbinnedName - The name of the unbinned data in the workspace.
   @param binnedName - The name of the binned data to create.
   @param resonanceMass - The mass value.
   @param cateIndex - The index of the category.
*/
void SigParam::binSingleDataSet(TString unbinnedName, TString binnedName,
				double resonanceMass, int cateIndex) {
  if (m_verbose) {
    std::cout << "SigParam: Rebin sample " << unbinnedName << ", category "
	      << cateIndex << " for mass " << resonanceMass << std::endl;
  }
  
  RooDataSet *unbinnedData = (RooDataSet*)m_ws->data(unbinnedName);
  RooRealVar *myyCurr
    = m_ws->var(Form("m_yy_%s", (getKey(resonanceMass, cateIndex)).Data()));
  
  // Create a histogram to automatically bin the points:
  int nBins = (int)(m_nBinsPerGeV*(myyCurr->getMax() - myyCurr->getMin()));
  TH1 *hist = unbinnedData
    ->createHistogram(Form("hist_%s",unbinnedName.Data()), *myyCurr, 
		      RooFit::Binning(nBins, myyCurr->getMin(),
				      myyCurr->getMax()));
  
  // Create a dataset to fill with binned entries:
  RooDataSet* binnedData = new RooDataSet(binnedName, binnedName,
					  RooArgSet(*myyCurr, *m_ws->var("wt")),
					  RooFit::WeightVar(*m_ws->var("wt")));
  for (int i_b = 1; i_b <= hist->GetNbinsX(); i_b++) {
    myyCurr->setVal(hist->GetXaxis()->GetBinCenter(i_b));
    m_ws->var("wt")->setVal(hist->GetBinContent(i_b));
    binnedData->add(RooArgSet(*myyCurr, *m_ws->var("wt")), 
		    hist->GetBinContent(i_b));
  }
  
  // Import binned dataset to the workspace:
  m_ws->import(*binnedData);
}

/**
   -----------------------------------------------------------------------------
   Calculate the width of the shape containing 68.2% of the events.
   @param resonanceMass - The mass value.
   @param cateIndex - The index of the category.
   @param mean - The peak value of the PDF.
*/
double SigParam::calculateStdDev(double resonanceMass, int cateIndex) {
  // The precision should be 0.001%:
  double precision = 0.00001;

  TString currKey = getKey(resonanceMass, cateIndex);
  
  // Load the mass variable and range:
  RooRealVar *observable = m_ws->var(Form("m_yy_%s", currKey.Data()));
  double minOrigin = observable->getMin();
  double maxOrigin = observable->getMax();
  
  // Use the parameterized PDF if it exists, otherwise the individual PDF:
  RooAbsPdf *currPdf = NULL;
  if ((m_ws->pdf(Form("sigPdf_%sc%d",m_signalType.Data(),cateIndex)))) {
    (*m_ws->var("mResonance")).setVal(resonanceMass);    
    currPdf = (m_ws->pdf(Form("sigPdf_%sc%d",m_signalType.Data(),cateIndex)));
  }
  else if((m_ws->pdf(Form("sigPdf_%s%s",m_signalType.Data(),currKey.Data())))){
    currPdf
      = (m_ws->pdf(Form("sigPdf_%s%s",m_signalType.Data(),currKey.Data())));
  }
  
  // Get the maximum possible signal width, and also the mean:
  double mean = getMeanOrStdDev("Mean", resonanceMass, cateIndex);
  double stdDev = TMath::Min(fabs(mean-minOrigin), fabs(maxOrigin-mean));
  
  // Calculate the initial integral using the full mass range:
  if (!m_verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  }
  observable->setRange("rangeFull", minOrigin, maxOrigin);
  RooAbsReal* intTot = (RooAbsReal*)currPdf
    ->createIntegral(RooArgSet(*observable), RooFit::NormSet(*observable), 
		     RooFit::Range("rangeFull"));
  double integralValTotal = intTot->getVal();
  
  // Initial values for the integral and step size:
  double integralVal = 1.0;
  double stepSize = stdDev;
  int nIterations = 0;
  double target = 0.682;
  int maxIterations = 30;
  // Iteratively find the STDDEV:
  while ((fabs(integralVal - target)/target) > precision && 
	 nIterations <= maxIterations) {
    
    // Decrease the step size by half each time:
    nIterations++;
    stepSize = 0.5 * stepSize;
    
    // STDDEV too large: decrease the width:
    if (integralVal > target) stdDev = stdDev - stepSize;
    
    // STDDEV too small: increase the width:
    else stdDev = stdDev + stepSize;
    
    // Then calculate the fractional integral value:
    if (!m_verbose) {
      RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    }
    observable->setRange(Form("range%f",stdDev), mean-stdDev, mean+stdDev);
    RooAbsReal* intCurr = (RooAbsReal*)currPdf
      ->createIntegral(RooArgSet(*observable), 
		       RooFit::NormSet(*observable), 
		       RooFit::Range(Form("range%f",stdDev)));
    integralVal = intCurr->getVal() / integralValTotal;
  }
  
  if (m_verbose) {
    std::cout << "SigParam: The STDDEV of the resonance is " << stdDev 
	      << " after " << nIterations << " iterations." << std::endl;
  }
  
  observable->setMin(minOrigin);
  observable->setMax(maxOrigin);
  return stdDev;
}

/**
   -----------------------------------------------------------------------------
   Get a list of categories corresponding to a single mass point.
   @param resonanceMass - The mass value.
   @return - A vector of category indices.
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
   Create a TF1 function object based on the results of a parameterized fit.
   @param varName - The name of the fit variable that has been parameterized.
   @param cateIndex - The index of the category.
   @param xMin - The minimum value of the function.
   @param xMax - The maximum value of the function.
   @return - A TF1 function with the fitted parameter values.
*/
TF1 *SigParam::createTF1FromParameterization(TString varName, int cateIndex,
					     double xMin, double xMax) {
  // Modify the formatting of the function as defined for RooFit:
  TString currFuncFormat = getVarParameterization(varName);
  currFuncFormat = currFuncFormat.ReplaceAll("obs", "((x-100.0)/100.0)");
  currFuncFormat = currFuncFormat.ReplaceAll("mRes", "x");
  for (int i_p = 0; i_p < getNParamsForVar(varName); i_p++) {
    currFuncFormat
      = currFuncFormat.ReplaceAll(Form("@%d",i_p), Form("[%d]",i_p));
  }
  // For variables that are not parameterized:
  if (currFuncFormat.EqualTo("")) currFuncFormat = "[0]";
  // Instantiate a new TF1 based on the function format above.
  TF1 *currFunc
    = new TF1(Form("fit_%s",varName.Data()), currFuncFormat, xMin, xMax);
  
  // Then set the parameter values:
  std::vector<TString> paramList = listParamsForVar(varName);
  for (int i_p = 0; i_p < (int)paramList.size(); i_p++) {
    currFunc->SetParameter(i_p, getParameterValue(paramList[i_p], cateIndex));
  }
  // And return the TF1:
  return currFunc;
}

/**
   -----------------------------------------------------------------------------
   Check if the dataset being requested has been instantiated (exists in map).
   @param massIndex - The index of the signal mass.
   @param cateIndex - The index of the category.
   @return - True iff the dataset has been defined.
*/
bool SigParam::dataExists(double resonanceMass, int cateIndex) {
  if (m_ws->data(Form("data_%s",(getKey(resonanceMass,cateIndex)).Data()))) {
    return true;
  }
  else return false;
}

/**
   -----------------------------------------------------------------------------
   This method generates and fits toy MC generated from a dataset, and returns
   the profiled Higgs mass and normalization as well as the truth values.
   @param resonanceMass - The floating value of the signal mass.
   @param cateIndex - The index of the category.
   @param dataType - "asimov", "mctoy", "pdftoy".
   @param seed - The random seed for pseudo-data generation (not for Asimov).
   @return - A vector: [0] = the truth resonance mass, [1] = the profiled 
   resonance mass, [2] = truth data normalization, [3] = profiled normalization.
*/
std::vector<double> SigParam::doBiasTest(double resonanceMass, int cateIndex,
					 TString dataType, int seed) {
  
  // First generate the desired data:
  generateData(resonanceMass, cateIndex, dataType, seed);
  
  // Free the parameters before fitting:
  TString currKey = getKey(resonanceMass, cateIndex);
  TString pdfName = Form("sigPdf_%sc%d", m_signalType.Data(), cateIndex);
  if (m_ws->pdf(pdfName)) {
    setParamsConstant(m_ws->pdf(pdfName), true);
    setResMassConstant(false, resonanceMass);
  }
  else {
    //pdfName = Form("sigPdf_%s%s", m_signalType.Data(), currKey.Data());
    std::cout << "SigParam: ERROR! doBiasTest() method currently only configured for parameterized shape, where mResonance can be free for the fit." 
	      << std::endl;
    exit(0);
  }
    
  // Then fit that dataset:
  RooFitResult* currFitResult
    = fitResult(resonanceMass, cateIndex, dataType, "ExtendedParameterized");
  
  // Return the fit bias information:
  std::vector<double> biasResult; biasResult.clear();
  
  // Bad fits have null pointer to status, nonzero value, or infinite/nan NLL
  if (currFitResult && currFitResult->status() == 0 &&
      std::isfinite(m_currNLL)) {
    
    biasResult.push_back(resonanceMass);
    biasResult.push_back(m_ws->var("mResonance")->getVal());
    biasResult.push_back(generatedDataNorm());
    biasResult.push_back(extendedTerm());
  }
  
  return biasResult;
}

/**
   -----------------------------------------------------------------------------
   Option to do a binned fit. Sets the private variable m_binned.
   @param binned - True iff the fit should be binned.
   @param nBinsPerGeV - The number of bins per GeV for the binned data.
*/
void SigParam::doBinnedFit(bool binned, double nBinsPerGeV) {
  m_binned = binned;
  m_nBinsPerGeV = nBinsPerGeV;
  if (m_verbose) {
    std::cout << "SigParam: Binned bool = " << m_binned << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Check if two doubles are equal.
   @param massValue1 - The first mass value to compare.
   @param massValue2 - The second mass value to compare.
   @return - True iff the masses are equal within 0.001 GeV.
*/
bool SigParam::equalMasses(double massValue1, double massValue2) {
  return (fabs(massValue1 - massValue2) <= 0.001);// mass precision in GeV
}

/**
   -----------------------------------------------------------------------------
   Retrieve the extended term value (normalization) from the most recent fit.
   @return - The value of the extended term from the last fit.
*/
double SigParam::extendedTerm() {
  return m_currExtendVal;
}

/**
   -----------------------------------------------------------------------------
   Perform a single or simultaneous fit.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category.
   @param dataType - The type of data being fitted. "" for nominal signal MC,
   "asimov", "mctoy", or "pdftoy" for generated data types.
   @param option - Fit options (e.g. "Extended", "Parameterized").
   @return - The RooFitResult, which gives fit status.
*/
RooFitResult* SigParam::fitResult(double resonanceMass, int cateIndex,
				  TString dataType, TString option) {
  if (m_verbose) {
    std::cout << "SigParam: Preparing to fit resonance" << std::endl;
  }
  
  // Clock the fit:
  clock_t time;
  time = clock();
  
  // Define the PDF and dataset names:
  TString currKey = getKey(resonanceMass,cateIndex);
  TString sigName = "";
  if (resonanceMass < 0.0) {
    sigName = Form("sigPdfTmp_%sc%d", m_signalType.Data(), cateIndex);
  }
  else if (option.Contains("Parameterized")) {
    sigName = Form("sigPdf_%sc%d", m_signalType.Data(), cateIndex);
  }
  else {
    sigName = Form("sigPdf_%s%s", m_signalType.Data(), currKey.Data());
  }
  
  // The default dataset (for "" dataType):
  TString dataName = (resonanceMass < 0.0) ? Form("data_c%d",cateIndex) :
    Form("data_%s",currKey.Data());  
  
  // Change the dataset name if the fit is binned:
  if (m_binned) {
    TString binnedDataName = dataName;
    binnedDataName.ReplaceAll("data", "dataBinned");
    
    // Also check that binned data have been created, otherwise create:
    if (m_ws->data(dataName) && !m_ws->data(binnedDataName)) {
      binTheData(dataName, resonanceMass, cateIndex);
    }
    dataName = binnedDataName;
  }
  
  // The section below modifies the PDF and dataset names for fitting generated
  // data such as asimov data or toy data:
  if (dataType.Contains("asimov") || dataType.Contains("mctoy") ||
      dataType.Contains("pdftoy")) {
    // Must specify a mass value for fitting generated datasets:
    if (resonanceMass < 0.0) {
      std::cout << "SigParam: ERROR! Must specify a mass to fit for generated data. Use the same mass that was used to generate Asimov or toy MC data." 
		<< std::endl;
      exit(0);
    }
    // Try to choose a parameterized PDF, then use individual if that fails:
    sigName = Form("sigPdf_%sc%d", m_signalType.Data(), cateIndex);
    if (!m_ws->pdf(sigName)) {
      sigName = Form("sigPdf_%s%s", m_signalType.Data(), currKey.Data());
    }
    dataName = Form("data_%s_%s", dataType.Data(), currKey.Data());
  }
  
  // Now create pointers to chosen signal and data ahead of fitting!
  RooDataSet *currData = (RooDataSet*)m_ws->data(dataName);
  RooAbsPdf *currSignal = m_ws->pdf(sigName);
  
  // In case an extended fit is requested:
  RooExtendPdf *currExtend = NULL; 
  RooRealVar *currNorm = NULL;
  if (option.Contains("Extended")) {
    currNorm = new RooRealVar("currNorm", "currNorm", 100.0, 0.0, 1000000.0);
    currExtend = new RooExtendPdf("extendedPdf", "extendedPdf", 
				  *currSignal, *currNorm);
  }
  
  // Make sure inputs are defined, else exit:
  if (!currSignal) {
    std::cout << "SigParam: ERROR! Signal for fit not defined: "
	      << sigName << std::endl;
    exit(0);
  }
  if (!currData) {
    std::cout << "SigParam: ERROR! Data for fit not defined: " 
	      << dataName << std::endl;
    exit(0);
  }

  int fitPrintLevel = m_verbose ? 0 : -1;
  RooFitResult *result = NULL;
  // Individual fit: apply the mass range requirement.
  if (resonanceMass > 0 && !option.Contains("Parameterized")) {
    double fitMin = m_ws->var(Form("m_yy_%s",currKey.Data()))->getMin();
    double fitMax = m_ws->var(Form("m_yy_%s",currKey.Data()))->getMax();
    if (option.Contains("Extended")) {
      result = currExtend->fitTo(*currData, RooFit::PrintLevel(fitPrintLevel),
				 RooFit::SumW2Error(kTRUE), RooFit::Save(true),
				 RooFit::Range(fitMin,fitMax));
      m_currExtendVal = currNorm->getVal();
    }
    else {
      
      result = currSignal->fitTo(*currData, RooFit::PrintLevel(fitPrintLevel),
				 RooFit::SumW2Error(kTRUE), RooFit::Save(true),
				 RooFit::Range(fitMin,fitMax));
    }
  }
  // Parameterized fit: no explicit mass range requirement.
  else {
    if (option.Contains("Extended")) {
      result = currExtend->fitTo(*currData, RooFit::PrintLevel(fitPrintLevel),
				 RooFit::SumW2Error(kTRUE), RooFit::Save(true));
      m_currExtendVal = currNorm->getVal();
    }
    else {
      result = currSignal->fitTo(*currData, RooFit::PrintLevel(fitPrintLevel),
				 RooFit::SumW2Error(kTRUE), RooFit::Save(true));
    }
  }
  
  m_currNLL = result->minNll();
  return result;
  
  // Fix the fit parameters:
  SigParam::setParamsConstant(currSignal, true);
  
  // Clock the fit:
  time = clock() - time;
  if (m_verbose) {
    std::cout << "SigParam: Fit procedure concluded." << std::endl;
    printf("\tFit required %d clock cycles (%f seconds).\n\n",
	   (int)time, ((float)time/CLOCKS_PER_SEC));
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Perform a simultaneous fit across multiple masses.
   @param cateIndex - The index of the category.
   @param dataType - The type of data being fitted. "" for nominal signal MC,
   "asimov", "mctoy", or "pdftoy" for generated data types.
   @param option - Fit options (e.g. "Extended").
   @return - The RooFitResult, which gives fit status.
*/
RooFitResult* SigParam::fitResult(int cateIndex, TString dataType,
				  TString option) {
  return SigParam::fitResult(-999.9, cateIndex, dataType, option);
}

/**
   -----------------------------------------------------------------------------
   Check whether the given function has been implemented in this class.
   @param function - The name of the function to test for definition.
   @return - True iff. the function has been defined.
*/
bool SigParam::functionIsDefined(TString function) {
  if (function.Contains("DoubleCB") || function.Contains("CBGA") ||
      function.Contains("GAx3") || function.Contains("BifurGA") ||
      function.Contains("Landau") || function.Contains("Voigt")) {
    return true;
  }
  else {
    std::cout << "SigParam: ERROR! " << function << " is undefined. Please use one of the following: DoubleCB, CBGA, GAx3, BifurGA, Landau, Voigt, CBPlusVoigt..." 
	      << std::endl;
    return false;
  }
}

/**
   -----------------------------------------------------------------------------
   This method generates and fits either Asimov data, toy data from the input 
   MC, or toy data from the fitted PDF. 
   @param resonanceMass - The floating value of the signal mass.
   @param cateIndex - The index of the category.
   @param type - "asimov", "mctoy", "pdftoy".
   @param seed - The random seed for pseudo-data generation (not for Asimov).
   @return - True iff. generation works and fit converges.
*/
bool SigParam::generateAndFitData(double resonanceMass, int cateIndex,
				  TString dataType, int seed) {
  
  // First generate the desired data:
  generateData(resonanceMass, cateIndex, dataType, seed);
  
  // Free the parameters before fitting:
  TString currKey = getKey(resonanceMass, cateIndex);
  TString pdfName = Form("sigPdf_%sc%d", m_signalType.Data(), cateIndex);
  if (m_ws->pdf(pdfName)) setResMassConstant(false, resonanceMass);
  else pdfName = Form("sigPdf_%s%s", m_signalType.Data(), currKey.Data());
  setParamsConstant(m_ws->pdf(pdfName), false);
    
  // Then fit that dataset:
  RooFitResult* result
    = fitResult(resonanceMass, cateIndex, dataType, "ExtendedParameterized");
  
  // Return the fit status:
  // Bad fits have null pointer to status, nonzero value, or infinite/nan NLL
  return (result && result->status() == 0 && std::isfinite(m_currNLL));
}

/**
   -----------------------------------------------------------------------------
   This method generates either Asimov data, toy data from the input MC, or toy
   data from the fitted PDF. 
   @param resonanceMass - The floating value of the signal mass.
   @param cateIndex - The index of the category.
   @param type - "asimov", "mctoy", "pdftoy".
   @param seed - The random seed for pseudo-data generation (not for Asimov).
   @return - A pointer to the new RooAbsData set. 
*/
RooDataSet* SigParam::generateData(double resonanceMass, int cateIndex, 
				   TString dataType, int seed) {
  if (m_verbose) {
    std::cout << "SigParam: Creating Asimov data mRes=" << resonanceMass
	      << ", cate=" << cateIndex << ", type=" << dataType << std::endl;
  }
  
  // The dataset to be generated and returned:
  RooDataSet *generatedData = NULL;
  
  // Get the name of the dataset:
  TString currKey = getKey(resonanceMass, cateIndex);
  TString dataName = Form("data_%s", currKey.Data());
    
  // Use simultaneous model if it exists, otherwise use individual model:
  TString obsName;
  TString pdfName = Form("sigPdf_%sc%d", m_signalType.Data(), cateIndex);
  if (m_ws->pdf(pdfName)) {
    // Need to set the resonance mass to the correct value:
    setResMassConstant(true, resonanceMass);
    // Also use the generic m_yy observable:
    obsName = "m_yy";
  }
  // Individual model must be used:
  else {
    pdfName = Form("sigPdf_%s%s", m_signalType.Data(), currKey.Data());
    obsName = Form("m_yy_%s", currKey.Data());
  }
  
  // Get the yield at this mass and category for generation:
  double nEventsToGenerate = getYieldInCategory(resonanceMass, cateIndex);
  
  //--------------------//
  // If Asimov data are requested:
  if (dataType.EqualTo("asimov")) {
    // Generate the Asimov data using the RooStats method:
    m_ws->factory(Form("SUM::sigPdfAsimov(nSigAsimov[%f]*%s)",
		       nEventsToGenerate, pdfName.Data()));
    generatedData = (RooDataSet*)RooStats::AsymptoticCalculator::
      GenerateAsimovData(*m_ws->pdf("sigPdfAsimov"),
			 RooArgSet(*m_ws->var(obsName)));
    //GenerateAsimovData(*m_ws->pdf(pdfName), RooArgSet(*m_ws->var(obsName)));
    
    // Then add ghost events:
    double min = m_ws->var(obsName)->getMin();
    double max = m_ws->var(obsName)->getMax();
    double ghostWt = 0.0000001;
    for (double massVal = min; massVal <=max; massVal += fabs((max-min)/10.0)) {
      m_ws->var(obsName)->setVal(massVal);
      generatedData->add(RooArgSet(*m_ws->var(obsName)), ghostWt);
    }
  }
  
  //--------------------//
  // If toy data from the fitted PDF are requested:
  else if(dataType.EqualTo("pdftoy")) {
    RooRandom::randomGenerator()->SetSeed(seed);
    generatedData = (RooDataSet*)m_ws->pdf(pdfName)
      ->generate(*m_ws->var(obsName),nEventsToGenerate,RooFit::Extended(true));
  }

  //--------------------//
  // If toy data from the input dataset are requested:
  else if (dataType.EqualTo("mctoy")) {
    // Check that the underlying dataset exists!
    if (!m_ws->data(dataName)) {
      std::cout << "SigParam: ERROR! dataset required but missing: " << dataName
		<< std::endl;
      exit(0);
    }
    
    // Instantiate the dataset to be returned:
    generatedData
      = new RooDataSet(Form("data_mctoy_%s",currKey.Data()),
		       Form("data_mctoy_%s",currKey.Data()),
		       RooArgSet(*m_ws->var(obsName),*m_ws->var("wt")),
		       RooFit::WeightVar(*m_ws->var("wt")));
    
    // Turn original dataset into histogram:
    double min = m_ws->var(obsName)->getMin();
    double max = m_ws->var(obsName)->getMax();
    int bins = (int)(m_nBinsPerGeV * (max - min));
    TString originObsName = Form("m_yy_%s", currKey.Data());
    TH1F *dataHist = (TH1F*)m_ws->data(dataName)
      ->createHistogram("histTmp", *m_ws->var(originObsName),
			RooFit::Binning(bins,min,max));
    
    // Randomly generate poisson normalization:
    RooRandom::randomGenerator()->SetSeed(seed);
    TRandom *random = new TRandom(); random->SetSeed(seed);
    int integerEventsToGenerate = random->Poisson(nEventsToGenerate);
    for (int i_r = 0; i_r < integerEventsToGenerate; i_r++) {
      // Generate random events from input data histogram:
      double currVal = dataHist->GetRandom();
      double currWt = 1.0;
      // Then add to new RooDataSet:
      m_ws->var(obsName)->setVal(currVal);
      m_ws->var("wt")->setVal(currWt);
      generatedData
	->add(RooArgSet(*m_ws->var(obsName), *m_ws->var("wt")), currWt);
    }
    
    // Then add ghost events:
    double ghostWt = 0.0000001;
    for (double massVal = min; massVal <=max; massVal += fabs((max-min)/10.0)) {
      m_ws->var(obsName)->setVal(massVal);
      m_ws->var("wt")->setVal(ghostWt);
      generatedData
	->add(RooArgSet(*m_ws->var(obsName), *m_ws->var("wt")), ghostWt);
    }
    delete dataHist;
  }
  
  //--------------------//
  // Option not recognized:
  else {
    std::cout << "SigParam: ERROR! GenerateData doesn't recognize type " 
	      << dataType << std::endl;
    exit(0);
  }
  
  if (m_verbose) {
    std::cout << "SigParam: generated data has " << generatedData->sumEntries() 
	      << " entries (weighted)" << std::endl;
  }
  
  // Rename the dataset:
  generatedData->SetNameTitle(Form("data_%s_%s", dataType.Data(), 
				   currKey.Data()),
			      Form("data_%s_%s", dataType.Data(),
				   currKey.Data()));
  // Save norm. so that it is retrievable via generatedDataNormalization():
  m_generatedDataNorm = generatedData->sumEntries();
  
  // Import and return new dataset:
  m_ws->import(*generatedData);
  return generatedData;
}

/**
   -----------------------------------------------------------------------------
   Get the weighted normalization of the dataset that was just generated.
   @return - The weighted sum of entries or the number of entries.
*/
double SigParam::generatedDataNorm() {
  return m_generatedDataNorm;
}

/**
   -----------------------------------------------------------------------------
   Retrieve a key string for the given mass and category index.
   @param resonanceMass - The floating value of the signal mass.
   @param cateIndex - The index of the category.
   @return - A key string specific to the mass and category.
*/
TString SigParam::getKey(double resonanceMass, int cateIndex) {
  TString key = Form("m%d_c%d", massDoubleToInt(resonanceMass), cateIndex);
  return key;
}

/**
   -----------------------------------------------------------------------------
   Get the mean value of the fit shape. If multiple functions are used, the 
   weighted mean is returned, where the weight is determined by the relative
   normalization of each component. Note: the StdDev is calculated via
   integration, and is not simply the resolution parameter for a given PDF.
   @param value - "Mean" or "StdDev".
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category for the PDF.
   @return - The weighted mean (mu) of the signal shape.
*/
double SigParam::getMeanOrStdDev(TString value, double resonanceMass,
				int cateIndex) {
  if (value.EqualTo("StdDev")) {
    return calculateStdDev(resonanceMass, cateIndex);
  }
  else if (value.EqualTo("Mean")) {
    std::vector<double> valList; valList.clear();
    std::vector<double> fracList; fracList.clear();
    std::vector<TString> currVars = variablesForFunction(m_currFunction);
    for (int i_v = 0; i_v < (int)currVars.size(); i_v++) {
      if (currVars[i_v].Contains("mu") && value.EqualTo("Mean")) {
	valList.push_back(getParameterValue(currVars[i_v], resonanceMass, 
					    cateIndex));
      }
      else if (currVars[i_v].Contains("frac")) {
	fracList.push_back(getParameterValue(currVars[i_v], resonanceMass,
					     cateIndex));
      }
    }
    
    // Then loop over the muList, and create a weighted sum of the mean values,
    // with the weight given by the "frac" variables:
    double result = 0.0;
    double priorFracs = 0.0;
    if (m_currFunction.Contains("BifurGA")) fracList.push_back(0.5);
    for (int i_m = 0; i_m < (int)valList.size(); i_m++) {
      if (i_m < (int)valList.size()-1) { 
	result += (valList[i_m] * fracList[i_m]);
	priorFracs += fracList[i_m];
      }
      else {
	result += (valList[i_m] * (1.0-priorFracs));
      }
    }
    return result;
  }
  else {
    std::cout << "SigParam: ERROR! value " << value << " undefined." 
	      << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get the mean value of the data used to fit the PDFs. 
   @param value - "Mean" or "StdDev".
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category for the PDF.
   @return - The weighted mean (mu) of the input data.
*/
double SigParam::getMeanOrStdDevInData(TString value, double resonanceMass,
				       int cateIndex) {
  TString currKey = getKey(resonanceMass, cateIndex);
  if (value.Contains("Mean")) {
    return m_ws->data(Form("data_%s",currKey.Data()))
      ->mean(*m_ws->var(Form("m_yy_%s",currKey.Data())));
  }
  else if (value.Contains("StdDev")) {
    return m_ws->data(Form("data_%s",currKey.Data()))
      ->sigma(*m_ws->var(Form("m_yy_%s",currKey.Data())));
  }
  else {
    std::cout << "SigParam: ERROR! value " << value << " undefined." 
	      << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get the number of categories contained in the datasets for fitting. Note:
   it is possible that there are different numbers of categories defined for 
   different mass points. This is up to the user to sort out. 
   @return - The total number of categories for the parameterization.
*/
int SigParam::getNCategories() {
  return m_nCategories;
}

/**
   -----------------------------------------------------------------------------
   Get the number of parameters used to parameterize the given variable.
   @param varName - The name of the fit variable that has been parameterized.
   @return - The number of parameters.
*/
int SigParam::getNParamsForVar(TString varName) {
  TString currFunction = getVarParameterization(varName);
  if (currFunction.EqualTo("")) return 0;
  // TString CountChar method requires a bug report for ROOT.
  //else return currFunction.CountChar("@");
  else {
    int counter = 0;
    for (Ssiz_t i_c = 0; i_c < currFunction.Length(); i_c++) {
      if (TString(currFunction[i_c]).EqualTo("@")) counter++;
    }
    return counter;
  }
}

/**
   -----------------------------------------------------------------------------
   Get the value of the fit error for a particular parameter of the signal PDF. 
   @param paramName - The name of the shape parameter of interest.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category.
   @return - The value of the specified signal parameter. 
*/
double SigParam::getParameterError(TString paramName, double resonanceMass,
				   int cateIndex) {
  if (paramName.Contains("nCB") || paramName.Contains("frac")){
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
   @return - The value of the specified signal parameter. 
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
   @return - The value of the specified signal parameter. 
*/
double SigParam::getParameterValue(TString paramName, double resonanceMass, 
				   int cateIndex) {
  if (paramName.Contains("nCB") || paramName.Contains("frac")) {
    return SigParam::getParameterValue(paramName, cateIndex);
  }
  else {
    RooRealVar *var = m_ws->var(Form("%s_%s%s", paramName.Data(), 
				     m_signalType.Data(),
				     (getKey(resonanceMass,cateIndex)).Data()));
    if (!var) {
      std::cout << "SigParam: requested parameter not found: param = "
		<< paramName << ", mass = " << resonanceMass << ", cate = "
		<< cateIndex << std::endl;
      exit(0);//return 0.0;
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
   @return - The value of the specified signal parameter. 
*/
double SigParam::getParameterValue(TString paramName, int cateIndex) {
  RooRealVar *var = m_ws->var(Form("%s_%sc%d", paramName.Data(),
				   m_signalType.Data(), cateIndex));
  if (!var) {
    std::cout << "SigParam: requested parameter not found: param = "
	      << paramName << ", cate = " << cateIndex << std::endl;
    exit(0);//return 0.0;
  }
  else {
    return var->getVal();
  }
}

/**
   -----------------------------------------------------------------------------
   Get the parameterized value of a particular parameter of the signal PDF. This
   method only works for parameterized fits.
   @param paramName - The name of the shape parameter of interest.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category.
   @return - The value of the specified signal parameter. 
*/

double SigParam::getParameterizedValue(TString paramName, double resonanceMass, 
				       int cateIndex) {
  m_ws->var("mResonance")->setVal(resonanceMass);
  TString funcName = Form("%s_%sc%d", paramName.Data(), m_signalType.Data(),
			  cateIndex);
  funcName.ReplaceAll("Nom","");
  if (m_ws->function(funcName)) return m_ws->function(funcName)->getVal();
  else {
    std::cout << "SigParam: ERROR! Unrecognized parameterization variable "
	      << funcName << ". Returning nominal variable value instead."
	      << std::endl;
    return getParameterValue(paramName, resonanceMass, cateIndex);
  }
}

/**
   -----------------------------------------------------------------------------
   Get the initial value [a] and range [b,c] of a parameter: "[a,b,c]"
   @param paramName - The name of the parameter of interest.
   @return - The initial value and range of the parameter. 
*/
TString SigParam::getParamState(TString paramName) {
  return m_paramState[paramName];
}

/**
   -----------------------------------------------------------------------------
   Retrieves the resonance parameterized as a function of mResonance.
   @param cateIndex - The index of the category for the PDF.
   @return - A pointer to the signal PDF.
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
   @return - A pointer to the signal PDF.
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
   Get the value of a test statistic from the most recent fit in the specified
   category. Chi2 is only available for binned fits, whereas NLL is available 
   for any fit. The fits are all done by minimizing -log(L).
   @param statistic - The name of the test statistic "Chi2" or "NLL". 
   @param cateIndex - The index of the category that was fit.
   @return - The value of the test statistic.
*/
double SigParam::getTestStat(TString statistic, int cateIndex) {
  if (statistic.EqualTo("Chi2") && !m_binned) {
    std::cout << "SigParam: Error! Chi^2 can only be calculated for binned fit."
	      << std::endl;
    exit(0);
  }
  return m_testStats[Form("%s_c%d", statistic.Data(), cateIndex)];
}

/**
   -----------------------------------------------------------------------------
   Get the value of a test statistic from the most recent fit in the specified
   category and resonance mass. Chi2 is only available for binned fits, whereas 
   NLL is available for any fit. The fits are all done by minimizing -log(L).
   @param statistic - The name of the test statistic "Chi2" or "NLL". 
   @param resonanceMass - The truth mass of the resonance that was fit.
   @param cateIndex - The index of the category that was fit.
   @return - The value of the test statistic.
*/
double SigParam::getTestStat(TString statistic, double resonanceMass,
			     int cateIndex) {
  if (statistic.EqualTo("Chi2") && !m_binned) {
    std::cout << "SigParam: Error! Chi^2 can only be calculated for binned fit."
	      << std::endl;
    exit(0);
  }
  TString currKey = getKey(resonanceMass,cateIndex);
  return m_testStats[Form("%s_%s", statistic.Data(), currKey.Data())];
}

/**
   -----------------------------------------------------------------------------
   Get the value of a test statistic from the most recent fit. Chi2 is only 
   available for binned fits, whereas NLL is available for any fit. The fits are
   all done by minimizing -log(L).
   @param statistic - The name of the test statistic "Chi2" or "NLL". 
   @return - The value of the test statistic from the single most recent fit.
*/
double SigParam::getTestStatLatest(TString statistic) {
  if (statistic.EqualTo("Chi2") && !m_binned) {
    std::cout << "SigParam: Error! Chi^2 can only be calculated for binned fit."
	      << std::endl;
    exit(0);
  }
  else if (statistic.EqualTo("Chi2")) {
    if (m_verbose) {
      std::cout << "SigParam: Most recent fit gives " << statistic << " = " 
		<< m_currChi2 << std::endl;
    }
    return m_currChi2;
  }
  else {
    if (m_verbose) {
      std::cout << "SigParam: Most recent fit gives " << statistic << " = " 
		<< m_currNLL << std::endl;
    }
    return m_currNLL;
  }
}

/**
   -----------------------------------------------------------------------------
   Get a list of parameters associated with the PDF in the given category.
   @param resonanceMass - The truth mass of the resonance
   @param cateIndex - The index of the category for which we want the PDF.
   @return - A vector of parameter names.
*/
std::vector<TString> SigParam::getVariableNames(double resonanceMass, 
						int cateIndex) {
  std::cout << "SigParam: Get PDF variables in category = " << cateIndex
	    << " and mass = " << resonanceMass << std::endl;
  std::vector<TString> result; result.clear();
  TString currKey = getKey(resonanceMass, cateIndex);
  TString pdfName = Form("sigPdf_%s%s", m_signalType.Data(), currKey.Data());
  RooArgSet* currSet = (*m_ws->pdf(pdfName)).getVariables();
  TIterator *iterArg = currSet->createIterator();
  RooRealVar *currArg = NULL;
  while ((currArg = (RooRealVar*)iterArg->Next())) {
    TString currArgName = currArg->GetName();
    currArgName = currArgName.ReplaceAll(currKey, "");
    currArgName = currArgName.ReplaceAll(Form("_%s",m_signalType.Data()), "");
    currArgName = currArgName.ReplaceAll(Form("c%d",cateIndex), "");
    if (!currArgName.Contains("m_yy")) result.push_back(currArgName);
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Retrieve the functional form of the parameterization in the resonance mass
   of the given fit variable.
   @param varName - The name of the fit variable that has been parameterized.
   @return - A string with the functional form of the parameterization.
*/
TString SigParam::getVarParameterization(TString varName) {
  std::map<TString,TString>::iterator varIter;
  for (varIter = m_varParameterization.begin(); 
       varIter != m_varParameterization.end(); varIter++) {
    if (varName.EqualTo(varIter->first)) {
      return varIter->second;
    }
  }
  if (m_verbose) {
    std::cout << "SigParam: WARNING No defined parameterization for "
	      << varName << std::endl;
  }
  return "";
}

/**
   -----------------------------------------------------------------------------
   Get the workspace that stores the SigParam class data.
   @return - A pointer to the workspace storing all datasets and PDFs.
*/
RooWorkspace* SigParam::getWorkspace() {
  return m_ws;
}

/**
   -----------------------------------------------------------------------------
   Get the signal yield for a particular mass in a particular category.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category for which we want the PDF.
   @return - The signal yield for the specified mass in the given category.
*/
double SigParam::getYieldInCategory(double resonanceMass, int cateIndex) {
  if (m_verbose) {
    std::cout << "SigParam: Get yield in category = " 
	      << cateIndex << " at mass = " << resonanceMass << std::endl;
  }
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
   Get the signal yield for a signal at a specified mass in a specified category
   inside a specified range of the observable.
   @param resonanceMass - The truth mass of the resonance.
   @param cateIndex - The index of the category for which we want the PDF.
   @param obsMin - The minimum of the observable range.
   @param obsMax - The maximum of the observable range. 
   @return - The signal yield in the window [obsMin,obsMax].
*/
double SigParam::getYieldInWindow(double resonanceMass, int cateIndex,
				  double obsMin, double obsMax) {
  TString pdfName = "";
  TString obsName = "";
  
  // Look for parameterized model:
  pdfName = Form("sigPdf_%sc%d", m_signalType.Data(), cateIndex);
  if (m_ws->pdf(pdfName)) {
    setResMassConstant(true, resonanceMass);
    obsName = Form("m_yy");
  }
  // If no parameterized model exists, get single point model:
  else {
    TString currKey = getKey(resonanceMass, cateIndex);
    pdfName = Form("sigPdf_%s%s", m_signalType.Data(), currKey.Data());
    obsName = Form("m_yy_%s", currKey.Data());
  }
  
  // Then check that both data and PDF exist:
  if (!m_ws->pdf(pdfName)) {
    std::cout << "SigParam: PDF for yield undefined." << std::endl;
    exit(0);
  }
  
  // Get the total integral of the PDF:
  double obsMinOrigin = m_ws->var(obsName)->getMin();
  double obsMaxOrigin = m_ws->var(obsName)->getMax();
  m_ws->var(obsName)->setRange("rangeTotal", obsMinOrigin, obsMaxOrigin);
  RooAbsReal* integralTotal = (RooAbsReal*)m_ws->pdf(pdfName)
    ->createIntegral(RooArgSet(*m_ws->var(obsName)), 
		     RooFit::NormSet(*m_ws->var(obsName)), 
		     RooFit::Range("rangeTotal"));
  double integralTotalVal = integralTotal->getVal();
  
  // Get the integral of the PDF in the window:
  m_ws->var(obsName)->setRange("rangeWindow", obsMin, obsMax);
  RooAbsReal* integralWindow = (RooAbsReal*)m_ws->pdf(pdfName)
    ->createIntegral(RooArgSet(*m_ws->var(obsName)), 
		     RooFit::NormSet(*m_ws->var(obsName)), 
		     RooFit::Range("rangeWindow"));
  double integralWindowVal = integralWindow->getVal();
  
  // Then reset observable range:
  m_ws->var(obsName)->setRange(obsMinOrigin, obsMaxOrigin);
  
  // Caculate and return yield:
  double yield = getYieldInCategory(resonanceMass, cateIndex);
  double fractionInWindow = integralWindowVal / integralTotalVal;
  double yieldInWindow = fractionInWindow * yield;
  return yieldInWindow;
}

/**
   -----------------------------------------------------------------------------
   Get the total inclusive signal yield for a signal at a specified mass in the
   specified observable window.
   @param resonanceMass - The truth mass of the resonance.
   @param obsMin - The minimum of the observable range.
   @param obsMax - The maximum of the observable range. 
*/
double SigParam::getYieldInWindow(double resonanceMass, double obsMin,	
				  double obsMax) {
  // Create a list of categories for this mass point:
  std::vector<int> currCategories = categoriesForMass(resonanceMass);
  // Loop through names of datasets, add components:
  double sum = 0.0;
  for (int i_c = 0; i_c < (int)currCategories.size(); i_c++) {
    sum += getYieldInWindow(resonanceMass, currCategories[i_c], obsMin, obsMax);
  }
  return sum;
}

/**
   -----------------------------------------------------------------------------
   Get the signal yield for a particular resonance mass in all categories.
   @param resonanceMass - The truth mass of the resonance.
   @return - The signal yield in all categories for the specified mass.
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
   Get a list of parameters required to parameterize the given variable as a 
   function of the resonance mass.
   @param varName - The name of the fit variable that has been parameterized.
   @return - A list of parameters corresponding to the variable.
*/
std::vector<TString> SigParam::listParamsForVar(TString varName) {
  std::vector<TString> result; result.clear();
  TString prefix[10] = {"a_","b_","c_","d_","e_","f_","g_","h_","i_","j_"};
  int varOrder = getNParamsForVar(varName);
  if (varOrder == 0) {
    result.push_back(varName);
  }
  else {
    for (int i_p = 0; i_p < varOrder; i_p++) {
      if (i_p == 10) {
	std::cout << "SigParam: Error! listParamsForVar() supports < 10 var."
		  << std::endl;
	exit(0);
      }
      result.push_back(Form("%s%s", (prefix[i_p]).Data(), varName.Data()));
    }
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Load the signal parameterization from file.
   @param directory - Name of the directory housing the input workspace.
   @return - True iff the file is successfully loaded.
*/
bool SigParam::loadParameterization(TString directory, TString signalType){
  std::cout << "SigParam: Load parameterization from " << directory
	    << " for signal type " << signalType << std::endl;
  setDirectory(directory);
  bool parameterizationExists = false;
  TFile inputFile(Form("%s/res_%sworkspace.root", m_directory.Data(), 
		       m_signalType.Data()));
  if (inputFile.IsOpen()) {
    m_ws = (RooWorkspace*)inputFile.Get("signalWS");
    if (m_ws) {
      parameterizationExists = true;
      setSignalType(signalType);
      m_currFunction = m_ws->var("functionName")->GetTitle();
      std::cout << "SigParam: Successfully loaded from file!" << std::endl;
    }
  }
  if (!parameterizationExists) {
    std::cout << "SigParam: ERROR loading from file." << std::endl;
    exit(0);
  }
  
  return parameterizationExists;
}

/**
   -----------------------------------------------------------------------------
   Parameterize the resonance shape in all categories.
   @param function - The functional form of the resonance.
   @return - True iff. all fits converge.
*/
bool SigParam::makeAllParameterizations(TString function) {
  std::cout << "SigParam: Engage full signal parameterization!" << std::endl;
  if (!functionIsDefined(function)) {
    std::cout << "SigParam: Error! Improper function " << function << std::endl;
    return false;
  }

  bool result = true;
  // Define models in each category independently:
  for (int i_c = 0; i_c < m_nCategories; i_c++) {
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
   @return - True iff. all fits converge.
*/
bool SigParam::makeCategoryParameterization(int cateIndex, TString function) {
  if (m_verbose) {
    std::cout << "SigParam: parameterizing category " << cateIndex << std::endl;
  }
  if (!functionIsDefined(function)) {
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
  
  // Define the RooFormulaVars which control mH parameterization:
  // Loop over mass points, define resonance model in each:
  for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
    TString currKey = getKey(currMassPoints[i_m],cateIndex);
    double mRegVal = regularizedMass(currMassPoints[i_m]);
    double mResVal = currMassPoints[i_m];
    parameterizeFunction(function, mRegVal, mResVal, cateIndex, false);
        
    // Create, individual resonance shapes, add to simultaneous PDF:
    resonanceCreator(currMassPoints[i_m], cateIndex, function);
    currSim->addPdf(*m_ws->pdf(Form("sigPdf_%s%s",m_signalType.Data(),
				    currKey.Data())), currKey);
  }
  
  // Import simultaneous PDF into workspace:
  m_ws->import(*currSim);
  
  // Set to store all mass observables and weight variable:
  RooArgSet *args = new RooArgSet();
    
  // Get the RooDataSets and observables in loop over categories:
  std::map<std::string,RooDataSet*> currDataMap; currDataMap.clear();
  for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
    // maybe add if statement here...
    TString currKey = getKey(currMassPoints[i_m], cateIndex);
    
    // Then add to the map:
    currDataMap[((std::string)currKey)]
      = (RooDataSet*)m_ws->data(Form("data_%s",currKey.Data()));
    
    // Reduce the dataset to within (rounded) +/-5 sigma of the mean
    double windowPercent = 0.1;
    double dataMin = round((1.0-windowPercent) * currMassPoints[i_m]);
    double dataMax = round((1.0+windowPercent) * currMassPoints[i_m]);
    //currDataMap[((std::string)currKey)]
    //->reduce(Form("m_yy_%s>%f && m_yy_%s<%f",
    //		    currKey.Data(), dataMin, currKey.Data(), dataMax));
    if (!m_verbose) {
      RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    }
    m_ws->var(Form("m_yy_%s",currKey.Data()))->setRange(dataMin, dataMax);
    
    // Add the current observable to the set:
    args->add(*m_ws->var(Form("m_yy_%s",currKey.Data())));
  }
  
  // Add weight variable to arg set:
  args->add(*(m_ws->var("wt")));
  
  // Create the combined dataset to be profiled:
  RooDataSet *obsData = new RooDataSet(Form("data_c%d",cateIndex),
				       Form("data_c%d",cateIndex), *args,
				       RooFit::Index(*m_cat),
				       RooFit::Import(currDataMap),
				       RooFit::WeightVar(*m_ws->var("wt")));
  m_ws->import(*obsData);
    
  // Then fit simultaneous PDF to combined dataset:
  RooFitResult *result = fitResult(cateIndex);
  
  // Save the test statistics:
  m_testStats[Form("NLL_c%d", cateIndex)] = getTestStatLatest("NLL");
  if (m_binned) {
    m_testStats[Form("Chi2_c%d", cateIndex)] = getTestStatLatest("Chi2");
  }
  
  // Then construct parametric resonance (function of mResonance):
  m_ws->factory("m_yy[10,10000]");
  parameterizeFunction(function, -999, -999, cateIndex, true);
  
  // Create the parameterized resonance shape (function of mResonance):
  resonanceCreator(-999, cateIndex, Form("%s_Parameterized",function.Data()));
  
  // Get the yield parameterization:
  makeYieldParameterization(cateIndex);
  
  // Return the fit status:
  // Bad fits have null pointer to status, nonzero value, or infinite/nan NLL
  return (result && result->status() == 0 && std::isfinite(m_currNLL));
}

/**
   -----------------------------------------------------------------------------
   Create the resonance for a single mass point and category. 
   @param resonanceMass - The truth mass of the resonance
   @param cateIndex - The index of the category to fit.
   @param function - The functional form of the resonance.
   @return - True iff. all fits converge.
*/
bool SigParam::makeSingleResonance(double resonanceMass, int cateIndex,
				   TString function) {
  if (m_verbose) {
    std::cout << "SigParam: individual fit in category " << cateIndex
	    << " and at mass " << resonanceMass << std::endl;
  }
  if (!functionIsDefined(function)) {
    std::cout << "SigParam: Error! Improper function " << function << std::endl;
    return false;
  }
  
  // Before calling resonanceCreator, need to define dependent variables.
  TString currKey = getKey(resonanceMass, cateIndex);
  
  // Get a list of variables for the current PDF and add them to the workspace:
  std::vector<TString> currVars = variablesForFunction(function);
  for (int i_v = 0; i_v < (int)currVars.size(); i_v++) {
    TString varSuffix = "";
    if (currVars[i_v].Contains("nCB") || currVars[i_v].Contains("frac")) {
      varSuffix = Form("%sc%d", m_signalType.Data(), cateIndex);
    }
    else {
      varSuffix = Form("%s%s", m_signalType.Data(), currKey.Data());
    }
    TString initialValAndRange = getParamState(currVars[i_v]);
    m_ws->factory(Form("%s_%s%s", currVars[i_v].Data(), varSuffix.Data(), 
		       initialValAndRange.Data()));
    // If it is a mean value, set to the resonance mass initially:
    if ((currVars[i_v].Contains("muCB") || currVars[i_v].Contains("muGA") ||
	 currVars[i_v].Contains("muVoigt") || currVars[i_v].Contains("muLA")) &&
	!(m_ws->var(Form("%s_%s", currVars[i_v].Data(), varSuffix.Data()))
	  ->isConstant())) {
      if (m_verbose) std::cout << "Setting mean." << std::endl;
      m_ws->var(Form("%s_%s", currVars[i_v].Data(), varSuffix.Data()))
	->setVal(resonanceMass);
      m_ws->var(Form("%s_%s", currVars[i_v].Data(), varSuffix.Data()))
	->setRange(0.9*resonanceMass,1.1*resonanceMass);
    }
  }
  
  // Create the single resonance PDF:
  resonanceCreator(resonanceMass, cateIndex, function);
  
  // Reduce the dataset to within (rounded) +/-5 sigma of the mean
  TString dataName = Form("data_%s",currKey.Data());
  double windowPercent = 0.1;
  double dataMin = round((1.0-windowPercent) * resonanceMass);
  double dataMax = round((1.0+windowPercent) * resonanceMass);
  //m_ws->data(dataName)->reduce(Form("m_yy_%s>%f && m_yy_%s<%f",currKey.Data(),
  //				    dataMin, currKey.Data(), dataMax));
  
  if (!m_verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  }
  m_ws->var(Form("m_yy_%s",currKey.Data()))->setRange(dataMin, dataMax);
  
  // Then fit single PDF to single dataset:
  RooFitResult *result = fitResult(resonanceMass, cateIndex);
  
  // Save the test statistics:
  m_testStats[Form("NLL_%s", currKey.Data())] = getTestStatLatest("NLL");
  if (m_binned) {
    m_testStats[Form("Chi2_%s", currKey.Data())] = getTestStatLatest("Chi2");
  }
  
  // Return the fit status. 
  // Bad fits have null pointer to status, nonzero value, or infinite/nan NLL
  return (result && result->status() == 0 && std::isfinite(m_currNLL));
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
  double mResErrors[100] = {0.0};
  double yieldErrors[100] = {0.0};
  
  // Retrieve dataset yield for each mass point that was imported:
  std::vector<double> currMassPoints = massPointsForCategory(cateIndex);
  for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
    TString dataName
      = Form("data_%s",(getKey(currMassPoints[i_m],cateIndex)).Data());
    if ((m_ws->data(dataName))) {
      mResValues[nResPoints] = regularizedMass(currMassPoints[i_m]);
      yieldValues[nResPoints] = (*m_ws->data(dataName)).sumEntries();
      yieldErrors[nResPoints] = normalizationError(m_ws->data(dataName));
      nResPoints++;
    }
  }
  
  // Use TF1 and TGraph to fit the yield:
  m_yieldFunc[cateIndex] = new TF1("yieldFunc", "pol3", 0.0, 0.5);
  m_yieldGraph[cateIndex] = new TGraphErrors(nResPoints, mResValues, 
					     yieldValues, mResErrors,
					     yieldErrors);
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
   @return - The floating value of the mass in GeV.
*/
double SigParam::massIntToDouble(int massInteger) {
  return ((double)massInteger) / 1000.0;
}

/**
   -----------------------------------------------------------------------------
   Convert the resonance mass value to an integer representation. 
   @param resonanceMass - The value of the mass.
   @return - The integer representation of the mass.
*/
int SigParam::massDoubleToInt(double resonanceMass) {
  resonanceMass = round(resonanceMass);
  return (int)(resonanceMass * 1000);
}

/**
   -----------------------------------------------------------------------------
   Get a list of mass points corresponding to a single category.
   @param cateIndex - The index of the category.
   @return - A vector of mass values.
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
   Import a list of names for the fit categories. Otherwise, categories will 
   simply be numbered. Sets the value of private variable m_cateNames. 
   @param cateNames - A vector of category names.
*/
void SigParam::nameTheCategories(std::vector<TString> cateNames) {
  if ((int)cateNames.size() != m_nCategories && m_verbose) {
    std::cout << "SigParam: WARNING! Number of defined categories different from number of category names provided." << std::endl;
  }
  m_cateNames = cateNames;
}

/**
   -----------------------------------------------------------------------------
   Get the normalization error for a weighted dataset.
   @param dataSet - The dataset for which we want the total normalization error.
   @return - The sumw2 normalization error for the dataset.
*/
double SigParam::normalizationError(RooAbsData *dataSet) {
  double normError = 0.0;
  // Loop over events stored in the input dataset:
  for (int i_d = 0; i_d < dataSet->numEntries(); i_d++) {
    dataSet->get(i_d);
    normError += dataSet->weightSquared();
  }
  return sqrt(normError);
}

/**
   -----------------------------------------------------------------------------
   Parameterize the specified function.
   @param function - The functional form of the resonance.
   @param mRegularized - The regularized mass.
   @param mResonance - The resonance mass.
   @param cateIndex - The event category index (starting at 0).
   @param parameterized - True iff for floating mResonance. 
*/
void SigParam::parameterizeFunction(TString function, double mRegularized,
				    double mResonance, int cateIndex,
				    bool parameterized) {
  std::vector<TString> currVariables = variablesForFunction(function);
  for (int i_v = 0; i_v < (int)currVariables.size(); i_v++) {
    parameterizeVar(currVariables[i_v], mRegularized, mResonance, cateIndex,
		    parameterized);
  }
}

/**
   -----------------------------------------------------------------------------
   Parameterize the specified variable.
   @param varName - The name of the shape variable of interest.
   @param mRegularized - The regularized mass.
   @param mResonance - The resonance mass.
   @param cateIndex - The event category index (starting at 0).
   @param parameterized - True iff for floating mResonance. 
*/
void SigParam::parameterizeVar(TString varName, double mRegularized,
			       double mResonance, int cateIndex,
			       bool parameterized) {
  if (m_verbose) {
    std::cout << "SigParam: parameterizeVar(" << varName << ", " << cateIndex 
	      << ")" << std::endl;
  }
  
  // Make sure all parameterization variables have been added:
  addVariable(varName, cateIndex);
  TString currKey = getKey(mResonance, cateIndex);
  
  // Then create a string of the parameters from the vector of parameters:
  TString paramString = "";
  std::vector<TString> paramList = listParamsForVar(varName);
  for (int i_p = 0; i_p < (int)paramList.size(); i_p++) {
    paramString
      += Form("%s_%sc%d",(paramList[i_p]).Data(),m_signalType.Data(),cateIndex);
    if (i_p < (int)paramList.size()-1) paramString += ",";
  }
  
  // Get the functional form of the variable:
  TString functionForm = getVarParameterization(varName);
  if (!functionForm.EqualTo("")) {
    if (parameterized) {
      int nParams = (int)paramList.size();
      functionForm = functionForm.ReplaceAll("obs", Form("@%d",nParams));
      functionForm = functionForm.ReplaceAll("mRes", Form("@%d",nParams+1));
      paramString += ",mRegularized,mResonance";
      m_ws->factory(Form("expr::%s_%sc%d('%s',{%s})", varName.Data(), 
			 m_signalType.Data(), cateIndex, functionForm.Data(), 
			 paramString.Data()));
    }
    else {
      functionForm = functionForm.ReplaceAll("obs", Form("%f", mRegularized));
      functionForm = functionForm.ReplaceAll("mRes", Form("%f", mResonance));
      m_ws->factory(Form("expr::%s_%s%s('%s',{%s})", varName.Data(), 
			 m_signalType.Data(),currKey.Data(),functionForm.Data(),
			 paramString.Data()));
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Create a RooDataSet with a new RooRealVar. 
*/
RooDataSet* SigParam::plotData(RooAbsData *data, RooRealVar *observable,
			       double xBins, double resonanceMass) {
  // Clone the dataset:
  TString cloneName = data->GetName();
  RooDataSet *dataClone 
    = (RooDataSet*)data->Clone(Form("%s_copy",cloneName.Data()));
  dataClone->changeObservableName(observable->GetName(), "m_yy");
  m_ws->import(*dataClone);
  return dataClone;
}

/**
   -----------------------------------------------------------------------------
   Create a ratio plot (or subtraction plot, for the moment...)
   @param data - The RooAbsData set for comparison.
   @param pdf - The PDF for comparison.
   @param observable - The mass observable for data and pdf. 
   @param xBins - The number of bins for the observable.
   @param chi2Prob - The chi^2 probability of the fit at this point.
   @return - A TGraphErrors to plot.
*/
TGraphErrors* SigParam::plotSubtraction(RooAbsData *data, RooAbsPdf *pdf, 
					RooRealVar *observable, double xBins,
					double &chi2Prob){
  double minOrigin = observable->getMin();
  double maxOrigin = observable->getMax();
  double nEvents = data->sumEntries();
  TH1F *originHist
    = (TH1F*)data->createHistogram("dataSub", *observable,
				   RooFit::Binning(xBins,minOrigin,maxOrigin));
  TGraphErrors *result = new TGraphErrors();
  double increment = (maxOrigin - minOrigin) / ((double)xBins);
  
  RooAbsReal* intTot
    = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable),
				       RooFit::NormSet(*observable), 
				       RooFit::Range("fullRange"));
  double valTot = intTot->getVal();
  
  int pointIndex = 0;
  for (double i_m = minOrigin; i_m < maxOrigin; i_m += increment) {
    if (!m_verbose) {
      RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    }
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
    double currWeight = currDataWeight - currPdfWeight;
    result->SetPoint(pointIndex, currMass, currWeight);
    
    double currError = originHist->GetBinError(pointIndex+1);
    result->SetPointError(pointIndex, 0.0, currError);
    pointIndex++;
  }
  observable->setMin(minOrigin);
  observable->setMax(maxOrigin);
  delete originHist;
  return result;
}

/**
   -----------------------------------------------------------------------------
   Create a ratio plot (or subtraction plot, for the moment...)
   @param data - The RooAbsData set for comparison.
   @param pdf - The PDF for comparison.
   @param observable - The mass observable for data and pdf. 
   @param xBins - The number of bins for the observable.
   @param chi2Val - The chi^2 probability of the fit at this point only.
   @return - A TGraphErrors to plot.
*/
TGraphErrors* SigParam::plotDivision(RooAbsData *data, RooAbsPdf *pdf, 
				     RooRealVar *observable, double xBins,
				     double &chi2Prob){
  // Store the original variable range:
  double minOrigin = observable->getMin();
  double maxOrigin = observable->getMax();
  
  double nEvents = data->sumEntries();
  TH1F *originHist
    = (TH1F*)data->createHistogram("dataSub", *observable,
				   RooFit::Binning(xBins,minOrigin,maxOrigin));
  TGraphErrors *result = new TGraphErrors();
  double increment = ((maxOrigin - minOrigin) / ((double)xBins));
  
  RooAbsReal* intTot
    = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable),
				       RooFit::NormSet(*observable), 
				       RooFit::Range("fullRange"));
  double valTot = intTot->getVal();
  int pointIndex = 0; int pointIndexNonZero = 0;
  for (double i_m = minOrigin; i_m < maxOrigin; i_m += increment) {
    if (!m_verbose) {
      RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    }
    observable->setRange(Form("range%2.2f",i_m), i_m, (i_m+increment));
    RooAbsReal* intCurr
      = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable), 
					 RooFit::NormSet(*observable), 
					 RooFit::Range(Form("range%2.2f",i_m)));
    double valCurr = intCurr->getVal();
    
    double currMass = i_m + (0.5*increment);
    double currPdfWeight = nEvents * (valCurr / valTot);
    TString varName = observable->GetName();
    double currDataWeight
      = data->sumEntries(Form("%s>%f&&%s<%f",varName.Data(),
			      i_m,varName.Data(), (i_m+increment)));
    double currWeight = currDataWeight / currPdfWeight;
    result->SetPoint(pointIndex, currMass, currWeight);
    
    double currError = originHist->GetBinError(pointIndex+1) / currPdfWeight;
    result->SetPointError(pointIndex, 0.0, currError);
    pointIndex++;
    
    double currChi2 = (((currDataWeight-currPdfWeight) * 
    			(currDataWeight-currPdfWeight)) / 
    		       ((originHist->GetBinError(pointIndex+1)) * 
    			(originHist->GetBinError(pointIndex+1))));
    if (std::isfinite(currChi2)) {
      chi2Prob += currChi2;
      pointIndexNonZero++;
    }
  }
  // Return to the original variable range, store the chi^2 value:
  observable->setMin(minOrigin);
  observable->setMax(maxOrigin);
  chi2Prob = TMath::Prob(chi2Prob, pointIndexNonZero);
  delete originHist;
  return result;
}

/**
   -----------------------------------------------------------------------------
   Plot a resonance PDF for all masses defined for one category. Note: this 
   method requires a simultaneous fit to run first, so that the m_yy object
   can be created in the workspace. 
   @param cateIndex - The index of the category.
*/
void SigParam::plotCategoryResonances(int cateIndex) {
  if (m_verbose) {
    std::cout << "SigParam: Plot resonances in cate. " << cateIndex 
	      << std::endl;
  }
  
  // Create a canvas with two pads (one main plot, one subtraction plot)
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
  
  // Also get the range of the plot:
  int xMin = 0; int xMax = 0;
  for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
    RooRealVar *currObs = m_ws
      ->var(Form("m_yy_%s",(getKey(currMassPoints[i_m],cateIndex).Data())));
    if (i_m == 0) { xMin = currObs->getMin(); xMax = currObs->getMax(); }
    else if (currObs->getMin() < xMin) xMin = currObs->getMin();
    else if (currObs->getMax() > xMax) xMax = currObs->getMax();
  }
  int xBins = (int)(m_nBinsPerGeV * (xMax - xMin));
  // Set this new limited range for the RooRealVar observable:
  if (!m_verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  }
  m_ws->var("m_yy")->setRange(xMin, xMax);
  
  // Create a RooPlot for the fit and data:
  RooPlot* frame
    = m_ws->var("m_yy")->frame(RooFit::Bins(xBins), RooFit::Range(xMin, xMax));
  frame->SetYTitle(Form("Events/%2.2f GeV", (1.0/((double)m_nBinsPerGeV))));
  frame->SetXTitle("Mass [GeV]");
  
  TH1F *medianHist = NULL;
  
  // Loop over mass points, drawing data and PDF for each:
  for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
    pad1->cd();
    RooAbsData *currData = NULL;
    RooAbsPdf *currPdf = NULL;
    RooRealVar *currObs = NULL;
    TString currKey = getKey(currMassPoints[i_m], cateIndex);
    
    // Plot the data:
    if ((m_ws->data(Form("data_%s",currKey.Data())))) {
      currData = (m_ws->data(Form("data_%s",currKey.Data())));
      currObs = m_ws->var(Form("m_yy_%s",currKey.Data()));
      currData = plotData(currData, currObs, xBins, currMassPoints[i_m]);
      currData->plotOn(frame);
    }
    else {
      std::cout << "SigParam: data for plotting undefined." << std::endl;
      return;
    }
    
    // Then plot the parameterized signal PDF:
    if ((m_ws->pdf(Form("sigPdf_%sc%d",m_signalType.Data(),cateIndex)))) {
      (*m_ws->var("mResonance")).setVal(currMassPoints[i_m]);
      currPdf = (m_ws->pdf(Form("sigPdf_%sc%d",m_signalType.Data(),cateIndex)));
      currPdf->plotOn(frame, RooFit::LineColor(4));
    }
    else {
      std::cout << "SigParam: ERROR! Parameterized shape not defined!" 
		<< std::endl;
      exit(0);
    }
    
    // Then draw the frame:
    if (i_m == 0) {
      frame->Draw();
      
      // Print ATLAS text on the plot:    
      TLatex l; l.SetNDC(); l.SetTextColor(kBlack);
      l.SetTextFont(72); l.SetTextSize(0.05); l.DrawLatex(0.7,0.88,"ATLAS");
      //l.SetTextFont(42); l.SetTextSize(0.05); l.DrawLatex(0.82,0.88,"Simulation");
      l.SetTextFont(42); l.SetTextSize(0.05); l.DrawLatex(0.82,0.88,"Internal");
      //l.DrawLatex(0.7, 0.81, Form("#scale[0.8]{#sqrt{s} = 13 TeV: #scale[0.7]{#int}Ldt = %2.1f fb^{-1}}",analysis_luminosity));
      l.DrawLatex(0.7, 0.82, "#sqrt{s} = 13 TeV");
      
      TLatex text; text.SetNDC(); text.SetTextColor(1);
      if ((int)m_cateNames.size() == m_nCategories) {
	text.DrawLatex(0.7, 0.76, m_cateNames[cateIndex]);
      }
      else {
	text.DrawLatex(0.7, 0.76, Form("category %d", cateIndex));
      }
    }
    else frame->Draw("SAME");

    // Switch to sub-plot:
    pad2->cd();
    if (i_m == 0) {
      medianHist = new TH1F("median", "median", xBins, xMin, xMax);
      for (int i_b = 1; i_b <= xBins; i_b++) medianHist->SetBinContent(i_b,1.0);
      medianHist->SetLineColor(kRed);
      medianHist->SetLineWidth(2);
      medianHist->GetXaxis()->SetTitle("Mass [GeV]");
      if (m_doRatioPlot) {
	medianHist->GetYaxis()->SetTitle("Data / Fit");
      }
      else medianHist->GetYaxis()->SetTitle("Data - Fit");
      medianHist->GetXaxis()->SetTitleOffset(0.95);
      medianHist->GetYaxis()->SetTitleOffset(0.7);
      medianHist->GetXaxis()->SetTitleSize(0.1);
      medianHist->GetYaxis()->SetTitleSize(0.1);
      medianHist->GetXaxis()->SetLabelSize(0.1);
      medianHist->GetYaxis()->SetLabelSize(0.1);
      medianHist->GetYaxis()->SetRangeUser(-0.2, 2.2);
      medianHist->GetYaxis()->SetNdivisions(5);
      medianHist->Draw();
      
      TLine *line = new TLine();
      line->SetLineStyle(1);
      line->SetLineWidth(2);
      line->SetLineColor(kRed);
      if (m_doRatioPlot) {
	line->SetLineWidth(1);
	line->SetLineStyle(2);
	line->DrawLine(xMin,((1.0+m_ratioMin)/2.0),xMax,((1.0+m_ratioMin)/2.0));
	line->DrawLine(xMin,((1.0+m_ratioMax)/2.0),xMax,((1.0+m_ratioMax)/2.0));
      }
      else line->DrawLine(xMin, 0.0, xMax, 0.0);
    }
    
    double currChi2Prob = 0.0;
    TGraphErrors* subData = (m_doRatioPlot) ? 
      plotDivision(currData,currPdf,m_ws->var("m_yy"),xBins,currChi2Prob) :
      plotSubtraction(currData,currPdf,m_ws->var("m_yy"),xBins,currChi2Prob);
    subData->Draw("EPSAME");
  }
  // Print the canvas:
  can->Print(Form("%s/plot_paramResonance_c%d%s", 
		  m_directory.Data(), cateIndex, m_fileFormat.Data()));
  delete can;
  delete medianHist;
}

/**
   -----------------------------------------------------------------------------
   Plot a resonance PDF for one value of the resonance mass in one category.
   @param resonanceMass - The mass value in GeV.
   @param cateIndex - The index of the category.
   @param dataType - "asimov", "mctoy", "pdftoy".
*/
void SigParam::plotSingleResonance(double resonanceMass, int cateIndex,
				   TString dataType) {
  if (m_verbose) {
    std::cout << "SigParam: Plotting resonance at mass " << resonanceMass 
	      << " in category " << cateIndex << " with dataset " 
	      << dataType << std::endl;
  }
  
  TString currKey = getKey(resonanceMass,cateIndex);
  bool parameterized = false;
  
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

  TString pdfName = "";
  TString dataName = "";
  TString obsName = "";
  // For fits to Asimov, MC toy, or PDF toy data:
  if (dataType.Contains("asimov") || dataType.Contains("mctoy") ||
      dataType.Contains("pdftoy")) {
    dataName = Form("data_%s_%s", dataType.Data(), currKey.Data());
    
    // Loof for parameterized model and then single mass point model:
    pdfName = Form("sigPdf_%sc%d", m_signalType.Data(), cateIndex);
    if (m_ws->pdf(pdfName)) {
      setResMassConstant(true, resonanceMass);
      obsName = "m_yy";
      parameterized = true;
    }
    else {
      pdfName = Form("sigPdf_%s%s", m_signalType.Data(), currKey.Data());
      obsName = Form("m_yy_%s", currKey.Data());
    }
  }
  
  // Else for nominal fits:
  else {
    dataName = Form("data_%s", currKey.Data());
    obsName = Form("m_yy_%s", currKey.Data());
    // Look for parameterized model and then single mass point model:
    pdfName = Form("sigPdf_%sc%d", m_signalType.Data(), cateIndex);
    if (m_ws->pdf(pdfName)) {
      setResMassConstant(true, resonanceMass);
      // Clone the dataset with a different observable:
      int rBins = (int)(m_nBinsPerGeV*(m_ws->var(obsName)->getMax() - 
				       m_ws->var(obsName)->getMin()));
      plotData(m_ws->data(dataName), m_ws->var(obsName), rBins, resonanceMass);
      dataName = Form("%s_copy", dataName.Data());
      parameterized = true;
    }
    else pdfName = Form("sigPdf_%s%s", m_signalType.Data(), currKey.Data());
  }
  
  // Then check that both data and PDF exist:
  if (!m_ws->data(dataName) || !m_ws->pdf(pdfName)) {
    std::cout << "SigParam: data or PDF for plotting undefined." << std::endl;
    exit(0);
  }
    
  // Load the appropriate observable and define the plot binning:
  double rMin = m_ws->var(obsName)->getMin();
  double rMax = m_ws->var(obsName)->getMax();
  int rBins = (int)(m_nBinsPerGeV*(rMax - rMin));
  
  // Create the RooPlot object using the observable, and set axis titles:
  RooPlot *frame 
    = m_ws->var(obsName)->frame(RooFit::Bins(rBins),RooFit::Range(rMin,rMax));
  frame->SetYTitle(Form("Events/%2.2f GeV", (1.0/((double)m_nBinsPerGeV))));
  frame->SetXTitle("Mass [GeV]");
    
  // Then add the data and PDF to the RooPlot:
  m_ws->data(dataName)->plotOn(frame);
  m_ws->pdf(pdfName)->plotOn(frame, RooFit::LineColor(2));
    
  // Draw the RooPlot:
  frame->Draw();
  // Special y-axis ranges for log-scale plots:
  if (m_useLogYAxis) {
    gPad->SetLogy();
    frame->GetYaxis()
      ->SetRangeUser(0.00001 *(*m_ws->data(dataName)).sumEntries(), 
		     (*m_ws->data(dataName)).sumEntries());
  }
    
  TLatex l; l.SetNDC(); l.SetTextColor(kBlack);
  l.SetTextFont(72); l.SetTextSize(0.05); l.DrawLatex(0.20,0.88,"ATLAS");
  //l.SetTextFont(42); l.SetTextSize(0.05); l.DrawLatex(0.32,0.88,"Simulation");
  l.SetTextFont(42); l.SetTextSize(0.05); l.DrawLatex(0.32,0.88,"Internal");
  //l.DrawLatex(0.2, 0.81, Form("#scale[0.8]{#sqrt{s} = 13 TeV: #scale[0.7]{#int}Ldt = %2.1f fb^{-1}}",analysis_luminosity));
  l.DrawLatex(0.2, 0.82, "#sqrt{s} = 13 TeV");
  
  TLatex text; text.SetNDC(); text.SetTextColor(1);
  if ((int)m_cateNames.size() == m_nCategories && m_nCategories > 0) {
    text.DrawLatex(0.2, 0.76, m_cateNames[cateIndex]);
  }
  else {
    text.DrawLatex(0.2, 0.76, Form("category %d", cateIndex));
  }
  
  double yVal = 0.88;
  TLatex varText; varText.SetNDC(); varText.SetTextColor(1);
  varText.SetTextSize(0.04);
  std::vector<TString> currVars = variablesForFunction(m_currFunction);
  for (int i_v = 0; i_v < (int)currVars.size(); i_v++) {
    TString currName = currVars[i_v];
    double currVal = 0.0;
    if (parameterized) getParameterizedValue(currName,resonanceMass,cateIndex);
    else getParameterValue(currName, resonanceMass, cateIndex);
    currName.ReplaceAll("frac", "fraction_{");
    currName.ReplaceAll("Nom","");
    currName.ReplaceAll("nCB","n_{CB");
    currName.ReplaceAll("width","w_{");
    currName.ReplaceAll("sigma","#sigma_{");
    currName.ReplaceAll("alpha","#alpha_{");
    currName.ReplaceAll("mu", "#mu_{");
    currName += "}";
    varText.DrawLatex(0.7, yVal, Form("%s\t = %2.2f",currName.Data(),currVal));
    yVal -= 0.06;
  }
    
  // Move to second pad for ratio or subtraction plot:
  pad2->cd();
  
  TH1F *medianHist = new TH1F("median", "median", rBins, rMin, rMax);
  medianHist->SetLineColor(kRed);
  medianHist->SetLineWidth(2);
  medianHist->GetXaxis()->SetTitle("Mass [GeV]");
  if (m_doRatioPlot) {
    for (int i_b = 1; i_b <= rBins; i_b++) medianHist->SetBinContent(i_b,1.0);
    medianHist->GetYaxis()->SetTitle("Data / Fit");
    medianHist->GetYaxis()->SetRangeUser(-0.2,2.2);
  }
  else {
    for (int i_b = 1; i_b <= rBins; i_b++) medianHist->SetBinContent(i_b,0.0);
    medianHist->GetYaxis()->SetTitle("Data - Fit");
  }
  medianHist->GetYaxis()->SetNdivisions(5);
  medianHist->GetXaxis()->SetTitleOffset(0.95);
  medianHist->GetYaxis()->SetTitleOffset(0.7);
  medianHist->GetXaxis()->SetTitleSize(0.1);
  medianHist->GetYaxis()->SetTitleSize(0.1);
  medianHist->GetXaxis()->SetLabelSize(0.1);
  medianHist->GetYaxis()->SetLabelSize(0.1);
  medianHist->Draw();
  
  double currChi2Prob = 0.0;
  TGraphErrors* subData = (m_doRatioPlot) ?
    plotDivision(m_ws->data(dataName), m_ws->pdf(pdfName), m_ws->var(obsName), 
		 rBins, currChi2Prob) : 
    plotSubtraction(m_ws->data(dataName), m_ws->pdf(pdfName), 
		    m_ws->var(obsName), rBins, currChi2Prob);
  subData->GetXaxis()->SetTitle("Mass [GeV]");
  
  // Draw lines showing ratio of 1.0 +/-0.5
  TLine *line = new TLine();
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  if (m_doRatioPlot) {
    line->SetLineWidth(1);
    line->SetLineStyle(2);
    line->DrawLine(rMin,((1.0+m_ratioMin)/2.0),rMax,((1.0+m_ratioMin)/2.0));
    line->DrawLine(rMin,((1.0+m_ratioMax)/2.0),rMax,((1.0+m_ratioMax)/2.0));
  }
  else line->DrawLine(rMin, 0.0, rMax, 0.0);
  subData->Draw("EPSAME");
  TLatex chiText; chiText.SetNDC(); chiText.SetTextColor(1);
  chiText.SetTextSize(0.1);
  //chiText.DrawLatex(0.7, 0.9, Form("p_{#chi} = %2.2f", currChi2Prob));
    
  can->Print(Form("%s/plot%s_singleRes_m%2.2f_c%d%s", m_directory.Data(), 
		  dataType.Data(), resonanceMass, cateIndex,
		  m_fileFormat.Data()));
  delete can;
  delete medianHist;
}

/**
   -----------------------------------------------------------------------------
   Graph the signal yields as function of resonance mass in the given category.
   @param cateIndex - The index of the category for yield plotting. 
*/
void SigParam::plotYields(int cateIndex) {
  std::cout << "SigParam: Plotting yields in category " << cateIndex
	    << std::endl;
  TCanvas *can = new TCanvas("can","can",800,800);
  can->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0.00, 0.33, 1.00, 1.00);
  TPad *pad2 = new TPad("pad2", "pad2", 0.00, 0.00, 1.00, 0.33);
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
  m_yieldGraph[cateIndex]->GetYaxis()->SetTitle("Signal yield");
  m_yieldGraph[cateIndex]->Draw("AEP");
  m_yieldFunc[cateIndex]->Draw("LSAME");
  TString formattedSignal = m_signalType.Data();
  formattedSignal.ReplaceAll("_","");

  TLatex text; text.SetNDC(); text.SetTextColor(1);
  text.DrawLatex(0.7, 0.76, Form("%s signal",formattedSignal.Data()));
  if ((int)m_cateNames.size() == m_nCategories) {
    text.DrawLatex(0.7, 0.70, m_cateNames[cateIndex]);
  }
  else {
    text.DrawLatex(0.7, 0.70, Form("Category %d", cateIndex));
  }
  
  // Print ATLAS text on the plot:    
  TLatex l; l.SetNDC(); l.SetTextColor(kBlack);
  l.SetTextFont(72); l.SetTextSize(0.05); l.DrawLatex(0.7,0.88,"ATLAS");
  //l.SetTextFont(42); l.SetTextSize(0.05); l.DrawLatex(0.82,0.88,"Simulation");
  l.SetTextFont(42); l.SetTextSize(0.05); l.DrawLatex(0.82,0.88,"Internal");
  //l.DrawLatex(0.7, 0.81, Form("#scale[0.8]{#sqrt{s} = 13 TeV: #scale[0.7]{#int}Ldt = %2.1f fb^{-1}}",analysis_luminosity));
  l.DrawLatex(0.7, 0.82, "#sqrt{s} = 13 TeV");
  
  pad2->cd();
  
  TGraphErrors *gRatio = new TGraphErrors();
  gRatio->SetNameTitle(Form("ratio_c%d",cateIndex),Form("ratio_c%d",cateIndex));
  for (int i_p = 0; i_p < m_yieldGraph[cateIndex]->GetN(); i_p++) {
    double currX = 0.0; double currY = 0.0;
    m_yieldGraph[cateIndex]->GetPoint(i_p, currX, currY);
    double errorX = m_yieldGraph[cateIndex]->GetErrorX(i_p);
    double errorY = m_yieldGraph[cateIndex]->GetErrorY(i_p);
    gRatio->SetPoint(i_p, currX, (currY/m_yieldFunc[cateIndex]->Eval(currX)));
    gRatio->SetPointError(i_p, (errorX / m_yieldFunc[cateIndex]->Eval(currX)),
			  (errorY / m_yieldFunc[cateIndex]->Eval(currX)));
  }
  gRatio->GetXaxis()->SetTitle("(Mass-100)/100 [GeV]");
  gRatio->GetYaxis()->SetTitle("Data / Fit");
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
  
  TLine *line = new TLine();
  line->SetLineStyle(1);
  line->SetLineWidth(1);
  line->SetLineColor(kBlue);
  line->DrawLine(gRatio->GetXaxis()->GetXmin(), 1,
		 gRatio->GetXaxis()->GetXmax(), 1);
  
  can->Print(Form("%s/plot_paramYield_%sc%d%s", m_directory.Data(),
		  m_signalType.Data(), cateIndex, m_fileFormat.Data()));
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Prints the mean and standard deviation for the signal in each category. 
   NOTE: The standard deviation is the actual interval around the peak 
   containing 68.2% of the signal events, and not the shape width parameter.
   @param resonanceMass - The value of the resonance mass for which the 
   parameter values should be printed. 
*/
void SigParam::printResTable(double resonanceMass) {
  TString tableLocation = Form("%s/latexTable.txt", m_directory.Data());
  std::ofstream latexTable(tableLocation);
  latexTable << "\\begin{table}[!htb]" << std::endl;
  latexTable << "\\caption{Mean and resolution for the PDF in each category.}"
	     << std::endl;
  latexTable << "\\label{tab:fitMeanAndRes}" << std::endl;
  latexTable << "\\centering" << std::endl;
  latexTable << "\\begin{tabular}{lcc}" << std::endl;
  latexTable << "\\hline" << std::endl;
  latexTable << "Category & Mean [GeV] & Resolution [GeV] \\\\" 
	     << std::endl;
  latexTable << "\\hline" << std::endl;
  latexTable << "\\hline" << std::endl;
  // loop over categories, print mean and resoultion for each. 
  for (int i_c = 0; i_c < m_nCategories; i_c++) {
    double currMean = getMeanOrStdDev("Mean", resonanceMass, i_c);
    double currSigma = getMeanOrStdDev("StdDev", resonanceMass, i_c);
    if ((int)m_cateNames.size() == m_nCategories) {
      latexTable << m_cateNames[i_c] << " & " << currMean << " & " << currSigma
		 << " \\\\" << std::endl;
    }
    else {
      latexTable << "category " << i_c << " & " << currMean << " & " 
		 << currSigma << " \\\\" << std::endl;
    }
  }
  latexTable << "\\hline" << std::endl;
  latexTable << "\\end{tabular}" << std::endl;
  latexTable << "\\end{table}" << std::endl;
  latexTable.close();
  std::cout << "SigParam: Printed LaTex table: " << tableLocation << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Create a regularized mass variable that helps minimizers converge.
   @param resonanceMass - The mass value in GeV.
   @return - The regularized mass mR = (m-100)/100;
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
  if (m_verbose) {
    std::cout << "SigParam: Adding resonance to workspace" 
	      << "\tresonanceMass: " << resonanceMass << "\tcategory: " 
	      << cateIndex << std::endl;
  }
  m_currFunction = function;
  m_ws->var("functionName")->SetTitle(m_currFunction);

  // Check that the dataset has been defined and is not empty.
  if (!dataExists(resonanceMass, cateIndex) &&
      !function.Contains("Parameterized")) {
    std::cout << "SigParam: Cannot fit resonance, no dataset." << std::endl;
    exit(0);
  }
  
  TString currKey = getKey(resonanceMass, cateIndex);
  TString suffix = (function.Contains("Parameterized")) ? 
    Form("%sc%d", m_signalType.Data(), cateIndex) :
    Form("%s%s", m_signalType.Data(), currKey.Data());
  TString obsName = (function.Contains("Parameterized")) ? 
    "m_yy" : Form("m_yy_%s",currKey.Data());
  
  // Define the Crystal Ball + Gaussian shape:
  if (function.Contains("CBGA")) {
    // Cystal Ball component:
    m_ws->factory(Form("RooCBShape::pdfCB_%s(%s, prod::muCB_%s(muCBNom_%s%s), prod::sigmaCB_%s(sigmaCBNom_%s%s), alphaCB_%s, nCB_%sc%d)", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data(), suffix.Data(), m_signalType.Data(), cateIndex));
    // Gaussian component:
    if (m_sameCBGAMean) {
      m_ws->factory(Form("RooGaussian::pdfGA_%s(%s, prod::muGA_%s(muCBNom_%s%s), prod::sigmaGA_%s(sigmaGANom_%s%s))", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
    }
    else {
      m_ws->factory(Form("RooGaussian::pdfGA_%s(%s, prod::muGA_%s(muGANom_%s%s), prod::sigmaGA_%s(sigmaGANom_%s%s))", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
    }
    
    // Crystal Ball + Gaussian:
    m_ws->factory(Form("SUM::sigPdf_%s(fracCB_%sc%d*pdfCB_%s,pdfGA_%s)", suffix.Data(), m_signalType.Data(), cateIndex, suffix.Data(), suffix.Data()));
  }
  
  // Define double-sided Crystal Ball shape:
  else if (function.Contains("DoubleCB")) {
    m_ws->factory(Form("HggTwoSidedCBPdf::sigPdf_%s(%s, prod::muCB_%s(muCBNom_%s%s), prod::sigmaCB_%s(sigmaCBNom_%s%s), alphaCBLo_%s, nCBLo_%sc%d, alphaCBHi_%s, nCBHi_%sc%d)", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data(), suffix.Data(), m_signalType.Data(), cateIndex, suffix.Data(), m_signalType.Data(), cateIndex));
  }
  else if (function.Contains("GAx3")) {
    m_ws->factory(Form("RooGaussian::pdfGA1_%s(%s, prod::muGA1_%s(muGA1Nom_%s%s), prod::sigmaGA1_%s(sigmaGA1Nom_%s%s))", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
    m_ws->factory(Form("RooGaussian::pdfGA2_%s(%s, prod::muGA2_%s(muGA2Nom_%s%s), prod::sigmaGA2_%s(sigmaGA2Nom_%s%s))", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
    m_ws->factory(Form("RooGaussian::pdfGA3_%s(%s, prod::muGA3_%s(muGA3Nom_%s%s), prod::sigmaGA3_%s(sigmaGA3Nom_%s%s))", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
  
    m_ws->factory(Form("SUM::sigPdf_%s(fracGA1_%sc%d*pdfGA1_%s,fracGA2_%sc%d*pdfGA2_%s,pdfGA3_%s)", suffix.Data(), m_signalType.Data(), cateIndex, suffix.Data(), m_signalType.Data(), cateIndex, suffix.Data(), suffix.Data()));
  }
  // Bifurcated Gaussian shape:
  else if (function.Contains("BifurGA")) {
    m_ws->factory(Form("RooBifurGauss::sigPdf_%s(%s, prod::muGA_%s(muGANom_%s%s), prod::sigmaGALow_%s(sigmaGALowNom_%s%s), prod::sigmaGAHi_%s(sigmaGAHiNom_%s%s))", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
  }
  // Landau shape:
  else if (function.Contains("Landau")) {
    m_ws->factory(Form("RooLandau::sigPdf_%s(%s, prod::muLA_%s(muLANom_%s%s), prod::sigmaLA_%s(sigmaLANom_%s%s))", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
  }
  else if (function.Contains("CBPlusVoigt")) {
    // Cystal Ball component:
    m_ws->factory(Form("RooCBShape::pdfCB_%s(%s, prod::muCB_%s(muCBNom_%s%s), prod::sigmaCB_%s(sigmaCBNom_%s%s), alphaCB_%s, nCB_%sc%d)", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data(), suffix.Data(), m_signalType.Data(), cateIndex));
    // Voigtian component:
    m_ws->factory(Form("RooVoigtian::pdfVoigt_%s(%s, prod::muVoigt_%s(muVoigtNom_%s%s), prod::widthVoigt_%s(widthVoigtNom_%s%s), prod::sigmaVoigt_%s(sigmaVoigtNom_%s%s))", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
    
    // Crystal Ball + Voigtian:
    m_ws->factory(Form("SUM::sigPdf_%s(fracCB_%sc%d*pdfCB_%s,pdfVoigt_%s)", suffix.Data(), m_signalType.Data(), cateIndex, suffix.Data(), suffix.Data()));
  }
  else if (function.Contains("Voigt")) {
    m_ws->factory(Form("RooVoigtian::sigPdf_%s(%s, prod::muVoigt_%s(muVoigtNom_%s%s), prod::widthVoigt_%s(widthVoigtNom_%s%s), prod::sigmaVoigt_%s(sigmaVoigtNom_%s%s))", suffix.Data(), obsName.Data(), suffix.Data(), suffix.Data(), m_listMSS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data(), suffix.Data(), suffix.Data(), m_listMRS.Data()));
  }
  
  // Define the yield:
  if (!function.Contains("Parameterized")) {
    double yieldValue
      = (*m_ws->data(Form("data_%s",currKey.Data()))).sumEntries();
    m_ws->factory(Form("sigYield_%s%s[%f]",m_signalType.Data(),currKey.Data(),yieldValue));
  }
  
  if (m_verbose) {
    std::cout << "SigParam: Resonance successfully added." << std::endl;
  }
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
  TString workspaceName = Form("%s/res_%sworkspace.root", m_directory.Data(),
			       m_signalType.Data());
  m_ws->importClassCode();// Import the PDF classes (HggTwoSidedCBPdf, etc.)
  m_ws->writeToFile(workspaceName);
  if (m_verbose) {
    std::cout << "SigParam: Saved workspace to " << workspaceName << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Save a list of parameter names and values.
*/
void SigParam::saveParameterList() {
  TString paramListName = Form("%s/resonance_paramList.txt",m_directory.Data());
  std::ofstream outFile(paramListName);
  RooArgSet args = m_ws->allVars();
  TIterator *iterArgs = args.createIterator();
  RooRealVar *currIter = NULL;
  while ((currIter = (RooRealVar*)iterArgs->Next())) {
    outFile << currIter->GetName() << " " << currIter->getVal() << std::endl;
  }
  outFile.close();
  if (m_verbose) {
    std::cout << "SigParam: Saved param list to " << paramListName << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Save a list of signal yields in all categories and at all mass points. 
*/
void SigParam::saveYieldList() {
  TString yieldListName = Form("%s/resonance_yieldList.txt",m_directory.Data());
  std::ofstream outFile(yieldListName);
  // Loop over datasets:
  for (int i_d = 0; i_d < (int)m_massCatePairs.size(); i_d++) {
    double currMass = (m_massCatePairs[i_d]).first;
    int currCate = (m_massCatePairs[i_d]).second;
    outFile << currMass << " " << currCate << " "
	    << getYieldInCategory(currMass, currCate) << std::endl;
  }
  outFile.close();
  if (m_verbose) {
    std::cout << "SigParam: saved yield list to " << yieldListName << std::endl;
  }
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
  if (m_verbose) {
    std::cout << "SigParam: I/O directory set to " << directory << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Use a logarithmic Y axis for the plots.
   @param useLogYAxis - True iff plots should have log axis.
*/
void SigParam::setLogYAxis(bool useLogYAxis) {
  m_useLogYAxis = useLogYAxis;
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
   -----------------------------------------------------------------------------
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
   @param setConstant - True iff the mass value should be constant in the model.
   @param resonanceMass - The constant resonance mass value.
*/
void SigParam::setResMassConstant(bool setConstant, double resonanceMass) {
  // Either set the resonance mass fixed or floating in the workspace:
  m_ws->var("mResonance")->setConstant(setConstant);
  m_ws->var("mResonance")->setVal(resonanceMass);
  if (m_verbose) {
    std::cout << "SigParam: setResMassConstant(" << setConstant << ", " 
	      << resonanceMass << ")" << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   @param setConstant - True iff the mass value should be constant in the model.
*/
void SigParam::setResMassConstant(bool setConstant) {
  setResMassConstant(setConstant, m_ws->var("mResonance")->getVal());
}

/**
   -----------------------------------------------------------------------------
   Set the signal type to avoid collisions when using many production modes.
   @param signalType - The type of signal. 
*/
void SigParam::setSignalType(TString signalType) {
  if (signalType == "") m_signalType = "";
  else m_signalType = Form("%s_", signalType.Data());
  if (m_verbose) {
    std::cout << "SigParam: Signal type set to " << signalType << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Define the functional form of the parameterization in the resonance mass of 
   the given fit variable.
   @param varName - The name of the fit variable that has been parameterized.
   @param function - The functional form for parameterizing the given variable.
                     The parameters should be listed as @0, @1, @2, while the 
		     observable mass should be written as "obs" and the truth 
		     mass of the resonance should be "mRes". Examples are shown
		     in the class initialization.
*/
void SigParam::setVarParameterization(TString varName, TString function) {
  m_varParameterization[varName] = function;
  if (m_verbose) {
    std::cout << "SigParam: Setting " << varName << " parameterization to " 
	      << function << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Tell the program whether or not to use the same value for the Crystal Ball
   and Gaussian means.
   @param sameCBGAMean - True iff the same mean value should be used. 
*/
void SigParam::useCommonCBGAMean(bool sameCBGAMean) {
  m_sameCBGAMean = sameCBGAMean;
  if (m_verbose) {
    std::cout << "SigParam: Use common mean for CB and GA = " << sameCBGAMean
	      << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Get the variables associated with the given resonance PDF.
   @param function - The functional form of the resonance.
   @return - A list of the variables for the PDF.
*/
std::vector<TString> SigParam::variablesForFunction(TString function) {
  // Check that the function is defined before retrieving variables:
  if (!functionIsDefined(function)) {
    std::cout << "SigParam: Undefined resonance PDF " << function << std::endl;
    exit(0);
  }
  std::vector<TString> result; result.clear();
  // Crystal Ball + Gaussian-specific parameters:
  if (function.Contains("CBGA")) {
    result.push_back("muCBNom");
    result.push_back("sigmaCBNom");
    result.push_back("alphaCB");
    result.push_back("sigmaGANom");
    result.push_back("nCB");
    result.push_back("fracCB");
    if (!m_sameCBGAMean) {
      result.push_back("muGANom");
    }
  }
  // Double Crystal Ball-specific parameters:
  else if (function.Contains("DoubleCB")) {
    result.push_back("muCBNom");
    result.push_back("sigmaCBNom");
    result.push_back("alphaCBLo");
    result.push_back("alphaCBHi");
    result.push_back("nCBLo");
    result.push_back("nCBHi");
  }
  // Triple Gaussian-specific parameters:
  else if (function.Contains("GAx3")) {
    result.push_back("muGA1Nom");
    result.push_back("muGA2Nom");
    result.push_back("muGA3Nom");
    result.push_back("sigmaGA1Nom");
    result.push_back("sigmaGA2Nom");
    result.push_back("sigmaGA3Nom");
    result.push_back("fracGA1");
    result.push_back("fracGA2");
  }
  // Bifurcated Gaussian-specific parameters:
  else if (function.Contains("BifurGA")) {
    result.push_back("muGANom");
    result.push_back("sigmaGALowNom");
    result.push_back("sigmaGAHiNom");
  }
  // Landau-specific parameters:
  else if (function.Contains("Landau")) {
    result.push_back("muLANom");
    result.push_back("sigmaLANom");
  }
  // Crystal-Ball plus Voigtian specific parameters:
  else if (function.Contains("CBPlusVoigt")) {
    result.push_back("muCBNom");
    result.push_back("sigmaCBNom");
    result.push_back("alphaCB");
    result.push_back("nCB");
    result.push_back("fracCB");
    result.push_back("muVoigtNom");
    result.push_back("widthVoigtNom");
    result.push_back("sigmaVoigtNom"); 
    
  }
  // Voigtian specific parameters:
  else if (function.Contains("Voigt")) {
    result.push_back("muVoigtNom");
    result.push_back("widthVoigtNom");
    result.push_back("sigmaVoigtNom"); 
  }
  else {
    std::cout << "SigParam: Undefined resonance PDF " << function << std::endl;
    exit(0);
  }

  return result;
}

/**
   -----------------------------------------------------------------------------
   Sets the std::cout print level. 
   @param beVerbose - True iff you want a lot of spammy text.
*/
void SigParam::verbosity(bool beVerbose) {
  m_verbose = beVerbose;
  if (m_verbose) {
    std::cout << "SigParam: m_verbosity = " << m_verbose << std::endl;
  }
}
