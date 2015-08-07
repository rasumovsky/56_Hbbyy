////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHWorkspace.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 06/07/2015                                                          //
//                                                                            //
//  This class builds the workspace for the resonant and non-resonant di-     //
//  Higgs search at 13 TeV.                                                   //
//                                                                            //
//  Signal types                                                              //
//    - "NonRes"                                                              //
//    - "ResMx300", "ResMx350", ...                                           //
//                                                                            //
//  Analysis types:                                                           //
//    - "NonRes", "Res"
//                                                                            //
//  options:                                                                  //
//    - "ResonantOnly": only generate the resonant search workspace.          //
//    - "NonResOnly": only generate the nonresonant search workspace.         //
//                                                                            //
//  Job options: "New", "FromFile" determine whether to create a new workspace//
//  or load a previously generated one. The systematics are controlled by the //
//  "nonorm", "nopes", "noper", "noss", "nobgm", "nomig", and "nosys" options.//
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHWorkspace.h"

using namespace RooFit;
using namespace RooStats;
using namespace CommonFunc;

/**
   -----------------------------------------------------------------------------
   Instantiate the class.
   @param newConfigFile - The name of the analysis config file.
   @param newAnalysisType - The name of the analysis type.
   @param newOptions - The job options ("New", "FromFile"), etc.
   @returns void
*/
DHWorkspace::DHWorkspace(TString newConfigFile, TString newAnalysisType,
			 TString newOptions) {
  m_configFile = newConfigFile;
  m_anaType = newAnalysisType;
  m_options = newOptions;
  m_allGoodFits = true;
  
  m_combinedWS = NULL;
  m_modelConfig = NULL;
  
  m_muNominalSH = 1;
  
  m_config = new Config(m_configFile);
  m_dataToPlot = (m_config->getBool("doBlind")) ? "asimovDataMu1" : "obsData";
  
  // Print workspace inputs:
  std::cout << "\nDHWorkspace: Initializing..."
	    << "\n\tconfigFile = " << m_configFile 
	    << "\n\tanaType = " << m_anaType
	    << "\n\toptions = " << m_options << std::endl;
  
  // Assign output directory, and make sure it exists:
  m_outputDir = Form("%s/%s/DHWorkspace", 
		     (m_config->getStr("masterOutput")).Data(),
		     (m_config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  system(Form("mkdir -vp %s/Plots/", m_outputDir.Data()));
  system(Form("mkdir -vp %s/rootfiles/", m_outputDir.Data()));
  system(Form("mkdir -vp %s/mu/", m_outputDir.Data()));
  
  // Set style for plots:
  CommonFunc::SetAtlasStyle();
    
  // Make new or load old workspace:
  if (m_options.Contains("FromFile")) loadWSFromFile();
  else createNewWS();
  
  if (fitsAllConverged()) {
    std::cout << "DHWorkspace: Successfully initialized!" << std::endl;
  }
  else {
    std::cout << "DHWorkspace: Fit failure during initialization." << std::endl;
  }
  return;
}

/**
   -----------------------------------------------------------------------------
   Checks whether all of the fits converged.
   @returns - true iff the fits all converged.
*/
bool DHWorkspace::fitsAllConverged() {
  return m_allGoodFits;
}

/**
   -----------------------------------------------------------------------------
   Retrieves the workspace created by this program.
*/
RooWorkspace* DHWorkspace::getCombinedWorkspace() {
  return m_combinedWS;
}

/**
   -----------------------------------------------------------------------------
   Retrieves a pointer to the model config.
*/
ModelConfig* DHWorkspace::getModelConfig() {
  return m_modelConfig;
}

/**
   -----------------------------------------------------------------------------
   Load a previously created workspace.
*/
void DHWorkspace::loadWSFromFile() {
  //Check to see if the workspace has actually been made.
  TFile inputFile(Form("%s/rootfiles/workspaceDH_%s.root", 
		       m_outputDir.Data(), m_anaType.Data()), "read");
  if (inputFile.IsOpen()) {
    std::cout << "DHWorkspace: Loading workspace from file..." << std::endl;
    // Load the single workspace file:
    m_combinedWS = (RooWorkspace*)inputFile.Get("combinedWS");
    m_modelConfig = (ModelConfig*)m_combinedWS->obj("modelConfig");
  }
  else {
    std::cout << "DHWorkspace: WARNING! Cannot locate requested workspace!" 
	      << std::endl;
    createNewWS();
  }
}

/**
   -----------------------------------------------------------------------------
   Create a workspace from scratch. 
*/
void DHWorkspace::createNewWS() {
  
  std::cout << "DHWorkspace: Create a new workspace from scratch." << std::endl;
  std::cout << "\n........................................" << std::endl;
  std::cout << "Luminosity at 13 TeV: "
	    << m_config->getStr("analysisLuminosity") << " pb-1." << std::endl;
  
  // Define and name analysis categories:
  std::vector<TString> cateNames
    = m_config->getStrV(Form("cateNames%s", m_anaType.Data()));
  int nCategories = cateNames.size();
  std::cout << "  Number of categories = " << nCategories << std::endl;
  std::cout << "........................................" << std::endl;
  
  // Read tables of PES and PER and store values:
  m_pes = new PESReader(m_config->getStr("fileNamePESValues"), nCategories);
  m_per = new PERReader(m_config->getStr("fileNamePERValues"), nCategories);
  
  //--------------------------------------//
  // Initialize classes relevant to workspace:
  // The combined workspace:
  m_combinedWS = new RooWorkspace("combinedWS");
  m_combinedWS->importClassCode();
  
  // Category workspaces, and RooCategory and Simultaneous PDF:
  RooWorkspace* cateWS[nCategories];
  RooCategory *categories = new RooCategory("categories", "categories");
  RooSimultaneous *combinedPdf
    = new RooSimultaneous("combinedPdf", "combinedPdf", *categories);
  
  // Instantiate parameter sets:
  RooArgSet *nuisanceParameters = new RooArgSet();
  RooArgSet *muSHConstants = new RooArgSet();
  RooArgSet *globalObservables = new RooArgSet();
  RooArgSet *observables = new RooArgSet();
  RooArgSet *constraints = new RooArgSet();
  
  // Maps for datasets:
  std::map<std::string, RooDataSet*> dm;
  map<string,RooDataSet*> dmAsimovMu0;
  map<string,RooDataSet*> dmAsimovMu1;
  
  //--------------------------------------//
  // Loop over channels:
  std::cout << "DHWorkspace: Looping over categories to define workspace."
	    << std::endl;
  for (m_currCateIndex = 0; m_currCateIndex < nCategories; m_currCateIndex++) {
    
    m_currCateName = cateNames[m_currCateIndex];
    
    // Create the workspace for a single category:
    cateWS[m_currCateIndex] = createNewCategoryWS();
    categories->defineType(m_currCateName);
    
    // Add PDFs and parameters: 
    TString namePdf = Form("model_%s", m_currCateName.Data());
    TString nameNP = Form("nuisanceParameters_%s", m_currCateName.Data());
    TString nameGlob = Form("globalObservables_%s", m_currCateName.Data());
    TString nameMuC = Form("muConstants_%s", m_currCateName.Data());
    TString nameObs = Form("observables_%s", m_currCateName.Data());
    combinedPdf->addPdf(*cateWS[m_currCateIndex]->pdf(namePdf), m_currCateName);
    nuisanceParameters->add(*cateWS[m_currCateIndex]->set(nameNP));
    globalObservables->add(*cateWS[m_currCateIndex]->set(nameGlob));
    muSHConstants->add(*cateWS[m_currCateIndex]->set(nameMuC));
    nuisanceParameters->add(*cateWS[m_currCateIndex]->set(nameMuC));
    observables->add(*cateWS[m_currCateIndex]->set(nameObs));

    // Add category datasets to combined workspace and combined datasets:
    TString nameOD = Form("obsData_%s", m_currCateName.Data());
    //TString nameAD0 = Form("asimovDataMu0_%s", m_currCateName.Data());
    //TString nameAD1 = Form("asimovDataMu1_%s", m_currCateName.Data());
    m_combinedWS->import(*(RooDataSet*)cateWS[m_currCateIndex]->data(nameOD));
    // m_combinedWS->import(*(RooDataSet*)cateWS[m_currCateIndex]->data(nameAD0));
    //m_combinedWS->import(*(RooDataSet*)cateWS[m_currCateIndex]->data(nameAD1));
    dm[(string)m_currCateName] = (RooDataSet*)m_combinedWS->data(nameOD);
    //dmAsimovMu0[(string)m_currCateName] 
    //  = (RooDataSet*)m_combinedWS->data(nameAD0);
    //dmAsimovMu1[(string)m_currCateName]
    //  = (RooDataSet*)m_combinedWS->data(nameAD1);
  }
  std::cout << "DHWorkspace: Beginning to combine all categories." << std::endl;
  
  // Define the combined dataset:
  RooRealVar wt("wt", "wt", 1);
  RooArgSet *args = new RooArgSet();
  args->add(*observables);
  args->add(wt);
  RooDataSet* obsData = new RooDataSet("obsData", "obsData", *args,
				       Index(*categories), Import(dm),
				       WeightVar(wt));
  //  RooDataSet* asimovDataMu0 = new RooDataSet("asimovDataMu0", "asimovDataMu0",
  //					     *args, Index(*categories), 
  //					     Import(dmAsimovMu0),WeightVar(wt));
  //RooDataSet* asimovDataMu1 = new RooDataSet("asimovDataMu1", "asimovDataMu1",
  //					     *args, Index(*categories), 
  //					     Import(dmAsimovMu1),WeightVar(wt));
  
  // Import PDFs, parameters, and dataset into workspace:
  m_combinedWS->import(*categories);
  m_combinedWS->import(*combinedPdf);
  m_combinedWS->defineSet("nuisanceParameters", *nuisanceParameters);
  m_combinedWS->defineSet("observables", *observables);
  m_combinedWS->defineSet("globalObservables", *globalObservables);
  m_combinedWS->defineSet("poi", RooArgSet(*m_combinedWS->var("mu_DH")));   
  m_combinedWS->defineSet("muSHConstants", *muSHConstants);
  m_combinedWS->import(*obsData);
  //m_combinedWS->import(*asimovDataMu0);
  //m_combinedWS->import(*asimovDataMu1);
  
  // Define the ModelConfig for the analysis and import to the workspace:
  m_modelConfig = new ModelConfig("modelConfig", m_combinedWS);
  m_modelConfig->SetPdf((*m_combinedWS->pdf("combinedPdf")));
  m_modelConfig->SetObservables((*m_combinedWS->set("observables")));
  m_modelConfig->SetParametersOfInterest((*m_combinedWS->set("poi")));
  m_modelConfig
    ->SetNuisanceParameters((*m_combinedWS->set("nuisanceParameters")));
  m_modelConfig
    ->SetGlobalObservables((*m_combinedWS->set("globalObservables")));
  m_combinedWS->import(*m_modelConfig);
  
  std::cout << "DHWorkspace: Printing the combined workspace." << std::endl;
  m_combinedWS->Print("v");
  
  
  //////////////
  (*m_combinedWS->var("mu_DH")).setVal(0.0);
  
  (*m_combinedWS->pdf("combinedPdf")).fitTo(*m_combinedWS->data("obsData"), Minos(RooArgSet(*m_combinedWS->set("nuisanceParameters"))), SumW2Error(kTRUE));
  
  //Minos(RooArgSet(*nuisCateWS))
  
  
  // Stupid secondary loop for Asimov data:
  for (m_currCateIndex = 0; m_currCateIndex < nCategories; m_currCateIndex++) {
    m_currCateName = cateNames[m_currCateIndex];
    
    (*m_combinedWS->var("nBkg_"+m_currCateName)).setVal((*m_combinedWS->data(Form("obsData_%s",m_currCateName.Data()))).sumEntries());
    dmAsimovMu0[(string)m_currCateName] = createAsimovData(0, m_muNominalSH);
    dmAsimovMu1[(string)m_currCateName] = createAsimovData(1, m_muNominalSH);
    
    //TString nameAD0 = Form("asimovDataMu0_%s", m_currCateName.Data());
    //TString nameAD1 = Form("asimovDataMu1_%s", m_currCateName.Data());
    //m_combinedWS->import(*(RooDataSet*)cateWS[m_currCateIndex]->data(nameAD0));
    //m_combinedWS->import(*(RooDataSet*)cateWS[m_currCateIndex]->data(nameAD1));
    //dmAsimovMu0[(string)m_currCateName] 
    //= (RooDataSet*)m_combinedWS->data(nameAD0);
    //dmAsimovMu1[(string)m_currCateName]
    //= (RooDataSet*)m_combinedWS->data(nameAD1);
  }    
  RooDataSet* asimovDataMu0 = new RooDataSet("asimovDataMu0", "asimovDataMu0",
  					     *args, Index(*categories), 
  					     Import(dmAsimovMu0),WeightVar(wt));
  RooDataSet* asimovDataMu1 = new RooDataSet("asimovDataMu1", "asimovDataMu1",
  					     *args, Index(*categories), 
  					     Import(dmAsimovMu1),WeightVar(wt));
  
  m_combinedWS->import(*asimovDataMu0);
  m_combinedWS->import(*asimovDataMu1);
  //////////
  
  // Save snapshot of original parameter values:
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*m_modelConfig->GetNuisanceParameters());
  poiAndNuis->add(*m_modelConfig->GetParametersOfInterest());
  m_combinedWS->saveSnapshot("paramsOrigin", *poiAndNuis);
  
  // Start profiling the data:
  std::cout << "DHWorkspace: Start profiling data" << std::endl;
    
  // Choose two models for the default test fits:
  std::vector<TString> sigDHModes = m_config->getStrV("sigDHModes");
  TString currDHSignal
    = m_anaType.Contains("NonRes") ? sigDHModes[0] : sigDHModes[1];
  
  DHTestStat *dhts 
    = new DHTestStat(m_configFile, currDHSignal, "FromFile", m_combinedWS);
  dhts->saveSnapshots(true);
  dhts->setPlotDirectory(Form("%s/Plots/", m_outputDir.Data()));
  double profiledMuDHVal = -999.0;
  // Mu = 0 fits:
  double nllMu0 = dhts->getFitNLL(m_dataToPlot, 0, true, profiledMuDHVal);
  if (!dhts->fitsAllConverged()) m_allGoodFits = false;
  // Mu = 1 fits:
  double nllMu1 = dhts->getFitNLL(m_dataToPlot, 1, true, profiledMuDHVal);
  if (!dhts->fitsAllConverged()) m_allGoodFits = false;
  // Mu free fits:
  double nllMuFree = dhts->getFitNLL(m_dataToPlot, 1, false, profiledMuDHVal);
  if (!dhts->fitsAllConverged()) m_allGoodFits = false;
  
  // Print summary of the fits:
  std::cout.precision(10);
  std::cout << "\nDHWorkspace: Printing likelihood results for " << m_anaType
	    << std::endl;
  std::cout << "\tnll(muDH = 1):  " << nllMu1 << std::endl;
  std::cout << "\tnll(muDH = 0):  " << nllMu0 << std::endl;
  std::cout << "\tnll(muDH free): " << nllMuFree << std::endl;
  std::cout << " " << endl;
  std::cout << "\tnll(S+B)/nll(B) " << nllMu1 - nllMu0 << std::endl;
  std::cout << "\tnll(muDH=1)/nll(muhat) = " << nllMu1-nllMuFree << std::endl;
  std::cout << "\tnll(muDH=0)/nll(muhat) = " << nllMu0-nllMuFree << std::endl;
  if (m_allGoodFits) std::cout << "all good fits = TRUE" << std::endl;
  else std::cout << "all good fits = FALSE" << std::endl;
  std::cout << "Profiled muDH value : " << profiledMuDHVal << std::endl;
  
  // Write the profiled mu value to file:
  ofstream fileMuProf;
  fileMuProf
    .open(Form("%s/mu/mu_%s.txt", m_outputDir.Data(), currDHSignal.Data()));
  fileMuProf << profiledMuDHVal << std::endl;
  fileMuProf.close();

  // Print the workspace before saving:
  std::cout << "DHWorkspace: Printing the workspace to be saved." << std::endl;
  m_combinedWS->Print("v");
  
  // Write workspace to file:
  m_combinedWS->writeToFile(Form("%s/rootfiles/workspaceDH_%s.root",
				 m_outputDir.Data(), m_anaType.Data()));
}

/**
   -----------------------------------------------------------------------------
   Create the workspace for a single analysis category. Note: m_anaType, 
   m_currCateIndex, and m_currCateName are defined with global scope and updated
   in each call to this method.
*/
RooWorkspace* DHWorkspace::createNewCategoryWS() {
    
  // The bools that control the systematic uncertainties:
  bool inclusive = (m_currCateName == "inclusive");
  bool channel_constraints_attached = (m_currCateIndex == 0);
  bool switch_norm = !m_options.Contains("nonorm");
  bool switch_pes = !m_options.Contains("nopes");
  bool switch_per = !m_options.Contains("noper");
  bool switch_ss  = !m_options.Contains("noss");
  bool switch_bgm = !m_options.Contains("nobgm");
  bool switch_mig = !m_options.Contains("nomig");
  bool switch_nosys = m_options.Contains("nosys");
  if (switch_nosys) {
    std::cout << "\tDHWorkspace: ALL systematics OFF" << endl;
    switch_norm = false;   switch_pes = false;   switch_per = false;
    switch_ss = false;     switch_bgm = false;   switch_mig = false;
  }
  std::cout << "\tNormalization systematics = " << switch_norm << std::endl;
  std::cout << "\tEnergy scale systematics  = " << switch_pes  << std::endl;
  std::cout << "\tResolution systematics    = " << switch_per  << std::endl;
  std::cout << "\tShape systematics         = " << switch_ss   << std::endl;
  std::cout << "\tBackground systematics    = " << switch_bgm  << std::endl;
  std::cout << "\tMigration systematics     = " << switch_mig  << std::endl;
  
  //--------------------------------------//
  // Create the individual channel workspace:
  RooWorkspace *tempWS
    = new RooWorkspace(Form("tempWS_%s", m_currCateName.Data()));
  
  // Nuisance parameters:
  RooArgSet *nuisParamsBkg = new RooArgSet();
  RooArgSet *nuisParamsCorr = new RooArgSet();
  RooArgSet *nuisParamsUncorr = new RooArgSet();
  // Constraint parameters:
  RooArgSet *constraints = new RooArgSet();
  RooArgSet *constraintsBias = new RooArgSet();
  // Global observables:
  RooArgSet *globalObs = new RooArgSet();
  RooArgSet *globalObsProc = new RooArgSet();
  // Expected terms:
  RooArgSet *expectedShape = new RooArgSet();
  RooArgSet *expectedBias = new RooArgSet();
  RooArgSet *expected = new RooArgSet();
  RooArgSet *expectedSH = new RooArgSet();
  RooArgSet *expectedDH = new RooArgSet();
  
  // array setup[5] is used to configure a nuisance parameter
  // [0]    [1]       [2]   [3]     
  // sigma, sigmalow, beta, nominal,
  
  //--------------------------------------//
  // Normalization systematics:
  if (switch_norm) {
    double setupLumi[4] = {0.036, 0, 1, 1};
    makeNP("Luminosity", setupLumi, *&nuisParamsCorr, *&constraints, 
	   *&globalObs, *&expected);
    double setupTrigger[4] = {0.005, 0, 1, 1};
    makeNP("Trigger", setupTrigger, *&nuisParamsCorr, *&constraints, 
	   *&globalObs, *&expected);
    double setupIsEM[4] = {0.0526, 0, 1, 1};
    makeNP("PhotonID", setupIsEM, *&nuisParamsCorr, *&constraints, *&globalObs,
	   *&expected);
    double setupIso[4] = {0.004, 0, 1, 1};
    makeNP("Isolation", setupIso, *&nuisParamsCorr, *&constraints, *&globalObs,
	   *&expected);
    double setupESCALE[4] = {0.003, 0, 1, 1};
    makeNP("ESCALE", setupESCALE, *&nuisParamsCorr, *&constraints, *&globalObs,
	   *&expected);
  }
  
  //--------------------------------------//
  // Migration systematics:
  if (switch_mig) {
    // Follow examples in other workspaces...
  }
  
  //--------------------------------------//
  // SYSTEMATICS: Spurious signal
  if (switch_bgm) {
    double ssEvents = 0.1;//need a spuriousSignal function;
    double setupBias[4] = {ssEvents, -999, 1, 0}; //Gaussian constraint
    makeNP("bias", setupBias, *&nuisParamsUncorr, *&constraintsBias,
	   *&globalObs, *&expectedBias);
    RooProduct sigBias("sigBias","sigBias",*expectedBias);
    tempWS->import(sigBias);
  }
  else tempWS->factory("sigBias[0]");//expectedBias
  
  //--------------------------------------//
  // SYSTEMATICS - Resolution:
  TString m_listMRS = "";
  std::vector<TString> perList; perList.clear();
  if (switch_per) {
    double setupPER[4] = {0.0, 0, 1, 1};
    // Loop over sources of resolution systematic uncertainty:
    for (int i_s = 0; i_s < m_per->getNumberOfSources(); i_s++) {
      TString currPERSource = m_per->getNameOfSource(i_s);
      TString currPERName = Form("EM_%s",currPERSource.Data());
      perList.push_back(currPERName);
      setupPER[0] = m_per->getValue(currPERSource, m_currCateIndex);
      setupPER[2] = m_per->getSign(currPERSource, m_currCateIndex);
      // Resolution on the inclusive shape:
      makeShapeNP(currPERName, "DH", setupPER, *&nuisParamsCorr, *&constraints,
		  *&globalObs, *&expectedShape);
      makeShapeNP(currPERName, "SH", setupPER, *&nuisParamsCorr, *&constraints,
		  *&globalObs, *&expectedShape);
    }
  }
  
  //--------------------------------------//
  // SYSTEMATICS - Energy-scale:
  TString m_listMSS = "";
  std::vector<TString> pesList; pesList.clear();
  if (switch_pes) {
    double setupPES[4] = {0.0, 0, 1, 1};
    // loop over sources of energy scale systematic uncertainty:
    for (int i_s = 0; i_s < m_pes->getNumberOfSources(); i_s++) {
      TString currPESSource = m_pes->getNameOfSource(i_s);
      TString currPESName = Form("EM_%s",currPESSource.Data());
      pesList.push_back(currPESName);
      setupPES[0] = m_pes->getValue(currPESSource, m_currCateIndex);
      setupPES[2] = m_pes->getSign(currPESSource, m_currCateIndex);
      makeNP(currPESName, setupPES, *&nuisParamsCorr, *&constraints, 
	     *&globalObs, *&expectedShape);
    }
  }
  
  //--------------------------------------//
  // Parameters of interest (POIs):
  double muMin = 0.0; double muMax = 1000.0;
  RooRealVar *mu_DH = new RooRealVar("mu_DH", "mu_DH", 1, muMin, muMax);
  RooRealVar *mu_SH = new RooRealVar("mu_SH", "mu_SH", 1, muMin, muMax);
  expectedDH->add(RooArgSet(*mu_DH));
  expectedSH->add(RooArgSet(*mu_SH));
  
  // Expectation values:
  RooProduct expectationDH("expectationDH", "expectationDH", *expectedDH);
  RooProduct expectationSH("expectationSH", "expectationSH", *expectedSH);
  RooProduct expectationCommon("expectationCommon","expectationCommon",
			       *expected);
  
  // Spurious signal term will assume the shape of "inclusive" pdf.
  tempWS->import(expectationCommon);
  tempWS->import(expectationSH);
  tempWS->import(expectationDH);
  tempWS->import(*expectedShape);
  tempWS->import(*expectedBias);
  
  // Declare the observable, and the observables set:
  TString obsVarName = m_anaType.EqualTo("NonRes") ? "m_yy" : "m_bbyy";
  double obsMin = m_anaType.EqualTo("NonRes") ?
    m_config->getNum("DHMyyRangeLo") : m_config->getNum("DHMyybbRangeLo");
  double obsMax = m_anaType.EqualTo("NonRes") ? 
    m_config->getNum("DHMyyRangeHi") : m_config->getNum("DHMyybbRangeHi");
  tempWS->factory(Form("%s[%f,%f]", obsVarName.Data(), obsMin, obsMax));
  tempWS->defineSet("obsprelim", obsVarName);
  
  //--------------------------------------//
  // Begin PDF construction.
    
  // Construct the background PDF (present in every category):
  std::cout << "DHWorkspace: Constructing background model" << std::endl;
  BkgModel *currBkgModel = new BkgModel(tempWS->var(obsVarName));
  vector<TString> bkgFunctionList
    = m_config->getStrV(Form("bkgFunctions%s",m_anaType.Data()));
  TString bkgFuncName = bkgFunctionList[m_currCateIndex];
  currBkgModel->addBkgToCateWS(tempWS, nuisParamsBkg, bkgFuncName);
    
  std::cout << "DHWorkspace: Printing corr. bkg. nuisance params." << std::endl;
  nuisParamsBkg->Print("v");
  
  // Start to construct the statistical model, add components one by one:
  TString modelForm = "SUM::modelSB(nBkg[100,0,1000000]*bkgPdf";
  
  // Construct SH and DH signals (depending on categories):
  std::cout << "DHWorkspace: Constructing signal model" << std::endl;
  if (m_anaType.EqualTo("NonRes")) {
    
    // SH (Single-Higgs) signal:
    if (DHAnalysis::cateHasComponent(m_config,m_currCateName,"BkgSingleHiggs")){
      std::cout << "DHWorkspace: Add single Higgs model component for category "
		<< m_currCateName << std::endl;
      // jj CR first, bb SR second
      double normSH[2] = {0.00868, 0.186};
      tempWS->factory(Form("RooCBShape::pdfCB_SH(%s, prod::muCB_SH(muCBNom_SH[124.96]%s), prod::sigmaCB_SH(sigmaCBNom_SH[1.53]%s), alphaCB_SH[1.56], nCB_SH[10.0])", obsVarName.Data(), m_listMSS.Data(), m_listMRS.Data()));
      tempWS->factory(Form("RooGaussian::pdfGA_SH(%s, prod::muGA_SH(muCBNom_SH%s), prod::sigmaGA_SH(sigmaGANom_SH[23.98]%s))", obsVarName.Data(), m_listMSS.Data(), m_listMRS.Data()));
      tempWS->factory(Form("SUM::sigPdfSH(fracCB_SH[1.00]*pdfCB_SH,pdfGA_SH)"));
      tempWS->factory(Form("nSH[%f]", normSH[m_currCateIndex]));
      
      // Normalization for SH (expectationCommon = mu*isEM*lumi*migr):
      tempWS->factory("prod::nSigSH(nSH,expectationCommon,expectationSH)");
      
      // Add to the statistical model:
      modelForm += ",nSigSH*sigPdfSH";
    }
    
    // Now BSM (Di-Higgs) signal:
    if (DHAnalysis::cateHasComponent(m_config,m_currCateName,"SigDiHiggs")){
      std::cout << "DHWorkspace: Add di-Higgs model component for category "
		<< m_currCateName << std::endl;
      tempWS->factory(Form("RooCBShape::pdfCB_DH(%s, prod::muCB_DH(muCBNom_DH[124.95]%s), prod::sigmaCB_DH(sigmaCBNom_DH[1.31]%s), alphaCB_DH[1.56], nCB_DH[10.0])", obsVarName.Data(), m_listMSS.Data(), m_listMRS.Data()));
      tempWS->factory(Form("RooGaussian::pdfGA_DH(%s, prod::muGA_DH(muCBNom_DH%s), prod::sigmaGA_DH(sigmaGANom_DH[2.84]%s))", obsVarName.Data(), m_listMSS.Data(), m_listMRS.Data()));
      tempWS->factory(Form("SUM::sigPdfDH(fracCB_DH[0.94]*pdfCB_DH,pdfGA_DH)"));
      //tempWS->factory("nDH[1,0,100]");
      tempWS->factory("nDH[1.0]");
      
      // Normalization for DH (expectationCommon = mu*isEM*lumi*migr):
      tempWS->factory("prod::nSigDH(nDH,expectationCommon,expectationDH)");
      
      // Add to the statistical model:
      modelForm += ",nSigDH*sigPdfDH,sigBias*sigPdfDH";
    }
  }
  else if (m_anaType.EqualTo("Resonant")) {
    std::cout << "DHWorkspace: Resonant analysis has not been implemented!"
	      << std::endl;
    exit(0);
  }
  
  // Finish constructing the model form and then add to workspace:
  modelForm += ")";
  tempWS->factory(modelForm);
  
  std::cout << "DHWorkspace: Model formed, now handling params." << std::endl;
  
  //--------------------------------------//
  // Attach constraint terms:
  
  // Only attach constraint term to first category. If constraint terms were
  // attached to each category, constraints would effectively be multiplied.
  if (m_currCateIndex == 0) {
    constraints->add(*constraintsBias);
    RooProdPdf constraint("constraint", "constraint", *constraints);
    tempWS->import(constraint);
    tempWS->factory("PROD::model(modelSB,constraint)");
  }
  
  // Except in the case where the constraints are uncorrelated between
  // categories, as with the spurious signal:
  else {
    RooProdPdf constraint("constraint", "constraint", *constraintsBias);
    tempWS->import(constraint);
    tempWS->factory("PROD::model(modelSB,constraint)");
  }
  
  //--------------------------------------//
  // Handle nuisance parameters, constraint term correlations and naming:
  
  /*
    Specify the group of nuisance parameters that are correlated between
    categories. Technically, this is done by sharing the same name for nuisance
    parameter between sub-channels. Their respective global observables should
    also share the same name. nuisParamsCorr should contain all correlated 
    nuisance parameters. All uncorrelated nuisance parameters should be included
    in nuisParamsUncorr.
  */
  TString corrNPNames = "mu_DH,mu_SH";
  
  // Iterate over nuisance parameters:
  TIterator *iterNuis = nuisParamsCorr->createIterator();
  RooRealVar* currNuis;
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    std::cout << "\t" << currNuis->GetName() << std::endl;
    corrNPNames += Form(",nuisPar_%s,globOb_%s",currNuis->GetName(),
			currNuis->GetName());
  }
  
  // Iterate over correlated nuisance parameters:
  TIterator *iterBkgNuis = nuisParamsBkg->createIterator();
  RooRealVar* currBkgNuis;
  while ((currBkgNuis = (RooRealVar*)iterBkgNuis->Next())) {
    std::cout << "\t" << currBkgNuis->GetName() << std::endl;
    corrNPNames += Form(",%s",currBkgNuis->GetName());
  }
  std::cout << "For category " << m_currCateName << ", correlate variables: "
	    << corrNPNames << std::endl;
  
  /*
    Sub-channel labeling
    Import the workspace tempWS to a new workspace and add m_currCateName as a 
    suffix to all nodes and variables. The correlated nuisance parameters
    and their respective global observables will not be renamed.
  */
  RooWorkspace* categoryWS = new RooWorkspace("workspace_"+m_currCateName);
  categoryWS->import((*tempWS->pdf("model")), RenameAllNodes(m_currCateName),
		     RenameAllVariablesExcept(m_currCateName,corrNPNames),
		     Silence());
  
  // Adding correlated nuisance parameters to nuisanceParameters:
  RooArgSet* nuisCateWS = new RooArgSet();
  iterNuis->Reset();
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    nuisCateWS->add(*(RooRealVar*)categoryWS->obj(currNuis->GetName()));
  }
  
  // Adding uncorrelated nuisance parameters to nuisanceParameters:
  TIterator *iterNuisUncorrelated = nuisParamsUncorr->createIterator();
  RooRealVar* currNuisUncorrelated;
  while ((currNuisUncorrelated = (RooRealVar*)iterNuisUncorrelated->Next())) {
    TString nuisName = (currNuisUncorrelated->GetName() 
			+ (TString)"_" + m_currCateName);
    nuisCateWS->add(*(RooRealVar*)categoryWS->obj(nuisName));
  }
  


  // Adding unconstrained NPs from the background pdf:
  RooArgSet* nuisBkgCateWS = new RooArgSet();
  TIterator *iterNuisBkg = nuisParamsBkg->createIterator();
  RooRealVar* currNuisBkg;
  while ((currNuisBkg = (RooRealVar*)iterNuisBkg->Next())) {
    TString parName = currNuisBkg->GetName()+(TString)"_"+m_currCateName;
    nuisCateWS->add(*categoryWS->var(parName));
    //nuisBkgCateWS->add(*categoryWS->var(parName));
    //nuisBkgCateWS->add(*categoryWS->var(currNuisBkg->GetName()));
  }
  
  /*
    Global observables:
    Global observables only appear in the constraint terms. All constraint terms
    of correlated nuisance parameters are attached to the pdf of the first
    subchannel. For those global observables, their names should be the same as
    those in the w. For other subchannels, only the bias constraint term is
    attached.
  */  
  RooArgSet *globsCateWS = new RooArgSet();
  TIterator *iterGlobs = globalObs->createIterator();
  RooRealVar *currGlobs;
  while ((currGlobs = (RooRealVar*)iterGlobs->Next())) {
    TString globName = currGlobs->GetName()+(TString)"_"+m_currCateName;
    if (categoryWS->obj(globName)) {
      globsCateWS->add(*(RooRealVar*)categoryWS->obj(globName));
      categoryWS->var(globName)->setConstant();
    }
    else if (categoryWS->obj(currGlobs->GetName())) {
      globsCateWS->add(*(RooRealVar*)categoryWS->obj(currGlobs->GetName()));
      categoryWS->var(currGlobs->GetName())->setConstant();
    }
  }
  
  /*
    Observables:
    Iterate over the observables in this category and add them to the new set.
  */
  RooArgSet *obsCateWS = new RooArgSet();
  TIterator *iterObs = tempWS->set("obsprelim")->createIterator();
  RooRealVar *currObs;
  while ((currObs = (RooRealVar*)iterObs->Next())) {
    TString currObsName = currObs->GetName()+(TString)"_"+m_currCateName;
    if (categoryWS->obj(currObsName)) {
      obsCateWS->add(*(RooRealVar*)categoryWS->obj(currObsName));
    }
    else {
      obsCateWS->add(*(RooRealVar*)categoryWS->obj(currObs->GetName()));
    }
  }
  
  // Set the single Higgs mu parameters to 1 and constant:
  std::cout << "DHWorkspace: Setting SH signals constant." << std::endl;
  RooArgSet* muConstCateWS = new RooArgSet();
  muConstCateWS->add(*categoryWS->var("mu_SH"));
  TIterator *iterMuConst = muConstCateWS->createIterator();
  RooRealVar *currMuConst;
  while ((currMuConst = (RooRealVar*)iterMuConst->Next())) {
    currMuConst->setVal(1.0);
    currMuConst->setConstant(true);
  }
  
  categoryWS->defineSet(Form("muConstants_%s",m_currCateName.Data()),
			*muConstCateWS);
  categoryWS->defineSet(Form("observables_%s",m_currCateName.Data()),
			*obsCateWS);
  categoryWS->defineSet(Form("nuisanceParameters_%s",m_currCateName.Data()),
			*nuisCateWS);
  categoryWS->defineSet(Form("globalObservables_%s",m_currCateName.Data()),
			*globsCateWS);
  
  //--------------------------------------//
  // Import the observed data set, and create a binned version:
  DHDataReader *dhdr = new DHDataReader(m_configFile, categoryWS->var(Form("%s_%s",obsVarName.Data(),m_currCateName.Data())));
  RooDataSet *obsData = NULL;
  if (m_anaType.EqualTo("NonRes")) {
    obsData = dhdr->loadNonResData(m_currCateName);
  }
  else {
    std::cout << "DHWorkspace: Error importing resonant data." << std::endl;
  }
  // Rename the imported dataset for consistency:
  TString obsDataName = Form("obsData_%s", m_currCateName.Data());
  obsData->SetNameTitle(obsDataName, obsDataName);
  categoryWS->import(*obsData);
  
  /////////////
  // Need to make sure shape of bkg always comes from CR
  /////////////

  // Fit background shape and set normalization:
  std::cout << "DHWorkspace: Fitting bkg. in " << m_currCateName << std::endl;
  (*categoryWS->var("nBkg_"+m_currCateName)).setVal(obsData->sumEntries());
  
  // Only fit if Control Region (CR)
  if (m_currCateName.Contains("CR")) {
    (*categoryWS->pdf("bkgPdf_"+m_currCateName)).fitTo(*categoryWS->data(obsDataName), Minos(RooArgSet(*nuisCateWS)), SumW2Error(kTRUE));
  }
  (*categoryWS->var("nBkg_"+m_currCateName)).setVal(obsData->sumEntries());
  
  //plotSingleCateFit(categoryWS, obsDataName, obsVarName);
  
  // Create Asimov data the old-fashioned way:
  //createAsimovData(categoryWS, 0, m_muNominalSH);
  //createAsimovData(categoryWS, 1, m_muNominalSH);
  
  // Print and return category workspace:
  std::cout << "DHWorkspace: Printing workspace for category:" << m_currCateName
	    << std::endl;
  categoryWS->Print("v");
  return categoryWS;
}

/**
   -----------------------------------------------------------------------------
   A private method for constructing a normalization systematic uncertainty.
   @param varName - the name of the systematic uncertainty.
   @param setup - the systematic uncertainty configuration parameters.
   @param nuisParams - the set of nuisance parameters to which this will add.
   @param constraints - the set of constraints to which this will add a term. 
   @param globalObs - the set of global observables, to which this will add.
   @param expected - the set of expected terms. 
   @returns - void. 
*/
void DHWorkspace::makeNP(TString varName, double setup[4],
			 RooArgSet *&nuisParams, RooArgSet *&constraints,
			 RooArgSet *&globalObs, RooArgSet *&expected) {
  
  // Settings for the systematic:
  double sigma    = setup[0];
  double sigmaLow = setup[1];
  double beta     = setup[2];
  double nominal  = setup[3];
  
  RooWorkspace* workspace = new RooWorkspace(varName);
  std::cout << "Creating a nuisance parameter " << varName << std::endl;

  // Create nuisance parameter with asymmetric uncertainty:
  if (sigmaLow > 0) {
    std::cout << "  parameter has an asymmetric uncertainty" << std::endl;
    
    RooRealVar* var = new RooRealVar(Form("nuisPar_%s",varName.Data()),
				     Form("nuisPar_%s",varName.Data()),
				     0,-5,5);
    RooRealVar* varBeta = new RooRealVar(Form("beta_%s",varName.Data()), 
					  Form("beta_%s",varName.Data()),
					  beta);
    RooProduct* varXBeta = new RooProduct(Form("%s_times_beta",varName.Data()),
					  Form("%s_times_beta",varName.Data()),
					  RooArgSet(*var,*varBeta));
    vector<double> sVarHi; sVarHi.clear(); sVarHi.push_back(1+sigma);
    vector<double> sVarLo; sVarLo.clear(); sVarLo.push_back(1-sigmaLow);
    RooArgList nuisList(*varXBeta);
    
    TString expName = Form("expected_%s",varName.Data());
    RooStats::HistFactory::FlexibleInterpVar expVar(expName, expName, nuisList,
						    nominal, sVarLo, sVarHi);
    workspace->import(expVar);
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varName.Data(), varName.Data(), varName.Data()));
  }
  
  // Create nuisance parameter with gaussian (symmetric) uncertainty:
  else if (sigmaLow == -999) {
    std::cout << "  parameter has a Gaussian constraint term" << std::endl;
    
    workspace->factory(Form("sum::expected_%s(nominal_%s[%f], prod::uncer_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]), sigma_%s[%f]))", varName.Data(), varName.Data(), nominal, varName.Data(), varName.Data(), varName.Data(), varName.Data(), beta, varName.Data(), sigma));
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varName.Data(), varName.Data(), varName.Data()));
  }
  
  // Create nuisance parameter with bifurcated Gaussian constraint term:
  else if (sigmaLow<0 && sigmaLow != -999) {
    
    TString valLogKappa = Form("%f",sqrt(log(1+pow(sigma,2))));
    TString valA = Form("%f",fabs(sigma/sigmaLow)); 
    
    std::cout << "  parameter has a Bif. Gauss constraint term" << std::endl;
    std::cout << "  asymmetric factor is " << valA << std::endl;
    
    workspace->factory(Form("valLogKappa_%s[%s]", varName.Data(), valLogKappa.Data()));
    workspace->factory(Form("RooExponential::expTerm_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]), valLogKappa_%s)", varName.Data(), varName.Data(), varName.Data(), varName.Data(), beta, varName.Data()));
    workspace->factory(Form("prod::expected_%s(expTerm_%s, nominal_%s[%f])", varName.Data(), varName.Data(), varName.Data(), nominal));
    workspace->factory(Form("RooBifurGauss::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1,%s)", varName.Data(), varName.Data(), varName.Data(), valA.Data()));
  }
  
  // Create a nuisance parameter with log-normal constraint term:
  else {
    std::cout << "  parameter has logNormal constraint term" << std::endl;
    TString valLogKappa = Form("%f",sqrt(log(1+pow(sigma,2))));
    
    workspace->factory(Form("valLogKappa_%s[%s]", varName.Data(), valLogKappa.Data()));
    workspace->factory(Form("RooExponential::expTerm_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]), valLogKappa_%s)", varName.Data(), varName.Data(), varName.Data(), varName.Data(), beta, varName.Data()));
    workspace->factory(Form("prod::expected_%s(expTerm_%s,nominal_%s[%f])", varName.Data(), varName.Data(), varName.Data(), nominal));
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varName.Data(), varName.Data(), varName.Data()));
  }
  
  // Add parameters and PDFs to relevant sets:
  nuisParams->add(*workspace->var(Form("nuisPar_%s",varName.Data())));
  constraints->add(*workspace->pdf(Form("constrPdf_%s",varName.Data())));
  globalObs->add(*workspace->var(Form("globOb_%s",varName.Data())));
  expected->add(*workspace->function(Form("expected_%s",varName.Data())));
}

/**
   -----------------------------------------------------------------------------
   A private method for constructing a shape systematic uncertainty. The shape
   NP maker will give variables used in parameterization process dependent name,
   but keep the same name for nuisance parameter, and global observables.
   @param varNameNP - the name of the systematic uncertainty.
   @param process - the corresponding signal process for the systematic.
   @param setup - the systematic uncertainty configuration parameters.
   @param nuisParams - the set of nuisance parameters to which this will add.
   @param constraints - the set of constraints to which this will add a term. 
   @param globalObs - the set of global observables, to which this will add.
   @param expected - the set of expected terms. 
   @returns - void. 
*/
void DHWorkspace::makeShapeNP(TString varNameNP, TString process,
			      double setup[4], RooArgSet *&nuisParams,
			      RooArgSet *&constraints, RooArgSet *&globalObs,
			      RooArgSet *&expected) {
  
  // Settings for the systematic:
  double sigma    = setup[0];
  double sigmaLow = setup[1];
  double beta     = setup[2];
  double nominal  = setup[3];
  TString varName = varNameNP + process;
  
  RooWorkspace* workspace = new RooWorkspace(varName);
  
  // Create nuisance parameter with asymmetric uncertainty:
  if (sigmaLow > 0) {
    std::cout << "  parameter for an asymmetric uncertainty" << std::endl;
    
    RooRealVar* var = new RooRealVar(Form("nuisPar_%s",varNameNP.Data()),
				     Form("nuisPar_%s",varNameNP.Data()),
				     0, -5, 5);
    RooRealVar* varBeta = new RooRealVar(Form("beta_%s",varName.Data()),
					  Form("beta_%s",varName.Data()),
					  beta);
    RooProduct* varXBeta = new RooProduct(Form("%s_times_beta",varName.Data()),
					  Form("%s_times_beta",varName.Data()),
					  RooArgSet(*var,*varBeta));
    
    vector<double> sVarHi; sVarHi.clear(); sVarHi.push_back(1+sigma);
    vector<double> sVarLo; sVarLo.clear(); sVarLo.push_back(1-sigmaLow);
    RooArgList nuisList(*varXBeta);
    RooStats::HistFactory::FlexibleInterpVar expVar(Form("expected_%s",varName.Data()), Form("expected_%s",varName.Data()), nuisList, nominal, sVarLo, sVarHi);
    
    workspace->import(expVar);
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varNameNP.Data(), varNameNP.Data(), varNameNP.Data()));
  }
  
  // Create a nuisance parameter with Gaussian constraint term:
  else if (sigmaLow == -999) {
    std::cout << "  parameter with a Gaussian constraint term" << std::endl;
    
    workspace->factory(Form("sum::expected_%s(nominal_%s[%f], prod::uncer_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]), sigma_%s[%f]))", varName.Data(), varName.Data(), nominal, varName.Data(), varName.Data(), varNameNP.Data(), varName.Data(), beta, varName.Data(), sigma));
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varNameNP.Data(), varNameNP.Data(), varNameNP.Data()));
  }
  
  // Create a nuisance parameter with log-normal constraint term:
  else {
    TString valLogKappa = Form("%f", sqrt( log( 1+pow(sigma,2)) ) );
    workspace->factory(Form("valLogKappa_%s[%s]", varName.Data(), valLogKappa.Data()));
    workspace->factory(Form("RooExponential::expTerm_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]),valLogKappa_%s)", varName.Data(), varName.Data(), varNameNP.Data(), varName.Data(), beta, varName.Data()));
    workspace->factory(Form("prod::expected_%s(expTerm_%s,nominal_%s[%f])", varName.Data(), varName.Data(), varName.Data(), nominal));
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varNameNP.Data(), varNameNP.Data(), varNameNP.Data()));
  }
  
  // Add parameters and PDFs to relevant sets:
  nuisParams->add(*workspace->var(Form("nuisPar_%s",varNameNP.Data())));
  constraints->add(*workspace->pdf(Form("constrPdf_%s",varNameNP.Data())));
  globalObs->add(*workspace->var(Form("globOb_%s",varNameNP.Data())));
  expected->add(*workspace->function(Form("expected_%s",varName.Data())));
}

/*
   -----------------------------------------------------------------------------
   Plot the data and fit in a single category

void DHWorkspace::plotSingleCateFit(RooWorkspace *cateWS, TString dataset, 
				    TString observableName) {
  std::cout << "DMWorkspace: Plot single category fit for "
	    << m_currCateName << std::endl;
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  RooPlot* frame = (*cateWS->var(observableName+"_"+m_currCateName)).frame(50);
  cateWS->data(dataset)->plotOn(frame);
  (*cateWS->pdf("model_"+m_currCateName)).plotOn(frame, LineColor(2));
  (*cateWS->pdf("model_"+m_currCateName)).plotOn(frame, Components((*cateWS->pdf("bkgPdf_"+m_currCateName))), LineColor(4), LineStyle(2));
  (*cateWS->pdf("model_"+m_currCateName)).plotOn(frame, Components((*cateWS->pdf("sigPdfDH_"+m_currCateName))), LineColor(3), LineStyle(2));
  (*cateWS->pdf("model_"+m_currCateName)).plotOn(frame, Components((*cateWS->pdf("sigPdfSH_"+m_currCateName))), LineColor(6), LineStyle(2));
  
  //double chi2 = frame->chiSquare();
  frame->SetYTitle("Events / GeV");
  frame->SetXTitle("Mass [GeV]");
  frame->Draw();
  
  TLatex text; text.SetNDC(); text.SetTextColor(1);
  text.DrawLatex(0.2, 0.81, Form("Category %d", m_currCateIndex));
  TH1F *histSH = new TH1F("histSH", "histSH", 1, 0, 1);
  TH1F *histDH = new TH1F("histDH", "histDH", 1, 0, 1);
  TH1F *histNR = new TH1F("histNR", "histNR", 1, 0, 1);
  TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
  histDH->SetLineColor(3);
  histSH->SetLineColor(6);
  histNR->SetLineColor(4);
  histSig->SetLineColor(2);
  TLegend leg(0.61, 0.6, 0.89, 0.77);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  leg.SetBorderSize(0);
  leg.AddEntry(histSig, "Sig. + bkg.", "l");
  leg.AddEntry(histDH, "Di-Higgs", "l");
  leg.AddEntry(histSH, "Single Higgs", "l");
  leg.AddEntry(histNR, "Non-resonant", "l");
  leg.Draw("SAME");
  can->Print(Form("%s/Plots/cateFit_%s.eps",m_outputDir.Data(),dataset.Data()));
  delete can;
}
*/


/**
   -----------------------------------------------------------------------------
   Create Asimov data for the statistical model, using a fit to observed data
   for the shape and normalizaiton of the background.
   @param cateWS - the current category workspace.
   @param valMuDH - the value of the di-Higgs signal strength to use.
   @param valMuSH - the value of the single-Higgs signal strength to use.
*/
void DHWorkspace::createAsimovData(RooWorkspace* cateWS, int valMuDH,
				   int valMuSH) {
  std::cout << "DHWorkspace: Creating Asimov data, mu=" << valMuDH << std::endl;
  
  // Set mu_DM and mu_SM to the specified values:
  RooRealVar *poi = NULL;
  double initialMuDH = 0.0;
  double initialMuSH = 0.0;
  if (cateWS->var("mu_DH")) {
    poi = cateWS->var("mu_DH");
    initialMuDH = poi->getVal();
    poi->setVal(valMuDH);
    poi->setConstant(true);  
  }
  if (cateWS->var("mu_SH")) {
    initialMuSH = cateWS->var("mu_SH")->getVal();
    cateWS->var("mu_SH")->setVal(valMuSH);
    cateWS->var("mu_SH")->setConstant(true);  
  }
    
  RooDataSet *asimov = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*cateWS->pdf(Form("model_%s", m_currCateName.Data())), *cateWS->set(Form("observables_%s",m_currCateName.Data())));
  asimov->SetNameTitle(Form("asimovDataMu%d_%s",valMuDH,m_currCateName.Data()),
		       Form("asimovDataMu%d_%s",valMuDH,m_currCateName.Data()));
  cateWS->import(*asimov);
  
  if (cateWS->var("mu_DH")) {
    cateWS->var("mu_DH")->setVal(initialMuDH);
    cateWS->var("mu_DH")->setConstant(false);
  }
  if (cateWS->var("mu_SH")) {
    cateWS->var("mu_SH")->setVal(initialMuSH);
    cateWS->var("mu_SH")->setConstant(true);
  }
  
  std::cout << "DHWorkspace: Asimov data has " << asimov->sumEntries() 
	    << " entries" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Create Asimov data for the statistical model, using a fit to observed data
   for the shape and normalizaiton of the background.
   @param cateWS - the current category workspace.
   @param valMuDH - the value of the di-Higgs signal strength to use.
   @param valMuSH - the value of the single-Higgs signal strength to use.
*/
RooDataSet* DHWorkspace::createAsimovData(int valMuDH, int valMuSH) {
  std::cout << "DHWorkspace: Creating Asimov data, mu=" << valMuDH << std::endl;
  
  // Set mu_DM and mu_SM to the specified values:
  RooRealVar *poi = NULL;
  double initialMuDH = 0.0;
  double initialMuSH = 0.0;
  poi = m_combinedWS->var("mu_DH");
  initialMuDH = poi->getVal();
  poi->setVal(valMuDH);
  poi->setConstant(true);  
  initialMuSH = m_combinedWS->var("mu_SH")->getVal();
  m_combinedWS->var("mu_SH")->setVal(valMuSH);
  m_combinedWS->var("mu_SH")->setConstant(true);  
  
  TString asimovName = Form("asimovDataMu%d_%s",valMuDH,m_currCateName.Data());
  //RooDataSet *asimov = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*m_combinedWS->pdf(Form("model_%s", m_currCateName.Data())), *m_combinedWS->set(Form("observables_%s",m_currCateName.Data())));
  RooDataSet *asimov = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*m_combinedWS->pdf(Form("model_%s",m_currCateName.Data())), *m_combinedWS->set("observables"));
  asimov->SetNameTitle(asimovName, asimovName);
  m_combinedWS->import(*asimov);
  
  m_combinedWS->var("mu_DH")->setVal(initialMuDH);
  m_combinedWS->var("mu_DH")->setConstant(false);
  m_combinedWS->var("mu_SH")->setVal(initialMuSH);
  m_combinedWS->var("mu_SH")->setConstant(true);
  
  std::cout << "DHWorkspace: Asimov data has " << asimov->sumEntries() 
	    << " entries" << std::endl;
  return (RooDataSet*)(m_combinedWS->data(asimovName));
}
