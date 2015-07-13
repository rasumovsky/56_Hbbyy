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
//    - "ResOnly": only generate the resonant search workspace.               //
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
   @param jobName - The name of the job
   @param cateScheme - The name of the event categorization
   @param options - The job options ("New", "FromFile"), etc.
   @returns void
*/
DHWorkspace::DHWorkspace(TString jobName, TString cateScheme,
			 TString options) {
  m_jobName = jobName;
  m_cateScheme = cateScheme;
  m_options = options;
  m_allGoodFits = true;
  
  m_combinedWS = NULL;
  m_mConfig.clear();
  
  std::cout << "\nDHWorkspace: Initializing..." << "\n\tjobName = " << m_jobName
	    << "\n\tcateScheme = " << m_cateScheme << "\n\toptions = "
	    << m_options << std::endl;
  
  // Assign output directory, and make sure it exists:
  m_outputDir = Form("%s/%s/DHWorkspace", DHAnalysis::masterOutput.Data(),
		     m_jobName.Data());
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
   @param analysisType - The type of analysis (Res, NonRes).
*/
ModelConfig* DHWorkspace::getModelConfig(TString analysisType) {
  return m_mConfig[analysisType];
}

/**
   -----------------------------------------------------------------------------
   Load a previously created workspace.
*/
void DHWorkspace::loadWSFromFile() {
  //Check to see if the workspace has actually been made.
  TFile inputFile(Form("%s/rootfiles/workspaceDH.root", m_outputDir.Data()),
		  "read");
  if (inputFile.IsOpen()) {
    std::cout << "DHWorkspace: Loading workspace from file..."<< std::endl;
    // Load the single workspace file:
    m_combinedWS = (RooWorkspace*)inputFile.Get("combinedWS");
    
    // Load the multiple models contained in the workspace:
    for (int i_t = 0; i_t < DHAnalysis::nAnalysisTypes; i_t++) {
      currAna = DHAnalysis::analysisTypes[i_t];
      if ((m_combinedWS->obj(Form("modelConfig_%s",currAna.Data())))) {
	m_mConfig[currAna] = (ModelConfig*)m_combinedWS
	  ->obj(Form("modelConfig_%s",currAna.Data()));
      }	
      else {
	std::cout << "DHWorkspace: could not retrieve model for " << currAna
		  << std::endl;
      }
    }
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
  
  // The combined workspace:
  m_combinedWS = new RooWorkspace("combinedWS");
  m_combinedWS->importClassCode();
  
  std::cout << "DHWorkspace: Create a new workspace from scratch." << std::endl;
  std::cout << "\n........................................" << std::endl;
  std::cout << "Luminosity at 13 TeV: " << DHAnalysis::analysisLuminosity 
	    << " pb-1." << std::endl;
  
  // Loop over the analysis types and define a statistical model for each:
  for (int i_t = 0; i_t < DHAnalysis::nAnalysisTypes; i_t++) {
    
    // Get the current analysis name:
    currAna = DHAnalysis::analysisTypes[i_t];
    currNCategories = DHAnalysis::getNumCategories(m_cateScheme, currAna); 
    
    // Check to see if the analysis should be implemented:
    if ((m_options.Contains("ResOnly") && currAna.EqualTo("NonRes")) ||
	(m_options.Contains("NonResOnly") && currAna.EqualTo("Res"))) continue;
    
    // Go ahead and built that analysis model!
    createNewModel();
  }

  // Write workspace to file:
  m_combinedWS->writeToFile(Form("%s/rootfiles/workspaceDH.root",
				 m_outputDir.Data()));
}

/**
   -----------------------------------------------------------------------------
   Create one of the analysis models from scratch. Called from createNewWS().
*/
void DHWorkspace::createNewModel() {
  std::cout << "DHWorkspace: Creating model for " << currAna << " analysis." 
	    << std::endl;
  // Define and name analysis categories:
  std::cout << "  Number of categories = " << currNCategories << std::endl;
  

  //--------------------------------------//
  // Initialize classes relevant to workspace:
  
  // Category workspaces, and RooCategory and Simultaneous PDF:
  RooWorkspace* cateWS[currNCategories];
  RooCategory *categories
    = new RooCategory(Form("categories_%s", currAna.Data()),
		      Form("categories_%s", currAna.Data()));
  RooSimultaneous *combinedPdf
    = new RooSimultaneous(Form("combinedPdf_%s", currAna.Data()),
			  Form("combinedPdf_%s", currAna.Data()), *categories);
  
  // Instantiate parameter sets:
  RooArgSet *nuisanceParameters = new RooArgSet();
  RooArgSet *muSHConstants = new RooArgSet();
  RooArgSet *globalObservables = new RooArgSet();
  RooArgSet *observables = new RooArgSet();
  RooArgSet *constraints = new RooArgSet();
  
  // Maps for datasets:
  std::map<std::string, RooDataSet*> dm;
  
  // Read tables of PES and PER and store values:
  pes = new PESReader(DHAnalysis::fileNamePESValues, currNCategories);
  per = new PERReader(DHAnalysis::fileNamePERValues, currNCategories);

  //--------------------------------------//
  // Loop over channels:
  std::cout << "DHWorkspace: Looping over " << currAna << " categories."
	    << std::endl;
  for (currCateIndex = 0; currCateIndex < currNCategories; currCateIndex++) {
    
    currCateName = DHAnalysis::cateIndexToName(m_cateScheme, currAna, 
					       currCateIndex);
    
    // Create the workspace for a single category:
    cateWS[currCateIndex] = createNewCategoryWS();
    categories->defineType(currCateName);
    
    // Add PDFs and parameters: 
    TString namePdf = Form("model_%s", currCateName.Data());
    TString nameNP = Form("nuisanceParameters_%s", currCateName.Data());
    TString nameGlob = Form("globalObservables_%s", currCateName.Data());
    TString nameMuC = Form("muConstants_%s", currCateName.Data());
    TString nameObs = Form("observables_%s", currCateName.Data());
    combinedPdf->addPdf(*cateWS[currCateIndex]->pdf(namePdf), currCateName);
    nuisanceParameters->add(*cateWS[currCateIndex]->set(nameNP));
    globalObservables->add(*cateWS[currCateIndex]->set(nameGlob));
    muSHConstants->add(*cateWS[currCateIndex]->set(nameMuC));
    nuisanceParameters->add(*cateWS[currCateIndex]->set(nameMuC));
    observables->add(*cateWS[currCateIndex]->set(nameObs));
    
    // Add category datasets to combined workspace and combined data maps:
    TString nameData = Form("obsData_%s", currCateName.Data());
    m_combinedWS->import(*(RooDataSet*)cateWS[currCateIndex]->data(nameData));
    dm[(string)currCateName] = (RooDataSet*)m_combinedWS->data(nameData);
  }
  std::cout << "DHWorkspace: Beginning to combine all categories." << std::endl;
  
  // Define the combined dataset:
  RooRealVar wt("wt", "wt", 1);
  RooArgSet *args = new RooArgSet();
  args->add(*observables);
  args->add(wt);
  RooDataSet* obsData = new RooDataSet(Form("obsData_%s", currAna.Data()),
				       Form("obsData_%s", currAna.Data()),
				       *args, Index(*categories), Import(dm),
				       WeightVar(wt));
  
  // Import PDFs, parameters, and dataset into workspace:
  m_combinedWS->import(*categories);
  m_combinedWS->import(*combinedPdf);
  m_combinedWS->defineSet(Form("nuisanceParameters_%s", currAna.Data()),
			  *nuisanceParameters);
  m_combinedWS->defineSet(Form("observables_%s", currAna.Data()),
			  *observables);
  m_combinedWS->defineSet(Form("globalObservables_%s", currAna.Data()),
			  *globalObservables);
  m_combinedWS->defineSet(Form("poi_%s",currAna.Data()),
			  RooArgSet(*m_combinedWS->var("mu_DH")));   
  m_combinedWS->defineSet(Form("muSHConstants_%s", currAna.Data()),
			  *muSHConstants);
  m_combinedWS->import(*obsData);
  
  // Define the ModelConfig for the analysis and import to the workspace:
  m_mConfig[currAna] = new ModelConfig(Form("modelConfig_%s",currAna.Data()),
				       m_combinedWS);
  m_mConfig[currAna]->SetPdf((*m_combinedWS->pdf(Form("combinedPdf_%s",currAna.Data()))));
  m_mConfig[currAna]->SetObservables((*m_combinedWS->set(Form("observables_%s",currAna.Data()))));
  m_mConfig[currAna]->SetParametersOfInterest((*m_combinedWS->set(Form("poi_%s",currAna.Data()))));
  m_mConfig[currAna]->SetNuisanceParameters((*m_combinedWS->set(Form("nuisanceParameters_%s",currAna.Data()))));
  m_mConfig[currAna]->SetGlobalObservables((*m_combinedWS->set(Form("globalObservables_%s",currAna.Data()))));
  m_combinedWS->import(*m_mConfig[currAna]);
  
  std::cout << "DHWorkspace: Printing the " << currAna << " combined workspace."
	    << std::endl;
  m_combinedWS->Print("v");
  
  // Start profiling the data:
  std::cout << "DHWorkspace: Start profiling data" << std::endl;
  
  // Save a snapshot of the parameters:
  m_combinedWS->saveSnapshot(Form("paramsOrigin_%s",currAna.Data()),
			     *nuisanceParameters);
  
  // Choose two models for the default test fits:
  TString currDHSignal = currAna.Contains("NonRes") ?
    DHAnalysis::sigDHModes[0] : DHAnalysis::sigDHModes[1];
  
  DHTestStat *dhts = new DHTestStat(m_jobName, currDHSignal, m_cateScheme,
				    "new", m_combinedWS);
  dhts->saveSnapshots(true);
  dhts->setPlotDirectory(Form("%s/Plots/", m_outputDir.Data()));
  
  m_dataToPlot = (DHAnalysis::doBlind) ?
    Form("asimovDataMu1_%s",currAna.Data()) : Form("obsData_%s",currAna.Data());

  double profiledMuDHVal = -999.0;
  TString currDataName = Form("%s_%s", m_dataToPlot.Data(), currAna.Data());
  // Mu = 0 fits:
  double nllMu0 = dhts->getFitNLL(currDataName, 0, true, profiledMuDHVal);
  if (!dhts->fitsAllConverged()) m_allGoodFits = false;
  // Mu = 1 fits:
  double nllMu1 = dhts->getFitNLL(currDataName, 1, true, profiledMuDHVal);
  if (!dhts->fitsAllConverged()) m_allGoodFits = false;
  // Mu free fits:
  double nllMuFree = dhts->getFitNLL(currDataName, 1, false, profiledMuDHVal);
  if (!dhts->fitsAllConverged()) m_allGoodFits = false;
  
  // Calculate profile likelihood ratios:
  double llrL1L0 = nllMu1 - nllMu0;
  double llrL1Lfree = profiledMuDHVal > 1.0 ? 0.0 : (nllMu1 - nllMuFree);
  double llrL0Lfree = profiledMuDHVal < 0.0 ? 0.0 : (nllMu0 - nllMuFree);
  
  // Print summary of the fits:
  std::cout.precision(10);
  std::cout << "\nDHWorkspace: Printing likelihood results for " << currAna
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
  fileMuProf.open(Form("%s/mu/mu_%s.txt",
		       m_outputDir.Data(), currDHSignal.Data()));
  fileMuProf << profiledMuDHVal << std::endl;
  fileMuProf.close();
}

/**
   -----------------------------------------------------------------------------
   Create the workspace for a single analysis category.
   @param currCategory
*/
RooWorkspace* DHWorkspace::createNewCategoryWS() {
    
  // The bools that control the systematic uncertainties:
  bool inclusive = (currCateName == "inclusive");
  bool channel_constraints_attached = (currCateIndex == 0);
  bool m_norm = !m_options.Contains("nonorm");
  bool m_pes = !m_options.Contains("nopes");
  bool m_per = !m_options.Contains("noper");
  bool m_ss  = !m_options.Contains("noss");
  bool m_bgm = !m_options.Contains("nobgm");
  bool m_mig = !m_options.Contains("nomig");
  bool m_nosys = m_options.Contains("nosys");
  if (m_nosys) {
    std::cout << "\tDHWorkspace: ALL systematics OFF" << endl;
    m_norm = false;   m_pes = false;   m_per = false;
    m_ss = false;     m_bgm = false;   m_mig = false;
  }
  std::cout << "\tNormalization systematics = " << m_norm << std::endl;
  std::cout << "\tEnergy scale systematics  = " << m_pes  << std::endl;
  std::cout << "\tResolution systematics    = " << m_per  << std::endl;
  std::cout << "\tShape systematics         = " << m_ss   << std::endl;
  std::cout << "\tBackground systematics    = " << m_bgm  << std::endl;
  std::cout << "\tMigration systematics     = " << m_mig  << std::endl;
  
  //--------------------------------------//
  // Create the individual channel workspace:
  RooWorkspace *tempWS = new RooWorkspace(Form("tmpWS_%s",currCateName.Data()));
  
  // nuispara:
  RooArgSet *nuisParams = new RooArgSet();
  RooArgSet *nuisParamsBkg = new RooArgSet();
  RooArgSet *nuisParamsUncorrelated = new RooArgSet();
  // constraints:
  RooArgSet *constraints = new RooArgSet();
  RooArgSet *constraintsBias = new RooArgSet();
  // globobs:
  RooArgSet *globalObs = new RooArgSet();
  RooArgSet *globalObsProc = new RooArgSet();
  // expected:
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
  if (m_norm) {
    double setupLumi[4] = {0.036, 0, 1, 1};
    makeNP("Luminosity", setupLumi, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
    double setupTrigger[4] = {0.005, 0, 1, 1};
    makeNP("Trigger", setupTrigger, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
    double setupIsEM[4] = {0.0526, 0, 1, 1};
    makeNP("PhotonID", setupIsEM, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
    double setupIso[4] = {0.004, 0, 1, 1};
    makeNP("Isolation", setupIso, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
    double setupESCALE[4] = {0.003, 0, 1, 1};
    makeNP("ESCALE", setupESCALE, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
  }
  
  //--------------------------------------//
  // Migration systematics:
  if (m_mig) {
    // Follow examples in other workspaces...
  }
  
  //--------------------------------------//
  // SYSTEMATICS: Spurious signal
  if (m_bgm) {
    double ssEvents = 0.1;//spuriousSignal function;
    double setupBias[4] = {ssEvents, -999, 1, 0}; //Gaussian constraint
    makeNP("bias", setupBias, *&nuisParamsUncorrelated, *&constraintsBias,
	   *&globalObs, *&expectedBias);
    RooProduct sigBias("sigBias","sigBias",*expectedBias);
    tempWS->import(sigBias);
  }
  else tempWS->factory("sigBias[0]");//expectedBias
  
  //--------------------------------------//
  // SYSTEMATICS - Resolution:
  TString m_listMRS = "";
  std::vector<TString> perList; perList.clear();
  if (m_per) {
    double setupPER[4] = {0.0, 0, 1, 1};
    // Loop over sources of resolution systematic uncertainty:
    for (int i_s = 0; i_s < per->getNumberOfSources(); i_s++) {
      TString currPERSource = per->getNameOfSource(i_s);
      TString currPERName = Form("EM_%s",currPERSource.Data());
      perList.push_back(currPERName);
      setupPER[0] = per->getValue(currPERSource, currCateIndex);
      setupPER[2] = per->getSign(currPERSource, currCateIndex);
      // Resolution on the inclusive shape:
      makeShapeNP(currPERName, "DH", setupPER, *&nuisParams, *&constraints,
		  *&globalObs, *&expectedShape);
      makeShapeNP(currPERName, "SH", setupPER, *&nuisParams, *&constraints,
		  *&globalObs, *&expectedShape);
    }
  }
  
  //--------------------------------------//
  // SYSTEMATICS - Energy-scale:
  TString m_listMSS = "";
  std::vector<TString> pesList; pesList.clear();
  if (m_pes) {
    double setupPES[4] = {0.0, 0, 1, 1};
    // loop over sources of energy scale systematic uncertainty:
    for (int i_s = 0; i_s < pes->getNumberOfSources(); i_s++) {
      TString currPESSource = pes->getNameOfSource(i_s);
      TString currPESName = Form("EM_%s",currPESSource.Data());
      pesList.push_back(currPESName);
      setupPES[0] = pes->getValue(currPESSource, currCateIndex);
      setupPES[2] = pes->getSign(currPESSource, currCateIndex);
      makeNP(currPESName, setupPES, *&nuisParams, *&constraints, *&globalObs,
	     *&expectedShape);
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
  RooProduct expectationDH("expectationDH","expectationDH", *expectedDH);
  RooProduct expectationSH("expectationSH","expectationSH", *expectedSH);
  RooProduct expectationCommon("expectationCommon","expectationCommon",
			       *expected);
  
  // Spurious signal term will assume the shape of "inclusive" pdf.
  tempWS->import(expectationCommon);
  tempWS->import(expectationSH);
  tempWS->import(expectationDH);
  tempWS->import(*expectedShape);
  tempWS->import(*expectedBias);
  
  // Declare the observable, and the observables set:
  TString obsName = currAna.EqualTo("NonRes") ? "m_yy" : "m_bbyy";
  double obsMin = currAna.EqualTo("NonRes") ?
    DHAnalysis::DHMyyRangeLo : DHAnalysis::DHMyybbRangeLo;
  double obsMax = currAna.EqualTo("NonRes") ? 
    DHAnalysis::DHMyyRangeHi : DHAnalysis::DHMyybbRangeHi;
  tempWS->factory(Form("%s[%f,%f]", obsName.Data(), obsMin, obsMax));
  tempWS->defineSet("obsprelim", obsName);
  
  //--------------------------------------//
  // Begin PDF construction.
  
  // Construct the background PDF (present in every category):
  std::cout << "DHWorkspace: Constructing background model" << std::endl;
  BkgModel *currBkgModel = new BkgModel(tempWS->var(Form("%s",obsName.Data())));
  TString bkgFuncName = DHAnalysis::cateToBkgFunc(currCateName);
  currBkgModel->addBkgToCateWS(tempWS, nuisParamsBkg, bkgFuncName);
  
  // Add bkg params to correlated collection (SR and CR have same shape):
  nuisParams->add(*nuisParamsBkg);
  
  // Start to construct the statistical model, add components one by one:
  TString modelForm = "SUM::modelSB(nBkg*bkgPdf";
  
  // Construct SH and DH signals (depending on categories):
  std::cout << "DHWorkspace: Constructing signal model" << std::endl;
  if (currAna.EqualTo("NonRes")) {
    
    // SH (Single-Higgs) signal:
    if (DHAnalysis::cateHasComponent(currCateName, "SH")) {
      // jj CR first, bb SR second
      double normSH[2] = {0.00868, 0.186};
      tempWS->factory(Form("RooCBShape::pdfCB_SH(%s, prod::muCB_SH(muCBNom_SH[124.96]%s), prod::sigmaCB_SH(sigmaCBNom_SH[1.53]%s), alphaCB_SH[1.56], nCB_SH[10.0])", obsName.Data(), m_listMSS.Data(), m_listMRS.Data()));
      tempWS->factory(Form("RooGaussian::pdfGA_SH(%s, prod::muGA_SH(muCBNom_SH%s), prod::sigmaGA_SH(sigmaGANom_SH[23.98]%s))", obsName.Data(), m_listMSS.Data(), m_listMRS.Data()));
      tempWS->factory(Form("SUM::sigPdf_SH(fracCB_SH[1.00]*pdfCB_SH,pdfGA_SH)"));
      tempWS->factory(Form("nSH[%f]", normSH[currCateIndex]));
      
      // Normalization for SH (expectationCommon = mu*isEM*lumi*migr):
      tempWS->factory("prod::nSigSH(nSH,expectationCommon,expectationSH)");
      
      // Add to the statistical model:
      modelForm += ",nSigSH*sigPdfSH";
    }
    
    // Now BSM (Di-Higgs) signal:
    if (DHAnalysis::cateHasComponent(currCateName, "DH")) {
      tempWS->factory(Form("RooCBShape::pdfCB_DH(%s, prod::muCB_DH(muCBNom_DH[124.95]%s), prod::sigmaCB_DH(sigmaCBNom_DH[1.31]%s), alphaCB_DH[1.56], nCB_DH[10.0])", obsName.Data(), m_listMSS.Data(), m_listMRS.Data()));
      tempWS->factory(Form("RooGaussian::pdfGA_DH(%s, prod::muGA_DH(muCBNom_DH%s), prod::sigmaGA_DH(sigmaGANom_DH[2.84]%s))", obsName.Data(), m_listMSS.Data(), m_listMRS.Data()));
      tempWS->factory(Form("SUM::sigPdf_DH(fracCB_DH[0.94]*pdfCB_DH,pdfGA_DH)"));
      tempWS->factory("nDH[5,0,100]");
      
      // Normalization for DH (expectationCommon = mu*isEM*lumi*migr):
      tempWS->factory("prod::nSigDH(nDH,expectationCommon,expectationDH)");
      
      // Add to the statistical model:
      modelForm += ",nSigDH*sigPdfDH,sigBias*sigPdfDH";
    }
  }
  else if (currAna.EqualTo("Res")) {
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
  if (currCateIndex == 0) {
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
    also share the same name. nuisParams should contain all correlated nuisance
    parameters. All uncorrelated nuisance parameters should be included in
    nuisParamsUncorrelated.
  */
  TString corrNPNames = "mu_DH,mu_SH";
  
  // Iterate over nuisance parameters:
  TIterator *iterNuis = nuisParams->createIterator();
  RooRealVar* currNuis;
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    std::cout << "\t" << currNuis->GetName() << std::endl;
    corrNPNames += Form(",nuisPar_%s,globOb_%s",currNuis->GetName(),
			currNuis->GetName());
  }
  std::cout << "For category " << currCateName << ", correlate variables: "
	    << corrNPNames << std::endl;
  
  /*
    Sub-channel labeling
    Import the workspace tempWS to another workspace and add currCateName as a 
    suffix to all nodes and variables of w. the correlated nuisance parameters
    and their respective global observables will not be renamed.
  */
  RooWorkspace* categoryWS = new RooWorkspace("workspace_"+currCateName);
  categoryWS->import((*tempWS->pdf("model")), RenameAllNodes(currCateName),
		     RenameAllVariablesExcept(currCateName,corrNPNames),
		     Silence());
  
  // Adding correlated nuisance parameters to nuisanceParameters:
  RooArgSet* nuisCateWS = new RooArgSet();
  iterNuis->Reset();
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    nuisCateWS->add(*(RooRealVar*)categoryWS->obj(currNuis->GetName()));
  }
  
  // Adding uncorrelated nuisance parameters to nuisanceParameters:
  TIterator *iterNuisUncorrelated = nuisParamsUncorrelated->createIterator();
  RooRealVar* currNuisUncorrelated;
  while ((currNuisUncorrelated = (RooRealVar*)iterNuisUncorrelated->Next())) {
    TString nuisName = (currNuisUncorrelated->GetName() 
			+ (TString)"_" + currCateName);
    nuisCateWS->add(*(RooRealVar*)categoryWS->obj(nuisName));
  }
  
  // Adding unconstrained NPs from the background pdf:
  RooArgSet* nuisBkgCateWS = new RooArgSet();
  TIterator *iterNuisBkg = nuisParamsBkg->createIterator();
  RooRealVar* currNuisBkg;
  while ((currNuisBkg = (RooRealVar*)iterNuisBkg->Next())) {
    TString parName = currNuisBkg->GetName()+(TString)"_"+currCateName;
    nuisBkgCateWS->add(*categoryWS->var(parName));
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
    TString globName = currGlobs->GetName()+(TString)"_"+currCateName;
    if ((bool)categoryWS->obj(globName)) {
      globsCateWS->add(*(RooRealVar*)categoryWS->obj(globName));
      categoryWS->var(globName)->setConstant();
    }
    else if ((bool)categoryWS->obj(currGlobs->GetName())) {
      globsCateWS->add(*(RooRealVar*)categoryWS->obj(currGlobs->GetName()));
      categoryWS->var(currGlobs->GetName())->setConstant();
    }
  }
  
  RooArgSet *obsCateWS = new RooArgSet();
  TIterator *iterObs = tempWS->set("obsprelim")->createIterator();
  RooRealVar *currObs;
  while ((currObs = (RooRealVar*)iterObs->Next())) {
    TString obsName = currObs->GetName()+(TString)"_"+currCateName;
    if ((bool)categoryWS->obj(obsName)) {
      obsCateWS->add(*(RooRealVar*)categoryWS->obj(obsName));
    }
    else {
      obsCateWS->add(*(RooRealVar*)categoryWS->obj(currObs->GetName()));
    }
  }
  
  // Set some of the mu values constant:
  std::cout << "DHWorkspace: Setting SH signals constant." << std::endl;
  RooArgSet* muConstCateWS = new RooArgSet();
  muConstCateWS->add(*categoryWS->var("mu_SH"));

  TIterator *iterMuConst = muConstCateWS->createIterator();
  RooRealVar *currMuConst;
  while ((currMuConst = (RooRealVar*)iterMuConst->Next())) {
    currMuConst->setVal(1.0);
    currMuConst->setConstant(true);
  }
  
  categoryWS->defineSet(Form("muConstants_%s",currCateName.Data()),
			*muConstCateWS);
  categoryWS->defineSet(Form("observables_%s",currCateName.Data()),
			*obsCateWS);
  categoryWS->defineSet(Form("nuisanceParameters_%s",currCateName.Data()),
			*nuisCateWS);
  categoryWS->defineSet(Form("globalObservables_%s",currCateName.Data()),
			*globsCateWS);
  
  //--------------------------------------//
  // Import the observed data set, and create a binned version:
  DHDataReader *dhdr = new DHDataReader(categoryWS->var(obsName));
  RooDataSet *obsData = NULL;
  if (currAna.EqualTo("NonRes")) {
    obsData = dhdr->loadNonResData(currCateName);
  }
  else {
    std::cout << "DHWorkspace: Error importing resonant data." << std::endl;
  }
  
  // Rename the imported dataset for consistency:
  TString obsDataName = Form("obsData_%s", currCateName.Data());
  obsData->SetNameTitle(obsDataName, obsDataName);
  
  // Set the background normalization parameter:
  std::cout << "DHWorkspace: Fitting bkg. in " << currCateName << std::endl;
  (*categoryWS->var("nBkg_"+currCateName)).setVal(obsData->sumEntries());
  (*categoryWS->pdf("bkgPdf_"+currCateName)).fitTo(*obsData, Minos(RooArgSet(*nuisBkgCateWS)), SumW2Error(kTRUE));//should be false?
  (*categoryWS->var("nBkg_"+currCateName)).setVal(obsData->sumEntries());
  
  categoryWS->import(*obsData);
  
  std::cout << "DHWorkspace: Printing workspace for category:" << currCateName
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
