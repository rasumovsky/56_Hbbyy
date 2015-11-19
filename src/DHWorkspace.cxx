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
//    - "NonRes", "Res"                                                       //
//                                                                            //
//  options:                                                                  //
//    - "ResonantOnly": only generate the resonant search workspace.          //
//    - "NonResOnly": only generate the nonresonant search workspace.         //
//                                                                            //
//  Job options: "New", "FromFile" determine whether to create a new workspace//
//  or load a previously generated one.                                       //
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
   @return void
*/
DHWorkspace::DHWorkspace(TString newConfigFile, TString newAnalysisType,
			 TString newOptions) {
  m_configFile = newConfigFile;
  m_anaType = newAnalysisType;
  m_options = newOptions;
  m_allGoodFits = true;
  
  m_ws = NULL;
  m_modelConfig = NULL;
    
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
   @return - true iff the fits all converged.
*/
bool DHWorkspace::fitsAllConverged() {
  return m_allGoodFits;
}

/**
   -----------------------------------------------------------------------------
   Retrieves the workspace created by this program.
*/
RooWorkspace* DHWorkspace::getCombinedWorkspace() {
  return m_ws;
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
    m_ws = (RooWorkspace*)inputFile.Get("combinedWS");
    m_modelConfig = (ModelConfig*)m_ws->obj("modelConfig");
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
  m_nCategories = cateNames.size();
  std::cout << "  Number of categories = " << m_nCategories << std::endl;
  std::cout << "........................................" << std::endl;
  
  //--------------------------------------//
  // Initialize classes relevant to workspace:
  // The combined workspace:
  m_ws = new RooWorkspace("combinedWS");
  
  // Category workspaces, and RooCategory and Simultaneous PDF:
  m_categories = new RooCategory("categories","categories");
  m_combinedPdf = new RooSimultaneous("combinedPdf","combinedPdf", *categories);
  
  // Instantiate parameter sets:
  m_nuisanceParameters = new RooArgSet();
  m_globalObservables = new RooArgSet();
  m_observables = new RooArgSet();
  m_poi = new RooArgSet();
  
  // Maps for datasets:
  m_combData.clear();
  m_combDataAsimov0.clear();
  m_combDataAsimov1.clear();
  
  // Create a map of names to expected value lists:
  m_expectedList.clear();
  
  // Create the parameter(s) of interest:
  std::vector<TString> listPOI = m_config->getStrV("MODEL_POI");
  for (int i_p = 0; i_p < (int)listPOI.size(); i_p++) {
    m_ws->factory(listPOI[i_p]);
    m_poi->add(*m_ws->var(varToName(listPOI[i_p])));
  }
  
  //--------------------------------------//
  // Loop over channels, create model for each:
  std::cout << "DHWorkspace: Looping over categories to define workspace."
	    << std::endl;
  for (m_currCateIndex = 0; m_currCateIndex < m_nCategories; m_currCateIndex++){
    m_currCateName = cateNames[m_currCateIndex];
    // Create the workspace for a single category:
    createNewCategoryWS();
  }
  
  //--------------------------------------//
  // Finished category loop, now building combined model:
  std::cout << "DHWorkspace: Beginning to combine all categories." << std::endl;
  
  // Import variable sets:  
  m_ws->defineSet("nuisanceParameters", *nuisanceParameters);
  m_ws->defineSet("globalObservables", *globalObservables);
  m_ws->defineSet("observables", *observables);
  m_ws->defineSet("poi", *m_poi);
  
  /*
    ADD SECTION CREATING PARAMETER SETS IN TH EWORKSPACE!
  */
  
  // Import the RooCategory and RooAbsPDF:
  m_ws->import(*categories);
  m_ws->import(*combinedPdf);
  
  // Define the combined dataset:
  RooRealVar wt("wt", "wt", 1);
  RooArgSet *dataArgs = new RooArgSet();
  dataArgs->add(*m_observables);
  dataArgs->add(wt);
  RooDataSet* obsData = new RooDataSet("obsData", "obsData", *dataArgs, 
				       Index(*m_categories), Import(m_combData),
				       WeightVar(wt));
  m_ws->import(*obsData);
  
  // Define the ModelConfig for the analysis and import to the workspace:
  m_modelConfig = new ModelConfig("modelConfig", m_ws);
  m_modelConfig->SetPdf((*m_ws->pdf("combinedPdf")));
  m_modelConfig->SetObservables((*m_ws->set("observables")));
  m_modelConfig->SetParametersOfInterest((*m_ws->set("poi")));
  m_modelConfig->SetNuisanceParameters((*m_ws->set("nuisanceParameters")));
  m_modelConfig->SetGlobalObservables((*m_ws->set("globalObservables")));
  m_ws->import(*m_modelConfig);
  
  std::cout << "DHWorkspace: Printing the combined workspace." << std::endl;
  m_ws->Print("v");
  
  //----------------------------------------//
  // Create Asimov data:
  
  // Set the POI to zero.
  m_ws->var(varToName(listPOI[0]))->setVal(0.0);
  
  // Do a simple background only fit before asimov data creation:
  m_ws->pdf("combinedPdf")
    ->fitTo(*m_ws->data("obsData"),
	    Minos(RooArgSet(*m_ws->set("nuisanceParameters"))),
	    SumW2Error(kTRUE));
  
  // Loop over categories for Asimov data following background-only fit:
  for (m_currCateIndex = 0; m_currCateIndex < m_nCategories; m_currCateIndex++){
    m_currCateName = cateNames[m_currCateIndex];
    
    //(*m_ws->var("nBkg_"+m_currCateName)).setVal((*m_ws->data(Form("obsData_%s",m_currCateName.Data()))).sumEntries());
  
    // Create background-only Asimov data:
    createAsimovData(0);
    // Create signal + background Asimov data:
    createAsimovData(1);
  }
  RooDataSet* asimovDataMu0 = new RooDataSet("asimovDataMu0", "asimovDataMu0",
  					     *dataArgs, Index(*categories), 
  					     Import(m_combDataAsimov0),
					     WeightVar(wt));
  RooDataSet* asimovDataMu1 = new RooDataSet("asimovDataMu1", "asimovDataMu1",
  					     *dataArgs, Index(*categories), 
  					     Import(m_combDataAsimov1),
					     WeightVar(wt));
  m_ws->import(*asimovDataMu0);
  m_ws->import(*asimovDataMu1);
  
   //----------------------------------------//
  // Save information
  
  // snapshot of original parameter values:
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*m_modelConfig->GetNuisanceParameters());
  poiAndNuis->add(*m_modelConfig->GetParametersOfInterest());
  m_ws->saveSnapshot("paramsOrigin", *poiAndNuis);
  
  // Print the workspace before saving:
  std::cout << "DHWorkspace: Printing the workspace to be saved." << std::endl;
  m_ws->Print("v");
  
  // Write workspace to file:
  m_ws->importClassCode();
  m_ws->writeToFile(Form("%s/rootfiles/workspaceDH_%s.root",
				 m_outputDir.Data(), m_anaType.Data()));
  
  //----------------------------------------//
  // Start profiling the data:
  std::cout << "DHWorkspace: Start profiling data" << std::endl;
    
  // Choose two models for the default test fits:
  std::vector<TString> sigDHModes = m_config->getStrV("sigDHModes");
  TString currDHSignal
    = m_anaType.Contains("NonRes") ? sigDHModes[0] : sigDHModes[1];
  
  DHTestStat *dhts 
    = new DHTestStat(m_configFile, currDHSignal, "FromFile", m_ws);
  dhts->saveSnapshots(true);
  dhts->setPlotDirectory(Form("%s/Plots/", m_outputDir.Data()));
  dhts->setPlotAxis(true, 0.005, 50);
  
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
  std::cout << "DHWorkspace: Printing likelihood results for " << m_anaType
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
  
}

/**
   -----------------------------------------------------------------------------
   Create the workspace for a single analysis category. Note: m_anaType, 
   m_currCateIndex, and m_currCateName are defined with global scope and updated
   in each call to this method.
*/
void DHWorkspace::createNewCategoryWS() {
  
  // The bools that control the systematic uncertainties:
  bool channel_constraints_attached = (m_currCateIndex == m_nCategories-1);

  // Add a new category:
  m_categories->defineType(m_currCateName);
  
  //--------------------------------------//
  // Add systematic uncertainties first:
  // (builds constraint term, global observables, nuisance parameters, expected)
  if (m_currCateIndex == 0) {
    m_constraints = new RooArgSet();
    std::vector<TString> systematics = m_config->getStrV("SysSources");
    for (int i_s = 0; i_s < (int)systematics.size(); i_s++) {
      addSystematic(m_config->getStr(Form("SysForm_%s",
					  (systematics[i_s]).Data())));
    }
    // Then create RooProduct of constraints:
    RooProduct constraint("constraint", "constraint", *m_constraints);
    m_ws->import(constraint);
    
    // Iterate over the expected sets and create products. 
    for (std::map<TString,RooArgSet*>::iterator iterEx = m_expectedList.begin();
	 iterEx != m_expectedList.end(); iterEx++) {
      TString currExpName = iterEx->first;
      RooArgSet *currExp = iterEx->second;
      RooProduct currExpectation(currExpName, currExpName, *currExp);
      m_ws->import(currExpectation);
    }
  }
  
  //--------------------------------------//
  // Add the observable:
  TString observable = m_config->getStr(Form("OBS_%s", m_currCateName.Data()));
  TString observableName = varToName(observable);
  m_ws->factory(observable);
  m_observables->add(m_ws->var(observableName));
  
  //--------------------------------------//
  // Add variables unrelated to systematics:
  std::vector<TString> variables
    = m_config->getStrV(Form("VARS_%s", m_currCateName.Data()));
  for (int i_v = 0; i_v < (int)variables.size(); i_v++) {
    TString variableName = varToName(variables[i_v]);
    // Make sure variable hasn't already been created:
    if (!m_ws->var(variableName)) {
      // Create in workspace:
      m_ws->factory(variables[i_v]);
      // Also add to list of nuisance parameters:
      m_nuisanceParameters->add(*m_ws->var(variableName));
    }
  }
  
  //--------------------------------------//
  // Add expressions:
  std::vector<TString> expressions
    = m_config->getStrV(Form("EXPRS_%s", m_currCateName.Data()));
  for (int i_e = 0; i_e < (int)expressions.size(); i_e++) {
    TString expressionName = funcToName(expressions[i_e]);
    // Create expression in workspace if not done already:
    if (!m_ws->function(expressionName)) m_ws->factory(expressions[i_e]);
  }
  
  //--------------------------------------//
  // Add component PDFs:
  std::vector<TString> pdfs
    = m_config->getStrV(Form("PDFS_%s", m_currCateName.Data()));
  for (int i_p = 0; i_p < (int)pdfs.size(); i_p++) {
    TString pdfName = funcToName(pdfs[i_p]);
    // Create PDF in workspace if it doesn't already exist:
    if (!m_ws->pdf(pdfName)) m_ws->factory(pdfs[i_p]);
  }
  
  //--------------------------------------//
  // Build the complete model:
  TString modelForm = m_config->getStr(Form("MODEL_%s",m_currCateName.Data()));
  m_ws->factory(modelForm);
  
  // Add model to combined PDF:
  TString modelName = funcToName(modelForm);
  m_combinedPdf->addPdf(*cateWS[m_currCateIndex]->pdf(namePdf), m_currCateName);
  
  // Attach constraint term to the first category:
  if (m_currCateIndex == 0) {
    tempWS->factory(Form("PROD::model_%s(%s,constraint)",
			 m_currCateName.Data(), modelName.Data()));
  }
  
  //--------------------------------------//
  // Import the observed data set:
  DHDataReader *dhdr = new DHDataReader(m_configFile,m_ws->var(observableName));
  RooDataSet *obsData = dhdr->loadData(m_currCateName);
  
  // Rename the imported dataset for consistency:
  TString dataName = Form("obsData_%s", m_currCateName.Data());
  obsData->SetNameTitle(dataName, dataName);
  m_ws->import(*obsData);
  m_combData[(string)m_currCateName] = (RooDataSet*)m_ws->data(dataName);
  
  std::cout << "DHWorkspace: Finished import for category " << m_currCateName 
	    << std::endl;
}

/**
   -----------------------------------------------------------------------------
   systematicForm will have:
   name[constr=gaus,center=1.0,,type=yield,comp=name~val,...,compLo=name~val...]
   
   constraint = "gaus", "logn", "asym" "bifur"
   central value = (1, 0)
   magnitude = percent, use sign to indicate corr/anti-corr
   magnitude = percent, use sign to indicate corr/anti-corr, for asym only
   type = "yield", "res", "scale"
   componentsLo = components to which this applies (can have many).
   components = components to which this applies (can have many).
   So that migration systematics are really just yield systematics that are 
   correlated between categories. 
   
   @param systematicForm - The form of the systematic, as outlined above.
*/
void DHWorkspace::addSystematic(TString systematicForm) {
  // Get the systematic name:
  TString systematicName = varToName(systematicForm);
  systematicForm.ReplaceAll(name,"");
  systematicForm.ReplaceAll("[","");
  systematicForm.ReplaceAll("]","");
  
  //----------------------------------------//
  // Retrieve systematic configuration data
  
  // Systematic settings that will be retrieved from the name:
  TString centralValue = "";
  TString constraint = "";
  TString type = "";
  std::map<TString,double> componentsLo; componentsLo.clear();
  std::map<TString,double> components; components.clear();
  
  // First split up the item into components:
  TObjArray *array = systematicForm.Tokenize(",");
  for (int i_t = 0; i_t < array->GetEntries(); i_t++) {
    TString currElement = ((TObjString*)array->At(i_e))->GetString();
    if (currElement.Contains("constr=")) {
      constraint = currElement;
      constraint.ReplaceAll("constr=","");
    }
    else if (currElement.Contains("center=")) {
      centralValue = currElement;
      centralValue.ReplaceAll("center=","");
    }
    else if (currElement.Contains("type=")) {
      type = currElement;
      type.ReplaceAll("type=","");
    }
    else if (currElement.Contains("compLo=")) {
      TString currComponentInfo = currElement;
      currComponentInfo.ReplaceAll("compLo=","");
      TString currComponentName = currComponentInfo;
      currComponentName.Remove(First("~"));
      TString currComponentVal = currComponentInfo;
      currComponentVal.Remove(0,name.First("~")+1);
      componentsLo[currComponentName] = currComponentVal.Atof();
    }
    else if (currElement.Contains("comp=")) {
      TString currComponentInfo = currElement;
      currComponentInfo.ReplaceAll("comp=","");
      TString currComponentName = currComponentInfo;
      currComponentName.Remove(First("~"));
      TString currComponentVal = currComponentInfo;
      currComponentVal.Remove(0,name.First("~")+1);
      components[currComponentName] = currComponentVal.Atof();
    }
  }
  
  // Check that all inputs have been provided:
  if ((centralValue.EqualTo("") || constraint.EqualTo("") || 
       type.EqualTo("") || (int)components.size() == 0) ||
      (constraint.EqualTo("asym") && (int)componentsLo.size() == 0)) {
    std::cout << "DHWOrkspace: ERROR loading systematic " << systematicForm 
	      << std::endl;
    exit(0);
  }
  
  //----------------------------------------//
  // Create the nuisance parameter, global observable, and constraint PDF:
  if (!m_ws->var(systematicName)) {
    // Bifurcated Gaussian constraint requires magnitude information:
    if (constraint.EqualTo("bifur")) {
      m_ws->factory(Form("%s[0,-5,5]", systematicName.Data()));
      m_ws->factory(Form("RNDM_%s[0,-5,5]", systematicName.Data()));
      std::map<TString,double>::iterator iterComponent = components.begin();
      TString componentName = iterComponent->firt;
      double magnitude = iterComponent->second;
      double magnitudeLo = (componentsLo.count(componentName) > 0) ?
	componentsLo[componentName] : 0.0;
      TString magnitudeRatio = (magnitude / magnitudeLo);
      m_ws->factory(Form("RooBifurGauss::constr_%s(RNDM_%s,%s,1,%f)",
			 systematicName.Data(), systematicName.Data(),
			 systematicName.Data(), magnitudeRatio));
    }
    // Simple Gaussian constraint is straightforward:
    else m_ws->factory(Form("RooGaussian::constr_%s(RNDM_%s,%s,1)"));
    
    // For now, only allow gaussian constraint term:
    m_ws->factory(Form("RooGaussian::constr_%s(RNDM_%s,%s,1)"));
    
    // Add the new objects to the relevant sets:
    m_nuisanceParameters->add(*m_ws_>var(systematicName));
    m_globalObservables->add(*m_ws->var(Form("RNDM_%s",systematicName.Data())));
    m_constraints->add(*m_ws->pdf(Form("constr_%s",systematicName.Data())));
  }
  
  //----------------------------------------//
  // Loop over components for this systematic and create expectation terms:
  for (std::map<TString,double>::iterator iterComponent = components.begin();
       iterComponent != components.end(); iterComponent++) {
    
    // Get the component name and magnitude:
    TString componentName = iterComponent->firt;
    double magnitude = iterComponent->second;
    double magnitudeLo = (componentsLo.count(componentName) > 0) ?
      componentsLo[componentName] : 0.0;
    
    // Define the variable name and expected name:
    TString varName = Form("%s_%s",systematicName.Data(),componentName.Data());
    TString expName = Form("expected_%s", varName.Data());
    
    // Asymmetric uncertainty:
    if (constraint.EqualTo("asym")) {
      std::cout << "DHWorkspace: " << systematicName << " has asym. uncertainty"
		<< std::endl;
      m_ws->factory(Form("prod::%_times_beta(%s,beta_%s[1]",
			 systematicName.Data(), varName.Data()));
      vector<double> magHi; magHi.clear(); magHi.push_back(1+magnitude);
      vector<double> magLo; magLo.clear(); magLo.push_back(1-magnitudeLo);
      RooArgList nuisList(*m_ws->factory(Form("%_times_beta",
					      systematicName.Data())));
      RooStats::HistFactory::FlexibleInterpVar
	expVar(expName, expName, nuisList, centralValue.Atof(), magLo, magHi);
      workspace->import(expVar);
    }
    
    // Gaussian constraint:
    else if (constraint.EqualTo("gaus")) {
      std::cout << "DHWorkspace: " << systematicName << " has gaus uncertainty"
		<< std::endl;
      m_ws->factory(Form("sum::%s(nominal_%s[%f], prod::uncer_%s(prod::%s_times_beta(%s, beta_%s[1]), sigma_%s[%f]))", expName.Data(), varName.Data(), centralValue.Atof(), varName.Data(), varName.Data(), systematicName.Data(), varName.Data(), varName.Data(), magnitude));
    }
    
    // Lognormal constraint:
    else if (constraint.EqualTo("logn")) {
      std::cout << "DHWorkspace: " << systematicName << " has gaus uncertainty"
		<< std::endl;
      double valLogKappa = sqrt(log(1+pow(magnitude,2)));
      m_ws->factory(Form("logKappa_%s[%f]", varName.Data(), valLogKappa));
      m_ws->factory(Form("RooExponential::expTerm_%s(prod::%s_times_beta(%s, beta_%s[1]), logKappa_%s[%f])", varName.Data(), varName.Data(), systematicName.Data(), varName.Data(), varName.Data(), valLogKappa));
      m_ws->factory(Form("prod::%s(expTerm_%s,nominal_%s[%f])", expName.Data(), varName.Data(), varName.Data(), centralValue.Atof()));
    }
    
    // Exit if improper constraint type:
    else {
      std::cout << "DHWorkspace: ERROR! Unknown constraint type " << constraint
		<< std::endl;
      exit(0);
    }
    
    // The key should depend on the systematic type and the component name:
    TString expKey = Form("expectation_%s_%s",type.Data(),componentName.Data());
    // Check if expected quantity exists, instantiate if not:
    if (m_expectedList.count(expKey) == 0) {
      m_expectedList[expKey] = new RooArgSet();
    }
    m_expectedList[expKey]->add(*m_ws->function(expName));
  }
}

/*
   -----------------------------------------------------------------------------
   Plot the data and fit in a single category
*/
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

/**
   -----------------------------------------------------------------------------
   Create Asimov data for the statistical model, using a fit to observed data
   for the shape and normalizaiton of the background.
   @param valMuDH - the value of the di-Higgs signal strength to use.
*/
void DHWorkspace::createAsimovData(int valMuDH) {
  std::cout << "DHWorkspace: Creating Asimov data, mu=" << valMuDH << std::endl;
  
  // Set mu_DM and mu_SM to the specified values:
  RooRealVar *poi = NULL;
  double initialMuDH = 0.0;
  poi = m_ws->var("mu_DH");
  initialMuDH = poi->getVal();
  poi->setVal(valMuDH);
  poi->setConstant(true);  
    
  TString asimovName = Form("asimovDataMu%d_%s",valMuDH,m_currCateName.Data());
  RooDataSet *currAsimov = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*m_ws->pdf(Form("model_%s",m_currCateName.Data())), *m_ws->set("observables"));
  currAsimov->SetNameTitle(asimovName, asimovName);
  m_ws->import(*currAsimov);
  
  m_ws->var("mu_DH")->setVal(initialMuDH);
  m_ws->var("mu_DH")->setConstant(false);
    
  std::cout << "DHWorkspace: Asimov data has " << currAsimov->sumEntries() 
	    << " entries" << std::endl;
  // Import the dataset to the workspace and add to data map:
  m_ws->import(*currAsimov);
  if (valMuDH == 1) {
    m_combDataAsimov1[(string)m_currCateName]
      = (RooDataSet*)m_ws->data(asimovName);
  }
  else if (valMuDH == 0) {
    m_combDataAsimov0[(string)m_currCateName]
      = (RooDataSet*)m_ws->data(asimovName);
  }
}

/**
   -----------------------------------------------------------------------------
   Takes in a variable declaration such as "name[1,0,5]" and returns "name".
   @param varForm - The form of the variable declared.
   @return - The name of the variable without the rest of the expression.
*/
TString DHWorkspace::varToName(TString varForm) {
  TString name = varForm;
  name.Remove(name.First("["));
  return name;
}

/**
   -----------------------------------------------------------------------------
   Takes in a function or pdf declaration such as "type::pdfname()" and returns
   "pdfname".
   @param funcForm - The form of the function declared.
   @return - The name of the function without the rest of the expression.
*/
TString DHWorkspace::funcToName(TString funcForm) {
  TString name = funcForm;
  name.Remove(name.First("("));
  name.Remove(0,name.First(":")+2);
  return name;
}
