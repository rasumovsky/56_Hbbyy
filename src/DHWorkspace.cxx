////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHWorkspace.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 01/14/2016                                                          //
//                                                                            //
//  This class builds the workspace for the resonant and non-resonant di-     //
//  Higgs search at 13 TeV. Most of the model configuration is done in the    //
//  config file, so that the code does not have to be recompiled each time a  //
//  change needs to be made to the fit model.                                 //
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
   @param newOptions - The job options ("New", "FromFile"), etc.
   @return void
*/
DHWorkspace::DHWorkspace(TString newConfigFile, TString newOptions) {
  m_configFile = newConfigFile;
  m_options = newOptions;
  m_allGoodFits = true;
  
  m_ws = NULL;
  m_modelConfig = NULL;
    
  m_config = new Config(m_configFile);
  m_dataToPlot = (m_config->getBool("DoBlind")) ? "asimovDataMu1" : "obsData";
  m_anaType = m_config->getStr("AnalysisType");
  
  // Print config file:
  if (m_config->getBool("Verbose")) {
    std::cout << "DHWorkspace: Printing config file" << std::endl;
    m_config->printDB();
  }
  
  // Print workspace inputs:
  std::cout << "\nDHWorkspace: Initializing..."
	    << "\n\tconfigFile = " << m_configFile 
	    << "\n\toptions = " << m_options << std::endl;
  
  // Assign output directory, and make sure it exists:
  m_outputDir = Form("%s/%s/DHWorkspace", 
		     (m_config->getStr("MasterOutput")).Data(),
		     (m_config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  system(Form("mkdir -vp %s/Plots/", m_outputDir.Data()));
  system(Form("mkdir -vp %s/rootfiles/", m_outputDir.Data()));
  system(Form("mkdir -vp %s/PoI/", m_outputDir.Data()));
  
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
   Create the workspace for a single analysis category. Note: m_currCateIndex,
   and m_currCateName are defined with global scope and updated in each call to 
   this method.
*/
void DHWorkspace::addCategory() {
  printer("DHWorkspace::addCategory()", false);
    
  // Add a new type to the RooCategory:
  m_categories->defineType(m_currCateName);
  
  //--------------------------------------//
  // Add systematic uncertainties first:
  // (constraint term, global observables, nuisance parameters, expectation)
  if (m_currCateIndex == 0) {
    m_constraints = new RooArgSet();
    if (m_config->isDefined("SysSources")) {
      std::vector<TString> systematics = m_config->getStrV("SysSources");
      for (int i_s = 0; i_s < (int)systematics.size(); i_s++) {
	addSystematic(m_config->getStr(Form("SysForm_%s",
					    (systematics[i_s]).Data())));
      }
    }
    
    // Create RooProduct of constraints:
    RooProdPdf constraint("constraint", "constraint", *m_constraints);
    m_ws->import(constraint);
        
    // Iterate over the expected sets and create products:
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
  TString observableName = nameOfVar(observable);
  m_ws->factory(observable);
  m_observables->add(*m_ws->var(observableName));
  
  //--------------------------------------//
  // Add variables unrelated to systematics:
  printer("DHWorkspace::addCategory(): Adding non-syst. variables...", false);
  
  if (m_config->isDefined(Form("VARS_%s", m_currCateName.Data()))) {
    std::vector<TString> variables
      = m_config->getStrV(Form("VARS_%s", m_currCateName.Data()));
    for (int i_v = 0; i_v < (int)variables.size(); i_v++) {
      TString variableName = nameOfVar(variables[i_v]);
      // Make sure variable hasn't already been created:
      if (!m_ws->var(variableName)) {
	// Create in workspace:
	m_ws->factory(variables[i_v]);
	// Also add to list of nuisance parameters:
	m_nuisanceParameters->add(*m_ws->var(variableName));
      }
    }
  }
  
  //--------------------------------------//
  // Add expressions:
  printer("DHWorkspace::addCategory(): Adding expressions...", false);
  
  if (m_config->isDefined(Form("EXPRS_%s", m_currCateName.Data()))) {
    std::vector<TString> expressions
      = m_config->getStrV(Form("EXPRS_%s", m_currCateName.Data()));
    for (int i_e = 0; i_e < (int)expressions.size(); i_e++) {
      TString expressionName = nameOfFunc(expressions[i_e]);
      // Create expression in workspace if not done already:
      if (!m_ws->function(expressionName)) m_ws->factory(expressions[i_e]);
    }
  }
  
  //--------------------------------------//
  // Add component PDFs:
  printer("DHWorkspace::addCategory(): Adding component PDFs...", false);
  
  if (m_config->isDefined(Form("PDFS_%s", m_currCateName.Data()))) {
    std::vector<TString> pdfs
      = m_config->getStrV(Form("PDFS_%s", m_currCateName.Data()));
    for (int i_p = 0; i_p < (int)pdfs.size(); i_p++) {
      TString pdfName = nameOfFunc(pdfs[i_p]);
      // Create PDF in workspace if it doesn't already exist:
      if (!m_ws->pdf(pdfName)) m_ws->factory(pdfs[i_p]);
    }
  }
  
  //--------------------------------------//
  // Build the complete model:
  printer("DHWorkspace::addCategory(): Building complete model...", false);
  
  TString modelForm
    = m_config->getStr(Form("MODEL_%s", m_currCateName.Data()), false);
  TString modelName = nameOfFunc(modelForm);
  m_ws->factory(modelForm);
  // Attach constraint term to the first category:
  if (m_currCateIndex == 0) {
    m_ws->factory(Form("PROD::model_%s(%s,constraint)",
		       m_currCateName.Data(), modelName.Data()));
  }
  else {
    m_ws->factory(Form("PROD::model_%s(%s)", 
		       m_currCateName.Data(), modelName.Data()));
  }
  // Add model including constraint to combined PDF:
  m_combinedPdf->addPdf(*m_ws->pdf(Form("model_%s", m_currCateName.Data())), 
			m_currCateName);
  
  //--------------------------------------//
  // Import the observed data set:
  printer("DHWorkspace::addCategory(): Importing observed dataset...", false);
    
  // Create RooDataSet object:
  TString dataName = Form("obsData_%s", m_currCateName.Data());
  RooRealVar wt("wt", "wt", 1.0);
  RooDataSet *obsData = new RooDataSet(dataName, dataName,
				       RooArgSet(*m_ws->var(observableName),wt),
				       RooFit::WeightVar(wt));
  // Open input text file to read mass points (and possibly weights...):
  TString textFileName = Form("%s/%s/DHMassPoints/masspoint_%s_%s.txt", 
			      (m_config->getStr("MasterOutput")).Data(),
			      (m_config->getStr("JobName")).Data(),
			      m_anaType.Data(), m_currCateName.Data());
  std::ifstream massInput(textFileName);
  if (massInput.is_open()) {
    std::cout << "DHWorkspace: Loading data points from file " << textFileName
	      << std::endl;
    double currMass; double currWeight;
    while (!massInput.eof()) {
      massInput >> currMass >> currWeight;
      m_ws->var(observableName)->setVal(currMass);
      wt.setVal(currWeight);
      obsData->add(RooArgSet(*m_ws->var(observableName),wt), currWeight);
    }
  }
  else printer(Form("DHWorkspace: ERROR! nonexistent input %s", 
		    textFileName.Data()), true);
  massInput.close();
  
  // Import RooDataSet into workspace, and add to dataset map (for combination):
  m_ws->import(*obsData);
  m_combData[(string)m_currCateName] = (RooDataSet*)m_ws->data(dataName);
  
  // Then make a quick plot before the end:
  plotCatePdfAndData(50);
  
  printer(Form("DHWorkspace: Added cate. %s", m_currCateName.Data()), false);
}

/**
   -----------------------------------------------------------------------------
   Add a systematic uncertainty to the analysis workspace. This is only done 
   once (for the first category, to which the constraint term is attached). 
   This method creates the nuisance parameter, global observable, expectation
   function, and constraint term. 
   
   systematicForm will have the following form (no spaces):
   name[constr=gaus,center=1.0,type=yield,comp=name~val,...,compLo=name~val...]
   
   constr = constraint type ("gaus", "logn", "asym" "bifur")
   center = central value (usually 1.0 or 0.0)
   type = the type of systematic (e.g. "yield", "res", "scale")
   comp = a component and its sys value, separated by ~
          Usually value in percent, use sign to indicate corr/anti-corr
   compLo = same as comp, but for asymmetric uncertainties. 
   Both comp and compLo tolerate multiple entries (they form a vector).
      
   @param systematicForm - The form of the systematic, as outlined above.
*/
void DHWorkspace::addSystematic(TString systematicForm) {
  printer(Form("DHWorkspace::addSystematic(%s)",systematicForm.Data()), false);
  
  // Get the systematic name:
  TString systematicName = nameOfVar(systematicForm);
  systematicForm.ReplaceAll(systematicName,"");
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
  
  // First split up the systematic form into components:
  TObjArray *array = systematicForm.Tokenize(",");
  for (int i_t = 0; i_t < array->GetEntries(); i_t++) {
    TString currElement = ((TObjString*)array->At(i_t))->GetString();
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
      currComponentName.Remove(currComponentName.First("~"));
      TString currComponentVal = currComponentInfo;
      currComponentVal.Remove(0,currComponentVal.First("~")+1);
      componentsLo[currComponentName] = currComponentVal.Atof();
    }
    else if (currElement.Contains("comp=")) {
      TString currComponentInfo = currElement;
      currComponentInfo.ReplaceAll("comp=","");
      TString currComponentName = currComponentInfo;
      currComponentName.Remove(currComponentName.First("~"));
      TString currComponentVal = currComponentInfo;
      currComponentVal.Remove(0,currComponentVal.First("~")+1);
      components[currComponentName] = currComponentVal.Atof();
      if (m_config->getBool("Verbose")) {
	std::cout << "\tcomponent name = " << currComponentName << std::endl;
      }
    }
  }
  
  // Check that all inputs have been provided:
  if ((centralValue.EqualTo("") || constraint.EqualTo("") || 
       type.EqualTo("") || (int)components.size() == 0) ||
      (constraint.EqualTo("asym") && (int)componentsLo.size() == 0)) {
    printer(Form("DHWorkspace: ERROR loading syst. %s", systematicForm.Data()),
	    true);
  }
  
  //----------------------------------------//
  // Create the nuisance parameter, global observable, and constraint PDF:
  if (!m_ws->var(systematicName)) {
    
    // Nuisance parameter and global observable:
    m_ws->factory(Form("%s[0,-5,5]", systematicName.Data()));
    m_ws->factory(Form("RNDM_%s[0,-5,5]", systematicName.Data()));
    
    // Bifurcated Gaussian constraint requires magnitude information:
    if (constraint.EqualTo("bifur")) {
      std::map<TString,double>::iterator iterComponent = components.begin();
      TString componentName = iterComponent->first;
      double magnitude = iterComponent->second;
      double magnitudeLo = (componentsLo.count(componentName) > 0) ?
	componentsLo[componentName] : 0.0;
      double magnitudeRatio = (magnitude / magnitudeLo);
      m_ws->factory(Form("RooBifurGauss::constr_%s(RNDM_%s,%s,1,%f)",
			 systematicName.Data(), systematicName.Data(),
			 systematicName.Data(), magnitudeRatio));
    }
   
    // Simple Gaussian constraint is straightforward:
    else m_ws->factory(Form("RooGaussian::constr_%s(RNDM_%s,%s,1)", 
			    systematicName.Data(), systematicName.Data(), 
			    systematicName.Data()));
    
    // Add the new objects to the relevant sets:
    m_nuisanceParameters->add(*m_ws->var(systematicName));
    m_globalObservables->add(*m_ws->var(Form("RNDM_%s",systematicName.Data())));
    m_constraints->add(*m_ws->pdf(Form("constr_%s",systematicName.Data())));
  }
  
  //----------------------------------------//
  // Loop over components for this systematic and create expectation terms:
  for (std::map<TString,double>::iterator iterComponent = components.begin();
       iterComponent != components.end(); iterComponent++) {
    
    // Get the component name and magnitude:
    TString componentName = iterComponent->first;
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
      m_ws->factory(Form("prod::%s_times_beta(%s,beta_%s[1]",
			 varName.Data(),systematicName.Data(),varName.Data()));
      vector<double> magHi; magHi.clear(); magHi.push_back(1+magnitude);
      vector<double> magLo; magLo.clear(); magLo.push_back(1-magnitudeLo);
      RooArgList nuisList(*m_ws->factory(Form("%s_times_beta",
					      systematicName.Data())));
      RooStats::HistFactory::FlexibleInterpVar
	expVar(expName, expName, nuisList, centralValue.Atof(), magLo, magHi);
      m_ws->import(expVar);
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
      m_ws->factory(Form("RooExponential::expTerm_%s(prod::%s_times_beta(%s,beta_%s[1]),logKappa_%s[%f])", varName.Data(), varName.Data(), systematicName.Data(), varName.Data(), varName.Data(), valLogKappa));
      m_ws->factory(Form("prod::%s(expTerm_%s,nominal_%s[%f])", expName.Data(), varName.Data(), varName.Data(), centralValue.Atof()));
    }
    
    // Exit if improper constraint type:
    else {
      printer(Form("DHWorkspace: Constraint type ERROR %s", constraint.Data()),
	      true);
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

/**
   -----------------------------------------------------------------------------
   Create Asimov data for the statistical model, using a fit to observed data
   for the shape and normalizaiton of the background.
   @param valPoI - The value of the parameter of interest.
*/
void DHWorkspace::createAsimovData(double valPoI) {
  printer(Form("DHWorkspace::createAsimovData(%2.2f)",valPoI), false);
  
  // Set the parameter of interest to specified value:
  //RooRealVar *poi = m_ws->var("npbBSM");
  //m_poi = new RooArgSet();
  //m_ws->defineSet("poi", *m_poi);
  

  RooRealVar *poi = (RooRealVar*)((RooArgSet*)m_ws->set("poi"))->first();





  double initialPoIVal = poi->getVal();
  poi->setVal(valPoI);
  poi->setConstant(true);  
  
  // Create the Asimov dataset using RooStats utility:
  TString asimovName = (valPoI > 0) ? 
    Form("asimovDataMu1_%s",m_currCateName.Data()) :
    Form("asimovDataMu0_%s",m_currCateName.Data());
  RooDataSet *currAsimov = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*m_ws->pdf(Form("model_%s",m_currCateName.Data())), *m_ws->set("observables"));
  currAsimov->SetNameTitle(asimovName, asimovName);
  m_ws->import(*currAsimov);
  // Import the dataset to the workspace and add to data map:
  m_ws->import(*currAsimov);
  if (valPoI > 0) {
    m_combDataAsimov1[(string)m_currCateName]
      = (RooDataSet*)m_ws->data(asimovName);
  }
  else {
    m_combDataAsimov0[(string)m_currCateName]
      = (RooDataSet*)m_ws->data(asimovName);
  }
  printer(Form("DHWorkspace: Asimov data entries=%f",currAsimov->sumEntries()),
	  false); 
  
  // Return PoI to original settings:
  poi->setVal(initialPoIVal);
  poi->setConstant(false);
}

/**
   -----------------------------------------------------------------------------
   Create a multi-category workspace from scratch. This workspace will include
   the PDFs as well as the dataset to fit and Asimov data for S+B and B-only
   hypotheses.
*/
void DHWorkspace::createNewWS() {
  printer("DHWorkspace::createNewWS()", false);
  
  if (m_config->getBool("Verbose")) {
    std::cout << "\n........................................" << std::endl;
    std::cout << "Luminosity at 13 TeV: "
	      << m_config->getStr("AnalysisLuminosity") << " pb-1" << std::endl;
  }
  
  // Define and name analysis categories:
  std::vector<TString> cateNames = m_config->getStrV("CateNames");
  m_nCategories = (int)cateNames.size();
  if (m_config->getBool("Verbose")) {
    std::cout << "  Number of categories = " << m_nCategories << std::endl;
    std::cout << "........................................" << std::endl;
  }
  
  //--------------------------------------//
  // Initialize objects necessary for workspace creation:
  // The combined workspace:
  m_ws = new RooWorkspace("combinedWS");
  
  // RooCategory and Simultaneous PDF:
  m_categories = new RooCategory("categories", "categories");
  m_combinedPdf = new RooSimultaneous("combinedPdf", "combinedPdf",
				      *m_categories);

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
  std::vector<TString> listPOI = m_config->getStrV("Model_PoI");
  for (int i_p = 0; i_p < (int)listPOI.size(); i_p++) {
    m_ws->factory(listPOI[i_p]);
    m_poi->add(*m_ws->var(nameOfVar(listPOI[i_p])));
  }

  //--------------------------------------//
  // Loop over channels, create model for each:
  std::cout << "DHWorkspace: Looping over categories to define workspace."
	    << std::endl;
  for (m_currCateIndex = 0; m_currCateIndex < m_nCategories; m_currCateIndex++){
    m_currCateName = cateNames[m_currCateIndex];
    // Create the workspace for a single category:
    addCategory(); // m_currCateName and m_currCateIndex globally scoped...
  }
  
  //--------------------------------------//
  // Finished category loop, now building combined model:
  std::cout << "DHWorkspace: Beginning to combine all categories." << std::endl;
  
  // Import variable sets:  
  m_ws->defineSet("nuisanceParameters", *m_nuisanceParameters);
  m_ws->defineSet("globalObservables", *m_globalObservables);
  m_ws->defineSet("observables", *m_observables);
  m_ws->defineSet("poi", *m_poi);
  
  // Import the RooCategory and RooAbsPDF:
  m_ws->import(*m_categories);
  m_ws->import(*m_combinedPdf);
  
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
  m_ws->var(nameOfVar(listPOI[0]))->setVal(0.0);
  
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
  					     *dataArgs, Index(*m_categories), 
  					     Import(m_combDataAsimov0),
					     WeightVar(wt));
  RooDataSet* asimovDataMu1 = new RooDataSet("asimovDataMu1", "asimovDataMu1",
  					     *dataArgs, Index(*m_categories), 
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
  
  //std::cout << "STOPPING SHORT FOR DEV." << std::endl;
  //exit(0);
  
  //----------------------------------------//
  // Start profiling the data (must do after writing workspace for some reason):
  printer("DHWorkspace: Start profiling data", false);
  
  DHTestStat *dhts = new DHTestStat(m_configFile, "FromFile", m_ws);
  dhts->saveSnapshots(true);
  dhts->setPlotDirectory(Form("%s/Plots/", m_outputDir.Data()));
  dhts->setPlotAxis(true, 0.005, 50);
  
  double profiledPOIVal = -999.0;
  // Mu = 0 fits:
  double nllMu0 = dhts->getFitNLL(m_dataToPlot, 0, true, profiledPOIVal);
  if (!dhts->fitsAllConverged()) m_allGoodFits = false;
  // Mu = 1 fits:
  double nllMu1 = dhts->getFitNLL(m_dataToPlot, 1, true, profiledPOIVal);
  if (!dhts->fitsAllConverged()) m_allGoodFits = false;
  // Mu free fits:
  double nllMuFree = dhts->getFitNLL(m_dataToPlot, 1, false, profiledPOIVal);
  if (!dhts->fitsAllConverged()) m_allGoodFits = false;
  
  // Print summary of the fits:
  std::cout.precision(10);
  std::cout << "DHWorkspace: Printing likelihood results." << std::endl;
  std::cout << "\tnll(muDH = 1):  " << nllMu1 << std::endl;
  std::cout << "\tnll(muDH = 0):  " << nllMu0 << std::endl;
  std::cout << "\tnll(muDH free): " << nllMuFree << std::endl;
  std::cout << " " << endl;
  std::cout << "\tnll(S+B)/nll(B) " << nllMu1 - nllMu0 << std::endl;
  std::cout << "\tnll(muDH=1)/nll(muhat) = " << nllMu1-nllMuFree << std::endl;
  std::cout << "\tnll(muDH=0)/nll(muhat) = " << nllMu0-nllMuFree << std::endl;
  if (m_allGoodFits) std::cout << "all good fits = TRUE" << std::endl;
  else std::cout << "all good fits = FALSE" << std::endl;
  std::cout << "Profiled muDH value : " << profiledPOIVal << std::endl;
  
  // Write the profiled mu value to file:
  ofstream filePoI;
  filePoI.open(Form("%s/PoI/poi_%s.txt", m_outputDir.Data(), m_anaType.Data()));
  filePoI << profiledPOIVal << std::endl;
  filePoI.close();
  
  delete dhts;
}

/**
   -----------------------------------------------------------------------------
   Checks whether all of the fits converged.
   @return - True iff the fits all converged.
*/
bool DHWorkspace::fitsAllConverged() {
  printer("DHWorkspace::fitsAllConverge()", false);
  return m_allGoodFits;
}

/**
   -----------------------------------------------------------------------------
   Retrieves the workspace created by this program.
   @return - The workspace object storing all class data.
*/
RooWorkspace* DHWorkspace::getCombinedWorkspace() {
  printer("DHWorkspace::getCombinedWorkspace()", false);
  return m_ws;
}

/**
   -----------------------------------------------------------------------------
   Retrieves a pointer to the model config.
   @return - The model configuration storing class data.
*/
ModelConfig* DHWorkspace::getModelConfig() {
  printer("DHWorkspace::getModelConfig()", false);
  return m_modelConfig;
}

/**
   -----------------------------------------------------------------------------
   Load an old workspace (basically just setting model configuration and 
   workspace pointers).
*/
void DHWorkspace::loadWSFromFile() {
  printer("DHWorkspace::loadWSFromFile()", false);
  
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
   Takes in a variable declaration such as "name[1,0,5]" and returns "name".
   @param varForm - The form of the variable declared.
   @return - The name of the variable without the rest of the expression.
*/
TString DHWorkspace::nameOfVar(TString varForm) {
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
TString DHWorkspace::nameOfFunc(TString funcForm) {
  TString name = funcForm;
  name.Remove(name.First("("));
  name.Remove(0,name.First(":")+2);
  return name;
}

/**
   -----------------------------------------------------------------------------
   Plot the data and fit in a single category
   @param nBins - The number of bins in the RooPlot.
*/
void DHWorkspace::plotCatePdfAndData(int nBins) {
  printer(Form("DHWorkspace::plotCatePdfAndData(%d)",nBins), false);
    
  // Observable information:
  TString observable = m_config->getStr(Form("OBS_%s", m_currCateName.Data()));
  TString observableName = nameOfVar(observable);
  
  // Create canvas and RooPlot:
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  RooPlot* frame = (*m_ws->var(observableName)).frame(nBins);
  m_ws->data(Form("obsData_%s", m_currCateName.Data()))->plotOn(frame);
  TString modelName = Form("model_%s", m_currCateName.Data());
  (*m_ws->pdf(modelName)).plotOn(frame, LineColor(2));
  frame->SetYTitle("Entries");
  frame->SetXTitle("Mass [GeV]");
  frame->Draw();
  
  // Print the category name:
  TLatex text; text.SetNDC(); text.SetTextColor(1);
  text.DrawLatex(0.2, 0.81, m_currCateName);
  
  // Print then delete the canvas:
  can->Print(Form("%s/Plots/cateFit_%s.eps",
		  m_outputDir.Data(), m_currCateName.Data()));
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Prints a statement (if verbose) and exits (if fatal).
   @param statement - The statement to print.
   @param isFatal - True iff. this should trigger an exit command.
*/
void DHWorkspace::printer(TString statement, bool isFatal) {
  if (m_config->getBool("Verbose") || isFatal) {
    std::cout << statement << std::endl;
  }
  if (isFatal) exit(0);
}
