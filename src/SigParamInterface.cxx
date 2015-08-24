////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParamInterface.cxx                                               //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/06/2015                                                          //
//                                                                            //
//  This is an interface class for the SigParam class. It allows one to load  //
//  preexisting workspaces or create new ones if they cannot be loaded.       //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "SigParamInterface.h"

/**
   -----------------------------------------------------------------------------
   Initialize the SigParamInterface class with a new RooCategory.
   @param newConfigFile - The name of the analysis config file.
   @param newOptions - The job options ("New", "FromFile")
*/
SigParamInterface::SigParamInterface(TString newConfigFile, TString newOptions){
  std::cout << "\nSigParamInterface::Initializing..."
	    << "\n\tconfigFile = " << newConfigFile
	    << "\n\toptions = " << newOptions << std::endl;
  
  // Assign member variables:
  m_configFile = newConfigFile;
    
  m_signalsOK = true;
  m_failedSigParam = "";
  m_sigMap.clear();
  
  // Set the ATLAS Style for plots:
  CommonFunc::SetAtlasStyle();
  
  // Load the configuration for the analysis:
  m_config = new Config(m_configFile);
  
  // Assign output directory, and make sure it exists:
  m_outputDir = Form("%s/%s/DMSigParam", 
		     (m_config->getStr("masterOutput")).Data(),
		     (m_config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Load the SM signal parameterization from file or start from scratch:
  std::vector<TString> sigSMModes = m_config->getStrV("sigSMModes");
  for (int i_SM = 0; i_SM < (int)sigSMModes.size(); i_SM++) {
    if ((newOptions.Contains("FromFile") && loadFile(sigSMModes[i_SM]))
	|| createNew(sigSMModes[i_SM])) {
      std::cout << "SigParamInterface: " << sigSMModes[i_SM] << " ready!"
		<< std::endl;
    }
    else m_signalsOK = false;
  }
  
  // Load the DM signal parameterization from file or start from scratch:
  std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
  for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
    if ((newOptions.Contains("FromFile") && loadFile(sigDMModes[i_DM]))
	|| createNew(sigDMModes[i_DM])) {
      std::cout << "SigParamInterface: " << sigDMModes[i_DM] << " ready!"
		<< std::endl;
    }
    else m_signalsOK = false;
  }
  
  // Also load or create the total SM parameterization:
  if ((newOptions.Contains("FromFile") && loadFile("SM")) || createNew("SM")) {
    std::cout << "SigParamInterface: Total SM signal ready!" << std::endl;
  }
  else m_signalsOK = false;
  
  if (m_signalsOK) {
    std::cout << "SigParamInterface: Successfully initialized!" << std::endl;
  }
  else {
    std::cout << "SigParamInterface: Problem initializing :(" << std::endl;
    std::cout << "\tFailed fits: " << m_failedSigParam << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Check if all of the signals were either loaded successfully or created.
*/
bool SigParamInterface::allSignalsReady() {
  if (!m_signalsOK) {
    std::cout << "SigParamInterface: Problems detected with following signals"
	      << m_failedSigParam << std::endl;
  }
  return m_signalsOK;
}

/**
   -----------------------------------------------------------------------------
   Create a new signal parameterization.
   @param signalType - The type of signal for parameterization.
   @returns - True iff successfully created.
*/
bool SigParamInterface::createNew(TString signalType) {
  std::cout << "SigParamInterface: Creating new signal for "
	    << signalType << std::endl;
  
  bool signalConverged = true;
  SigParam *sp = new SigParam(signalType, m_outputDir);
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    RooDataSet *currDataSet = getData(signalType, i_c);
    sp->addDataSet(m_config->getNum("higgsMass"), i_c, currDataSet, "m_yy");
    
    if (sp->makeSingleResonance(m_config->getNum("higgsMass"), i_c, 
				m_config->getStr("resonancePDF"))) {
      sp->saveAll();
      sp->plotSingleResonance(m_config->getNum("higgsMass"), i_c);
      m_sigMap[signalType] = sp;
    }
    else {
      signalConverged = false;
      m_failedSigParam += signalType + ", ";
    }
  }
  return signalConverged;
}

/**
   -----------------------------------------------------------------------------
   Retrieve the relevant datasets for the signal and category of interest.
   @param signalType - The type of signal for parameterization.
   @param cateIndex - The category index.
   @returns - A RooDataSet with all data points necessary for fitting.
*/
RooDataSet* SigParamInterface::getData(TString signalType, int cateIndex) {
  if (signalType.EqualTo("SM")) {
    RooDataSet *currData = NULL;
    std::vector<TString> sigSMModes = m_config->getStrV("sigSMModes");
    for (int i_SM = 0; i_SM < (int)sigSMModes.size(); i_SM++) {
      DMMassPoints *mp = new DMMassPoints(m_configFile, sigSMModes[i_SM],
					  "FromFile", NULL);
      if (i_SM == 0) currData = mp->getCateDataSet(cateIndex);
      else currData->append(*(mp->getCateDataSet(cateIndex)));
    }
    return currData;
  }
  else {
    DMMassPoints *mp
      = new DMMassPoints(m_configFile, signalType, "FromFile", NULL);
    return mp->getCateDataSet(cateIndex);
  }
}

/**
   -----------------------------------------------------------------------------
   Get a pointer to the SigParam class for the specified signal type.
   @param signalType - The type of signal for parameterization.
*/
SigParam* SigParamInterface::getSigParam(TString signalType) {
  return m_sigMap[signalType];
}

/**
   -----------------------------------------------------------------------------
   Load signal parameterization from file. If that fails, create new.
   @param signalType - The type of signal for parameterization.
   @returns - True iff successfully loaded or created.
*/
bool SigParamInterface::loadFile(TString signalType) {
  SigParam *sp = new SigParam(signalType, m_outputDir);
  if (sp->loadParameterization(m_outputDir, signalType)) {
    m_sigMap[signalType] = sp;
    std::cout << "SigParamInterface: Successful load from file for "
	      << signalType << std::endl;
    return true;
  }
  else {
    std::cout << "SigParamInterface: Failed to load from file for "
	      << signalType << std::endl;
    return createNew(signalType);
  }
}
