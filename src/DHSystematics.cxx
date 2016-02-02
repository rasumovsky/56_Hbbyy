////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DHSystematics.cxx                                                         //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 01/01/2016                                                          //
//                                                                            //
//  This program loads text files with systematic uncertainties and produces  //
//  outputs for use with the workspaces.                                      //
//                                                                            //
//  Steps:                                                                    //
//    1 Instantiate class DHSystematics()                                     //
//    2 loadSystematicsFile()                                                 //
//    3 setSysToDefaults()                                                    //
//    4 setConstrCenterTypeIncl()                                             //
//    5 printWorkspaceInput()                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHSystematics.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the DHSystematics class.
   @param newConfigFile - The name of the analysis config file.
   @param options - Options for the toy analysis: ForcePlot, CLScan
*/
DHSystematics::DHSystematics(TString newConfigFile, TString options) {
  
  // Load the config file:
  m_config = new Config(newConfigFile);
  TString jobName = m_config->getStr("JobName");
  TString anaType = m_config->getStr("AnalysisType");
  
  printer(Form("DHSystematics::DHSystematics(%s)",newConfigFile.Data()),false);
  
  // set input and output directories:
  m_outputDir = Form("%s/%s/DHSystematics",
		     (m_config->getStr("MasterOutput")).Data(),
		     jobName.Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Clear the sys values and names storage:
  m_maxCategories = 0;
  m_cateNames.clear();
  m_missingSys = false;
  m_sampleNames.clear();
  m_sysValues.clear();
  m_sysNames.clear();
  
  // Systematic information:
  m_constraintType.clear();
  m_centralValue.clear();
  m_sysType.clear();
  m_inclusive.clear();
}

/**
   -----------------------------------------------------------------------------
   Set the names of categories in the analysis.
   @param cateNames - A vector containing the ordered category names.
*/
void DHSystematics::addCategoryNames(std::vector<TString> cateNames) {
  m_cateNames = cateNames;
}

/**
   -----------------------------------------------------------------------------
   Add a systematic uncertainty to the class data.
   @param sample - The name of the sample.
   @param systematic - The name of the systematic variation.
   @param category - The category for the systematic.
   @param up_down - "up" or "down" variation.
   @param systematicValue - The value of the systematic.
*/
void DHSystematics::addSystematic(TString sample, TString systematic,
				  int category, TString up_down, 
				  double systematicValue) {
  //printer(Form("DHSystematics::addSystematic(%s, %s, %d, %s, %f)",
  //	       sample.Data(), systematic.Data(), category, up_down.Data(),
  //	       systematicValue), false);
  
  // Store systematic value:
  m_sysValues[mKey(sample, systematic, category, up_down)] = systematicValue;
  
  // Increase the total number of categories if necessary:
  if (category >= m_maxCategories) m_maxCategories = category + 1;
  
  // Add the name to the list of unique systematics (if not already done):
  addSystematicName(systematic);
}

/**
   -----------------------------------------------------------------------------
   Add a systematic uncertainty name to the class data, avoids duplicates.
   @param systematic - The name of the systematic variation.
*/
void DHSystematics::addSystematicName(TString systematic) {
  bool matched = false;
  for (int i_n = 0; i_n < (int)m_sysNames.size(); i_n++) {
    if (systematic.EqualTo(m_sysNames[i_n])) {
      matched = true;
      break;
    }
  }
  if (!matched) m_sysNames.push_back(systematic);
}

/**
   -----------------------------------------------------------------------------
   Get the inclusive symmetric value of a systematic:
   @param sample - The name of the sample.
   @param systematic - The name of the systematic variation.
   @return - The inclusive, symmetric value of a systematic:
*/
double DHSystematics::getInclSymmetricSystematic(TString sample, 
						 TString systematic) {
  double sysValue = 0.0;
  for (int i_c = 0; i_c < m_maxCategories; i_c++) {
    sysValue += getSymmetricSystematic(sample, systematic, i_c);
  }

  // Then average:
  return (sysValue / ((double)m_maxCategories));
}

/**
   -----------------------------------------------------------------------------
   Get the inclusive value of a systematic uncertainty.
   @param sample - The name of the sample.
   @param systematic - The name of the systematic variation.
   @param up_down - "up" or "down" variation of the systematic.
   @return - The inclusive value of a systematic uncertainty.
*/				  
double DHSystematics::getInclSystematic(TString sample, TString systematic, 
					TString up_down) {
  double sysValue = 0.0;
  for (int i_c = 0; i_c < m_maxCategories; i_c++) {
    sysValue += getSystematic(sample, systematic, i_c, up_down);
  }
  
  // Then average:
  return (sysValue / ((double)m_maxCategories));
}

/**
   -----------------------------------------------------------------------------
   Get a symmetrized systematic uncertainty's value from the class data.
   @param sample - The name of the sample.
   @param systematic - The name of the systematic variation.
   @param category - The category for the systematic.
   @return - The symmetric value of the systematic uncertainty.
*/
double DHSystematics::getSymmetricSystematic(TString sample,
					     TString systematic, int category) {
  int nMissing = 0;
  double valueUp = getSystematic(sample, systematic, category, "up");
  if (m_missingSys) nMissing++;
  double valueDown = getSystematic(sample, systematic, category, "down");
  if (m_missingSys) nMissing++;
  double valueSymmetric;
  if (nMissing == 0) valueSymmetric = ((fabs(valueUp) + fabs(valueDown)) / 2.0);
  else if (nMissing == 1) valueSymmetric = (fabs(valueUp) + fabs(valueDown));
  else valueSymmetric = 0;
  return valueSymmetric;
}

/**
   -----------------------------------------------------------------------------
   Get a systematic uncertainty's value from the class data.
   @param sample - The name of the sample.
   @param systematic - The name of the systematic variation.
   @param category - The category for the systematic.
   @param up_down - "up" or "down" variation.
   @return - The value of the systematic uncertainty.
*/
double DHSystematics::getSystematic(TString sample, TString systematic ,
				    int category, TString up_down) {
  m_missingSys = false;
  double systematicValue = 0.0;
  TString other = (up_down == "up") ? "down" : "up";
  if (m_sysValues.count(mKey(sample, systematic, category, up_down)) > 0) {
    systematicValue = m_sysValues[mKey(sample, systematic, category, up_down)];
  }
  //else if (m_sysValues.count(mKey(sample, systematic, category, other)) > 0) {
  //systematicValue = m_sysValues[mKey(sample, systematic, category, other)];
  //systematicValue *= -1.0;
  //}
  else {
    m_missingSys = true;
    systematicValue = 0.0;
    
    //std::cout << "ERROR! Exiting after call to:" << std::endl;
    //printer(Form("DHSystematics::getSystematic(%s, %s, %d, %s) = %f",
    //		 sample.Data(), systematic.Data(), category, up_down.Data(),
    //		 systematicValue), true);
  }

  return systematicValue;
}

/**
   -----------------------------------------------------------------------------
   Groups several systematic uncertainties together for a given sample. 
   WARNING! The systematics are all assumed to be of the same constraint type,
   with the same systematic type and central value and categorization.
   @param sample - The name of the sample.
   @param groupName - The name of the combined systematic uncertainty.
   @param sysComponents - The component systematic uncertainties that will be 
   added together in quadrature. 
*/
void DHSystematics::groupSyst(TString sample, TString groupName,
			      std::vector<TString> sysComponents) {
  
  // Variables to store information about the grouped systematic:
  TString constraintType = "asym";
  double centralValue = 1.0;
  TString sysType = "yield";
  bool inclusive = false;
  
  // Store the sum of squares systematics values:
  std::map<int,double> sumSqrValsUp; sumSqrValsUp.clear();
  std::map<int,double> sumSqrValsDown; sumSqrValsDown.clear();
  for (int i_c = 0; i_c < m_maxCategories; i_c++) {
    sumSqrValsUp[i_c] = 0.0;
    sumSqrValsDown[i_c] = 0.0;
  }
  
  // Loop over the systematic uncertainties:
  for (int i_s = 0; i_s < (int)sysComponents.size(); i_s++) {
    
    // Grab the systematic data from the first entry into the group:
    if (i_s == 0) {
      constraintType = m_constraintType[sysComponents[i_s]];
      centralValue = m_centralValue[sysComponents[i_s]];
      sysType = m_sysType[sysComponents[i_s]];
      inclusive = m_inclusive[sysComponents[i_s]];
    }
    
    // Loop over categories:
    for (int i_c = 0; i_c < m_maxCategories; i_c++) {
      double currValueUp
	= getSystematic(sample, sysComponents[i_s], i_c, "up");
      double currValueDown
	= getSystematic(sample, sysComponents[i_s], i_c, "down");
      sumSqrValsUp[i_c] += (currValueUp * currValueUp);
      sumSqrValsDown[i_c] += (currValueDown * currValueDown);
    }
  }
  
  // Then at the end, get square root of systematics and add to class data:
  for (int i_c = 0; i_c < m_maxCategories; i_c++) {
    sumSqrValsUp[i_c] = sqrt(sumSqrValsUp[i_c]);
    sumSqrValsDown[i_c] = sqrt(sumSqrValsDown[i_c]);
    addSystematic(sample, groupName, i_c, "up", sumSqrValsUp[i_c]);
    addSystematic(sample, groupName, i_c, "down", sumSqrValsDown[i_c]);
  }
  
  // Save the settings for the grouped systematic:
  setConstrCenterTypeIncl(groupName, constraintType, centralValue, sysType,
			  inclusive);
}

/**
   -----------------------------------------------------------------------------
   Groups several systematic uncertainties together for all samples. 
   WARNING! The systematics are all assumed to be of the same constraint type,
   with the same systematic type and central value and categorization.
   @param groupName - The name of the combined systematic uncertainty.
   @param sysComponents - The component systematic uncertainties that will be 
   added together in quadrature. 
*/
void DHSystematics::groupSystAllSamples(TString groupName,
					std::vector<TString> sysComponents) {
  for (int i_s = 0; i_s < (int)m_sampleNames.size(); i_s++) {
    groupSyst(m_sampleNames[i_s], groupName, sysComponents);
  }
}

/**
   -----------------------------------------------------------------------------
   List the names of samples in this class
*/
std::vector<TString> DHSystematics::listSampleNames() {
  return m_sampleNames;
}

/**
   -----------------------------------------------------------------------------
   List the names of systematic uncertainties in this class
*/
std::vector<TString> DHSystematics::listSystematicNames() {
  return m_sysNames;
}

/**
   -----------------------------------------------------------------------------
   Load the systematics for one sample from file into the class data.
   @param fileName - The name of the file containing systematic uncertainties.
   @param sample - The name of the sample described in the file (ggH, VBF, etc.)
*/
void DHSystematics::loadSystematicsFile(TString fileName, TString sample) {
  printer(Form("loadSystematicsFile(%s, %s)",fileName.Data(), sample.Data()),
	  false);
  
  // Read in file with systematic variations:
  std::ifstream inputFile(fileName);
  if (inputFile.is_open()) {
    int lineIndex = 0;
    std::string line = "";
    
    // Add this sample to the list of stored samples:
    m_sampleNames.push_back(sample);
    
    // Each line represents one systematic shift for one sample:
    while (std::getline(inputFile, line)) {
      lineIndex++;
      if (lineIndex < 4) continue;

      TString systematicName = "";
      TString up_down = "";
      TObjArray *array = ((TString)line).Tokenize(",");
      for (int i_t = 0; i_t < array->GetEntries(); i_t++) {
	TString currElement = ((TObjString*)array->At(i_t))->GetString();
	
	// First element in the string is the syst. name:
	if (i_t == 0) {
	  systematicName = currElement;
	  if (systematicName.Contains("__1up")) {
	    systematicName.ReplaceAll("__1up", "");
	    up_down = "up";
	  }
	  else if (systematicName.Contains("__1down")) {
	    systematicName.ReplaceAll("__1down", "");
	    up_down = "down";
	  }
	}
	// Other elements contain value of syst. in one category (jj, bj, bb):
	else {
	  double sysValue = currElement.Atof();
	  addSystematic(sample, systematicName, i_t-1, up_down, sysValue);
	}
      }
    }
    inputFile.close();
  }
  // Exit and specify error if systematics don't load properly:
  else {
    printer(Form("DHSystematics: Error opening file %s",fileName.Data()), true);
  }
  
}

/**
   -----------------------------------------------------------------------------
   Get the key to systematics map values.
   @param sample - The name of the sample.
   @param systematic - The name of the systematic.
   @param category - The index of the category.
   @return - The key to the m_sysValues map.
*/
TString DHSystematics::mKey(TString sample, TString systematic, int category,
			    TString up_down) {
  return Form("%s_%s_%d_%s", sample.Data(), systematic.Data(), category,
	      up_down.Data());
}

/**
   -----------------------------------------------------------------------------
   Prints a statement (if verbose) and exits (if fatal).
   @param statement - The statement to print.
   @param isFatal - True iff. this should trigger an exit command.
*/
void DHSystematics::printer(TString statement, bool isFatal) {
  if (m_config->getBool("Verbose") || isFatal) {
    std::cout << statement << std::endl;
  }
  if (isFatal) exit(0);
}

/**
   -----------------------------------------------------------------------------
   Prints the input for Andrew's workspace config file.
*/
void DHSystematics::printWorkspaceInput() {
  TString fileName = Form("%s/workspaceConfig.cfg",m_outputDir.Data());
  printer(Form("DHSystematics::printWorkspaceInput() --> %s",fileName.Data()),
	  false);

  ofstream outputFile(fileName);
  
  // First list the nuisance parameters:
  outputFile << "SysSources:		";
  for (int i_n = 0; i_n < (int)m_sysNames.size(); i_n++) {
    outputFile << m_sysNames[i_n] << " ";
  }
  outputFile << "\n" << std::endl;
  
  // Then go over the nuisance parameters:
  for (int i_n = 0; i_n < (int)m_sysNames.size(); i_n++) {
    outputFile << "SysForm_" << m_sysNames[i_n] << ": " << m_sysNames[i_n] 
	       << "[constr=" << m_constraintType[m_sysNames[i_n]] 
	       << ",center=" << m_centralValue[m_sysNames[i_n]] 
	       << ",type=" << m_sysType[m_sysNames[i_n]]
	       << ",";
    
    // Then loop over the samples:
    for (int i_s = 0; i_s < (int)m_sampleNames.size(); i_s++) {
      
      if (i_s != 0)outputFile<< ",";
      
      // For inclusive case, do one overall value (not per category):
      if (m_inclusive[m_sysNames[i_n]]) {
	if ((m_constraintType[m_sysNames[i_n]]).Contains("asym")) {
	  double valueUp = getInclSystematic(m_sampleNames[i_s],
					     m_sysNames[i_n], "up");
	  double valueDown = getInclSystematic(m_sampleNames[i_s], 
					       m_sysNames[i_n], "down");
	  outputFile << "comp=" << m_sampleNames[i_s] << "~" << valueUp << "," 
		     << "compLo=" << m_sampleNames[i_s] << "~" << valueDown;
	}
	else {
	  double value
	    = getInclSymmetricSystematic(m_sampleNames[i_s], m_sysNames[i_n]);
	  outputFile << "comp=" << m_sampleNames[i_s] << "~" << value;
	}
      }
      
      // Otherwise specify different value in each category:
      else {
	// Loop over the categories:
	for (int i_c = 0; i_c < m_maxCategories; i_c++) {
	  if (i_c != 0) outputFile << ",";
	  if ((m_constraintType[m_sysNames[i_n]]).Contains("asym")) {
	    double valueUp = getSystematic(m_sampleNames[i_s],
					   m_sysNames[i_n], i_c, "up");
	    double valueDown = getSystematic(m_sampleNames[i_s], 
					     m_sysNames[i_n], i_c, "down");
	    outputFile << "comp=" << m_sampleNames[i_s] << "_" 
		       << m_cateNames[i_c] << "~" << valueUp << "," 
		       << "compLo=" << m_sampleNames[i_s] << "_" 
		       << m_cateNames[i_c] << "~" << valueDown;
	  }
	  
	  else {
	    double value = getSymmetricSystematic(m_sampleNames[i_s],
						  m_sysNames[i_n], i_c);
	    outputFile << "comp=" << m_sampleNames[i_s] << "_"
		       << m_cateNames[i_c] << "~" << value;
	  }
	}
      }
    }
    outputFile << "]" << std::endl;
  }
  outputFile.close();
}

/**
   -----------------------------------------------------------------------------
   Use default settings for each defined systematic uncertainty (to overwrite 
   by user later... Set to asymmetric uncertainty on the signal yield with a 
   different value in each category with a central value of 1.
*/
void DHSystematics::setSysToDefaults() {
  for (int i_s = 0; i_s < (int)m_sysNames.size(); i_s++) {
    setConstrCenterTypeIncl(m_sysNames[i_s], "asym", 1.0, "yield", false);
  }
}

/**
   -----------------------------------------------------------------------------
   Set the general systematic uncertainty settings.
   @param systematic - The name of the systematic uncertainty.
   @param constraintType - The type of constraint term (asym, gaus, logn)
   @param central value - The central value of the systematic (0 or 1 usually)
   @param systematicType - "yield", "scale", or "res"
   @param inclusive - True iff inclusive uncertainty.
*/
void DHSystematics::setConstrCenterTypeIncl(TString systematic,
					    TString constraintType,
					    double centralValue,
					    TString systematicType,
					    bool inclusive) {
  m_constraintType[systematic] = constraintType;
  m_centralValue[systematic] = centralValue;
  m_sysType[systematic] = systematicType;
  m_inclusive[systematic] = inclusive;
}
