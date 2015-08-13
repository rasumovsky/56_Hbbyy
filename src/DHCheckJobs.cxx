////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHCheckJobs.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 11/08/2015                                                          //
//                                                                            //
//  This class checks to see whether jobs of a particular type have finished. //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHCheckJobs.h"

/**
   Initialize the class for checking job statuses.
   @param configFileName - The name of the config file for the analysis.
*/
DHCheckJobs::DHCheckJobs(TString configFileName) {
  m_config = new Config(configFileName);
  updateJobStatus("DHTestStat");
  updateJobStatus("DHMuLimit");
  return;
}

/**
   Get the number of jobs that need to be resubmitted.
   @param jobType - the type of job (DHTestStat, DHMuLimit).
   @returns - The number of jobs that failed the first attempt.
*/
int DHCheckJobs::getNumberToResubmit(TString jobType) {
  std::vector<TString> currList = getResubmitList(jobType);
  return (int)currList.size();
}

/**
   Get the list of signal points that must be resubmitted.
   @param jobType - the type of job (DHTestStat, DHMuLimit).
   @returns - A vector of the signal point names.
*/
std::vector<TString> DHCheckJobs::getResubmitList(TString jobType) {
  std::vector<TString> result; result.clear();
  if (jobType.EqualTo("DHTestStat")) result = m_listDHTestStat;
  else if (jobType.EqualTo("DHMuLimit")) result = m_listDHMuLimit;
  return result;
}

/**
   Update the status of jobs from a particular program.
   @param jobType - the type of job (DHTestStat, DHMuLimit)
   @returns - void. Updates the lists of failed jobs.
*/
void DHCheckJobs::updateJobStatus(TString jobType) {
  
  // Save names of failed jobs:
  m_listDHTestStat.clear();
  m_listDHMuLimit.clear();
  
  // Then loop over submissions to see whether output files exist:
  std::vector<TString> sigDHModes = m_config->getStrV("sigDHModes");
  for (int i_DH = 0; i_DH < (int)sigDHModes.size(); i_DH++) {
    TString currDHSignal = sigDHModes[i_DH];
    
    TString fileName;
    TString fullName;
    
    if (jobType.EqualTo("DHTestStat")) {
      fileName = Form("CL_values_%s.txt", currDHSignal.Data());
      fullName = Form("%s/%s/DHTestStat/CL/%s", 
		      (m_config->getStr("masterOutput")).Data(), 
		      (m_config->getStr("jobName")).Data(), fileName.Data());
    }
    
    if (jobType.EqualTo("DHMuLimit")) {
      fileName = Form("text_CLs_%s.txt", currDHSignal.Data());
      fullName = Form("%s/%s/DHMuLimit/single_files/%s", 
		      (m_config->getStr("masterOutput")).Data(), 
		      (m_config->getStr("jobName")).Data(), fileName.Data());
    }

    // Then test the existence of the file.
    std::ifstream testFile(fullName);
    if (!testFile) {
      if (jobType.EqualTo("DHTestStat")) {
	m_listDHTestStat.push_back(currDHSignal);
      }
      if (jobType.EqualTo("DHMuLimit")) {
	m_listDHMuLimit.push_back(currDHSignal);
      }
    }
  }
}

/**
   Print the list and number of failed jobs.
   @param jobType - the type of job (DHTestStat, DHMuLimit).
   @returns void. Prints to the terminal.
*/
void DHCheckJobs::printResubmitList(TString jobType) {
  // Get the relevant list of failed jobs:
  std::vector<TString> currList = getResubmitList(jobType);
  std::cout << "Failed to make the following " << jobType << " files ("
	    << currList.size() << " total)" << std::endl;
  for (int i_f = 0; i_f < (int)currList.size(); i_f++) {
    std::cout << "\t" << currList[i_f] << std::endl;
  }
}
