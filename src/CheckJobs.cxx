////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: CheckJobs.cxx                                                       //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/06/2015                                                          //
//                                                                            //
//  This class checks to see whether jobs of a particular type have finished. //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CheckJobs.h"

using namespace DMAnalysis;
/**
   Initialize the class for checking job statuses.
   @param newJobName - the name of the analysis job.
*/
CheckJobs::CheckJobs(TString newJobName) {
  jobName = newJobName;
  updateJobStatus("DMWorkspace");
  updateJobStatus("DMTestStat");
  updateJobStatus("DMMuLimit");
  return;
}

/**
   Get the number of jobs that need to be resubmitted.
   @param jobType - the type of job (DMWorkspace, DMTestStat).
   @returns - the number of jobs that failed the first attempt.
*/
int CheckJobs::getNumberToResubmit(TString jobType) {
  std::vector<TString> currList = getResubmitList(jobType);
  return (int)currList.size();
}

/**
   Get the list of signal points that must be resubmitted.
   @param jobType - the type of job (DMWorkspace, DMTestStat).
   @returns - a vector of the signal point names.
*/
std::vector<TString> CheckJobs::getResubmitList(TString jobType) {
  std::vector<TString> result; result.clear();
  if (jobType.EqualTo("DMWorkspace")) {
    result = listDMWorkspace;
  }
  else if (jobType.EqualTo("DMTestStat")) {
    result = listDMTestStat;
  }
  else if (jobType.EqualTo("DMMuLimit")) {
    result = listDMMuLimit;
  }
  return result;
}

/**
   Update the status of jobs from a particular program.
   @param jobType - the type of job (DMWorkspace, DMTestStat)
   @returns - void. Updates the lists of failed jobs.
*/
void CheckJobs::updateJobStatus(TString jobType) {
  
  // Save names of failed jobs:
  listDMWorkspace.clear();
  listDMTestStat.clear();
  listDMMuLimit.clear();
  
  // Then loop over submissions to see whether output files exist:
  for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
    
    TString currDMSignal = sigDMModes[i_DM];
    
    TString fileName;
    TString fullName;
    if (jobType.EqualTo("DMWorkspace")) {
      fileName = Form("workspaceDM_%s.root", currDMSignal.Data());
      fullName = Form("%s/%s/DMWorkspace/rootfiles/%s", masterOutput.Data(),
		      jobName.Data(), fileName.Data());
    }
    if (jobType.EqualTo("DMTestStat")) {
      fileName = Form("CL_values_%s.txt", currDMSignal.Data());
      fullName = Form("%s/%s/DMTestStat/CL/%s", masterOutput.Data(), 
		      jobName.Data(), fileName.Data());
    }
    
    if (jobType.EqualTo("DMMuLimit")) {
      fileName = Form("text_CLs_%s.txt", currDMSignal.Data());
      fullName = Form("%s/%s/DMMuLimit/single_files/%s", masterOutput.Data(), 
		      jobName.Data(), fileName.Data());
    }

    // Then test the existence of the file.
    std::ifstream testFile(fullName);
    if (!testFile) {
      if (jobType.EqualTo("DMWorkspace")) {
	listDMWorkspace.push_back(currDMSignal);
      }
      if (jobType.EqualTo("DMTestStat")) {
	listDMTestStat.push_back(currDMSignal);
      }
      if (jobType.EqualTo("DMMuLimit")) {
	listDMMuLimit.push_back(currDMSignal);
      }
    }
  }
}

/**
   Print the list and number of failed jobs.
   @param jobType - the type of job (DMWorkspace, DMTestStat).
   @returns void. Prints to the terminal.
*/
void CheckJobs::printResubmitList(TString jobType) {
  // Get the relevant list of failed jobs:
  std::vector<TString> currList = getResubmitList(jobType);
  std::cout << "Failed to make the following " << jobType << " files ("
	    << currList.size() << " total)" << std::endl;
  for (int i_f = 0; i_f < (int)currList.size(); i_f++) {
    std::cout << "\t" << currList[i_f] << std::endl;
  }
}
