////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: CheckJobs.h                                                         //
//  Class: CheckJobs.cxx                                                      //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/06/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef CheckJobs_h
#define CheckJobs_h

#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include "TString.h"

#include "DMAnalysis.h"

class CheckJobs {

 public:
  
  // Constructor and destructor:
  CheckJobs(TString newJobName);
  virtual ~CheckJobs() {};
  
  // Mutators:
  int getNumberToResubmit(TString jobType);
  std::vector<TString> getResubmitList(TString jobType);
  void updateJobStatus(TString jobType);
  
  // Accessors:
  void printResubmitList(TString jobType);
  
 private:
    
  TString jobName;
  std::vector<TString> listDMWorkspace;
  std::vector<TString> listDMTestStat;
  std::vector<TString> listDMMuLimit;
  
};

#endif
