////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHCheckJobs.h                                                       //
//  Class: DHCheckJobs.cxx                                                    //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 11/08/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DHCheckJobs_h
#define DHCheckJobs_h

// C++ libraries:
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>

// ROOT libraries:
#include "TString.h"

// Package libraries:
#include "Config.h"
#include "DHAnalysis.h"

class DHCheckJobs {

 public:
  
  // Constructor and destructor:
  DHCheckJobs(TString configFileName);
  virtual ~DHCheckJobs() {};
  
  // Mutators:
  int getNumberToResubmit(TString jobType);
  std::vector<TString> getResubmitList(TString jobType);
  void updateJobStatus(TString jobType);
  
  // Accessors:
  void printResubmitList(TString jobType);
  
 private:
  
  Config *m_config;
  
  std::vector<TString> m_listDHTestStat;
  std::vector<TString> m_listDHMuLimit;
  
};

#endif
