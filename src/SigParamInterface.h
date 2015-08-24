////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParamInterface.h                                                 //
//  Class: SigParamInterface.cxx                                              //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/06/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef SigParamInterface_h
#define SigParamInterface_h

// C++ libraries:
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

// ROOT libraries:
#include "TString.h"

// Package libraries:
#include "CommonFunc.h"
#include "Config.h"
#include "DMAnalysis.h"
#include "DMMassPoints.h"
#include "SigParam.h"

class SigParamInterface 
{
  
 public:
  
  // Constructor / destructor:
  SigParamInterface(TString newConfigFile, TString newOptions);
  virtual ~SigParamInterface() {};
  
  // Accessors:
  bool allSignalsReady();
  bool createNew(TString signalType);
  RooDataSet *getData(TString signalType, int cateIndex);
  SigParam* getSigParam(TString signalType);
  bool loadFile(TString signalType);

 private:
  
  // Member variables:
  TString m_configFile;
  TString m_outputDir;
  
  TString m_failedSigParam;
  bool m_signalsOK;
  
  // Access to analysis settings:
  Config *m_config;
 
  std::map<TString,SigParam*> m_sigMap;

};

#endif
