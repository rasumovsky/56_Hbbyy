////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHDataReader.h                                                      //
//  Class: DHDataReader.cxx                                                   //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 09/07/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DHDataReader_h
#define DHDataReader_h

// C++ includes:
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>

// ROOT includes:
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"

// Package includes:
#include "Config.h"
#include "DHAnalysis.h"
#include "RooFitHead.h"

class DHDataReader {

 public:
  
  // Constructor and destructor:
  DHDataReader(TString configFile, RooRealVar *observable);
  virtual ~DHDataReader() {};
  
  // Mutators:
  RooDataSet* loadNonResData(TString cateName);
  RooDataSet* loadResData(TString cateName);
  void setMassObservable(RooRealVar *observable);
  
 private:
  
  Config *m_config;
  RooRealVar *m_observable;
  
};

#endif
