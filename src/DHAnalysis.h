////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHAnalysis.h                                                        //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 06/08/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef _DHAnalysis_h_
#define _DHAnalysis_h_

// C++ libraries:
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

// ROOT libraries:
#include "TString.h"

// Package libraries:
#include "Config.h"

namespace DHAnalysis {
  
  TString getAnalysisType(Config *config, TString sampleName);
  bool cateHasComponent(Config *config, TString cateName, TString component);
  bool isSMSample(Config *config, TString sampleName);
  bool isDHSample(Config *config, TString sampleName);
  bool isSignalSample(Config *config, TString sampleName);
  bool isWeightedSample(Config *config, TString sampleName);
  int getMediatorMass(TString sampleName);

};

#endif
