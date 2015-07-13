////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: BkgModel.h                                                          //
//  Class: BkgModel.cxx                                                       //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/06/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef BkgModel_h
#define BkgModel_h

// C++ includes:
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// ROOT includes:
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TString.h"

// Package includes:
#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooBernsteinM.h"
#include "RooLandau.h"

class BkgModel {
  
 public:
  
  BkgModel(RooRealVar *newObs);
  virtual ~BkgModel() {};
  
  // Accessors:
  void addBkgToCateWS(RooWorkspace *&workspace, RooArgSet *&nuisParams,
		      TString function);
  RooAbsPdf* getBkgPDF(TString function);
  RooRealVar* getMassObservable();
  
  // Mutators:
  void setObservable(RooRealVar *newObs);
  
 private:
  
  // Member methods:
  int getOrderFromFunc(TString function);
  
  // Member variables:
  RooRealVar *m_obs;
  TString m_obsName;
  double m_obsMin;
  double m_obsMax;
};

#endif
