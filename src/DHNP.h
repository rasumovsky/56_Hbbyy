////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DHNP.h                                                                    //
//  Class: DHNP.cxx                                                           //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 23/02/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DHNP_h
#define DHNP_h

// Package libraries:
#include "CommonHead.h"
#include "RooFitHead.h"

class DHNP {

 public:
  
  DHNP(TString nuisName);
  virtual ~DHNP() {};
  
  // Public accessor methods:
  TString getName();
  RooRealVar *getNuis(TString fit);
  double getPoIUncertainty();

  // Public mutator methods:
  void setNuis(TString fit, RooRealVar *nuisPar);
  void setPoIUncertainty(double uncertainty);
  
 private:
  
  // Private member variables:
  TString m_nuisName;
  std::map<TString,RooRealVar*> m_nuis;
  double m_uncertaintyOnPoI;

};

#endif

