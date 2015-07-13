////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMPseudoExp.h                                                       //
//  Main method: DMPseudoExp.cxx                                              //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 23/06/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Package includes:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "DMAnalysis.h"
#include "DMTestStat.h"
#include "RooBernsteinM.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "statistics.h"

TString options;
std::vector<double> numEventsPerCate;

void createPseudoData(RooWorkspace *w, ModelConfig *mc, int seed);
