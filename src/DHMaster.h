////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHMaster.h                                                          //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 01/01/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "Config.h"
#include "DHSystematics.h"
#include "DHTestStat.h"
#include "DHToyAnalysis.h"

bool m_isFirstJob;

Config *m_config;

void submitTSViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal);

void submitMLViaBsub(TString exeConfigFile, TString exeOption, 
		     TString exeSignal);

void submitPEViaBsub(TString exeConfigFile, TString exeOption,
		     int exeSeed, int exeToysPerJob, int resonanceMass);
