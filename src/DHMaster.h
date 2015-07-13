////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHMaster.h                                                          //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/07/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "BkgModel.h"
#include "CheckJobs.h"
#include "DHAnalysis.h"
#include "DHTestStat.h"
#include "DHToyAnalysis.h"
#include "SigParam.h"

bool isFirstJob;

void makeExe(TString exeName);

void submitWSViaBsub(TString exeJobName, TString exeOption, TString exeSignal,
		     TString exeCateScheme);

void submitTSViaBsub(TString exeJobName, TString exeOption, TString exeSignal);

void submitMLViaBsub(TString exeJobName, TString exeOption, TString exeSignal);

void submitPEViaBsub(TString exeJobName, TString exeOption, TString exeSignal,
		     int exeSeed, int exeToysPerJob);
