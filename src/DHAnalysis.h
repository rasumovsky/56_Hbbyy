////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHAnalysis.h                                                        //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/06/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef _DHAnalysis_h_
#define _DHAnalysis_h_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

#include "TString.h"

namespace DHAnalysis {
  
  ////////////////////////////////////////
  //         GLOBAL PARAMETERS          //
  ////////////////////////////////////////
  
  // Set True for final analysis on data:
  const bool doBlind = true;
  
  const TString analysisName = "DHAnalysis";
  
  // Luminosity in pb-1:
  const double analysisLuminosity = 10000;
  
  const double higgsMass = 125.09;// GeV
  
  const double DHMyyRangeLo = 105.0;// GeV
  const double DHMyyRangeHi = 160.0;// GeV

  const double DHMyybbRangeLo = 200;// GeV
  const double DHMyybbRangeHi = 900;// GeV

  const int nAnalysisTypes = 2;
  const TString analysisTypes[nAnalysisTypes] = {"NonRes", "Res"};
  
  const int nSMModes = 6;
  const TString sigSMModes[nSMModes] = {"ggH","VBF","WH","ZH","bbH","ttH"};
  
  const int nDHModes = 2;
  const TString sigDHModes[nDHModes] = {"NonRes","ResMx300"};
  
  const int nMCProcesses = 1;
  const TString MCProcesses[nMCProcesses] = {"gg_gjet"};
  
  ////////////////////////////////////////
  //    INPUT AND OUTPUT DIRECTORIES    //
  ////////////////////////////////////////
  
  // Location of global input files:
  const TString masterInput
    = "/afs/cern.ch/work/a/ahard/files_Hbbgg/GlobalInputs";
  // Location of output directory:
  const TString masterOutput
    = "/afs/cern.ch/work/a/ahard/files_Hbbgg/FullAnalysis";
  
  // Location of this software package:
  const TString packageLocation
    = "/afs/cern.ch/user/a/ahard/analysis/56_Hbbgg";
  
  // Holding location of cluster job files:
  const TString clusterFileLocation = "/afs/cern.ch/work/a/ahard/jobfiles";
  
  // Location of input TTrees:
  const TString resInput = "/afs/cern.ch/user/a/ahard/work_directory/files_Hbbgg/GlobalInputs/inputs_8TeV/femtotuple.root";
  
  const TString nonResInput = "/afs/cern.ch/user/a/ahard/work_directory/files_Hbbgg/GlobalInputs/inputs_8TeV/femtotuple_gg_mass.root";
  
  const TString singleHiggsList = "/afs/cern.ch/user/a/ahard/work_directory/files_Hbbgg/GlobalInputs/inputs_8Te/MCList_8TeV_singleHiggs.txt";
  
  const TString resData = "/afs/cern.ch/user/a/ahard/work_directory/files_Hbbgg/GlobalInputs/inputs_8TeV/data12_paper.root";
  
  // Locations of systematic uncertainty files:
  const TString fileNamePESValues = "";
  const TString fileNamePERValues = "";

  ////////////////////////////////////////
  //        FOR JOB SUBMISSION          //
  ////////////////////////////////////////
  
  const TString exeWorkspace = "DHWorkspaceWrapper";
  const TString jobScriptWorkspace = "scripts/jobFileWorkspace.sh";
  
  const TString exeTestStat = "DHTestStatWrapper";
  const TString jobScriptTestStat = "scripts/jobFileTestStat.sh";
  
  const TString exeMuLimit = "DHMuLimit";
  const TString jobScriptMuLimit = "scripts/jobFileMuLimit.sh";
  
  const TString exePseudoExp = "DHPseudoExp";
  const TString jobScriptPseudoExp = "scripts/jobFilePseudoExp.sh";
  
  ////////////////////////////////////////
  //           MEMBER FUNCTIONS         //
  ////////////////////////////////////////
  
  bool cateHasComponent(TString cateName, TString component);
  TString cateIndexToName(TString cateScheme, TString analysisType, 
			  int cateIndex);
  int cateNameToIndex(TString cateName);
  TString cateToBkgFunc(TString cateName);
  std::vector<TString> getFitComponents(TString cateName);
  int getNumCategories(TString cateScheme, TString analysisType);  
  bool isSMSample(TString sampleName);
  bool isDHSample(TString sampleName);
  bool isSignalSample(TString sampleName);
  bool isWeightedSample(TString sampleName);
  bool isNonResCate(TString cateName);
  bool isResCate(TString cateName);
  
  TString getAnalysisType(TString signalName);
  int getMediatorMass(TString sampleName);

};

#endif
