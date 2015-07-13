////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DMToyAnalysis.h                                                           //
//  Class: DMToyAnalysis.cxx                                                  //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 24/06/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMToyAnalysis_h
#define DMToyAnalysis_h

#include "CommonHead.h"
#include "CommonFunc.h"
#include "DMToyTree.h"
#include "DMTestStat.h"
#include "RooFitHead.h"
#include "statistics.h"

class DMToyAnalysis {

 public:
  
  DMToyAnalysis(TString newJobName, TString newDMSignal, TString newCateScheme,
		TString newOptions);
  virtual ~DMToyAnalysis() {};
  
  void fillToyHistograms(int muValue, DMToyTree *toyTree);
  void getAsymptoticForm(TString statistic);
  TH1F* getAsymptoticHist();
  TH1F* getGlobsHist(TString paramName, TString fitType, int toyMu);
  TH1F* getNuisHist(TString paramName, TString fitType, int toyMu);
  TH1F* getMuHist(int toyMu);
  TH1F* getStatHist(TString statistic, int toyMu);
  void plotParameter(TString paramName, TString paramType, int toyMu);
  void plotProfiledMu(); 
  void plotTestStat(TString statistic);
  void plotTestStatComparison(TString statistic);
    
 private:
  
  TString printStatName(TString statistic);
  
  // Private member variables:
  TString jobName;
  TString DMSignal;
  TString cateScheme;
  TString options;
  TString outputDir;
  
  // Classes for statistics access:
  DMTestStat *dmts;
  RooWorkspace *workspace;
  
  // Test statistic binning:
  int nBins;
  int binMin;
  int binMax;
  
  // Histograms:
  TH1F *hAsymptotic;
  TH1F *hMuProfiled[2];
  TH1F *hQ0[2];
  TH1F *hQMu[2];
  TH1F *hQMuTilde[2];
  TH1F *hNuisMu0[20][2];
  TH1F *hNuisMu1[20][2];
  TH1F *hNuisMuFree[20][2];
  TH1F *hGlobsMu0[20][2];
  TH1F *hGlobsMu1[20][2];
  TH1F *hGlobsMuFree[20][2];
  
  // Parameter data:
  std::vector<std::string> namesGlobs;
  std::vector<std::string> namesNuisParams;
  int nGlobs;
  int nNuisParams;
  
};

#endif

