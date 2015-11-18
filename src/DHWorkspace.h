////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHWorkspace.h                                                       //
//  Class: DHWorkspace.cxx                                                    //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 06/07/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DHWorkspace_h
#define DHWorkspace_h

// Package libraries:
#include "BkgModel.h"
#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "DHAnalysis.h"
#include "DHDataReader.h"
#include "DHTestStat.h"
#include "PERReader.h"
#include "PESReader.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "statistics.h"

class DHWorkspace {

 public:
  
  DHWorkspace(TString newConfigFile, TString newAnalysisType, TString options);
  virtual ~DHWorkspace() {};
  
  bool fitsAllConverged();
  RooWorkspace* getCombinedWorkspace();
  ModelConfig* getModelConfig();
  
 private:
  
  void loadWSFromFile();
  void createNewWS();
  RooWorkspace* createNewCategoryWS();
  void makeNP(TString varName, double setup[4], RooArgSet *&nuisParams,
	      RooArgSet *&constraints, RooArgSet *&globalObs,
	      RooArgSet *&expected);
  void makeShapeNP(TString varnameNP, TString process, double setup[4],
		   RooArgSet *&nuisParams, RooArgSet *&constraints,
		   RooArgSet *&globalObs, RooArgSet *&expected);
  void createAsimovData(RooWorkspace *cateWS, int valMuDH, int valMuSH);
  RooDataSet* createAsimovData(int valMuDH, int valMuSH);
  
  void getTemporarySingleHiggs(RooWorkspace* workspace);
  void plotSingleCateFit(RooWorkspace *cateWS, TString dataset, 
			 TString observableName);
  TString varToName(TString varForm);
  TString funcToName(TString funcForm);

  // Member variables:
  TString m_configFile;
  TString m_anaType;
  TString m_options;
  TString m_outputDir;
  int m_nCategories;
  int m_muNominalSH;
  TString m_dataToPlot;
  
  // Helper classes:
  Config *m_config;
  
  // Updated for each call to createNewCategoryWS():
  int m_currCateIndex;
  TString m_currCateName;
  
  // The Final RooWorkspace and ModelConfig:
  RooWorkspace *m_ws;
  ModelConfig *m_modelConfig;
  
  RooArgSet *m_nuisanceParameters;
  RooArgSet *m_globalObservables;
  RooArgSet *m_observables;
  RooArgSet *m_poi;
  
  RooCategory *m_categories;
  RooSimultaneous *m_combinedPdf;
  
  std::map<std::string, RooDataSet*> m_combData;
  map<string,RooDataSet*> m_combDataAsimov0;
  map<string,RooDataSet*> m_combDataAsimov1;
  
  // Track whether fits converge:
  bool m_allGoodFits;
  
};

#endif

