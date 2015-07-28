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

#include "BkgModel.h"
#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "CommonFunc.h"
#include "statistics.h"
#include "PESReader.h"
#include "PERReader.h"
#include "DHAnalysis.h"
#include "DHTestStat.h"
#include "DHDataReader.h"

class DHWorkspace {

 public:
  
  DHWorkspace(TString jobName, TString cateScheme, TString options);
  virtual ~DHWorkspace() {};
  
  bool fitsAllConverged();
  RooWorkspace* getCombinedWorkspace();
  ModelConfig* getModelConfig(TString analysisType);
  
 private:
  
  void loadWSFromFile();
  void createNewModel();
  void createNewWS();
  RooWorkspace* createNewCategoryWS();
  void makeNP(TString varName, double setup[4], RooArgSet *&nuisParams,
	      RooArgSet *&constraints, RooArgSet *&globalObs,
	      RooArgSet *&expected);
  void makeShapeNP(TString varnameNP, TString process, double setup[4],
		   RooArgSet *&nuisParams, RooArgSet *&constraints,
		   RooArgSet *&globalObs, RooArgSet *&expected);
  void plotSingleCateFit(RooWorkspace *cateWS, TString dataset, 
			 TString observableName);
  
  // Member variables:
  TString m_jobName;
  TString m_cateScheme;
  TString m_options;
  TString m_outputDir;

  int m_muNominalSH;
  TString m_dataToPlot;
  
  // Helper classes:
  PESReader *pes;
  PERReader* per;
  
  // Updated for each call to createNewModel():
  TString currAna;
  int currNCategories;
  
  // Updated for each call to createNewCategoryWS():
  int currCateIndex;
  TString currCateName;
  
  // The Final RooWorkspace and ModelConfig:
  RooWorkspace *m_combinedWS;
  std::map<TString, ModelConfig*> m_mConfig;
  
  // Track whether fits converge:
  bool m_allGoodFits;
  
};

#endif

