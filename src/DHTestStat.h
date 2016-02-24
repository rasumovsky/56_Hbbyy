////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHTestStat.h                                                        //
//  Class: DHTestStat.cxx                                                     //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/11/2015                                                          //
//                                                                            //
//  This class allows the user to calculate p0 and CLs based on an input      //
//  workspace.                                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DHTestStat_h
#define DHTestStat_h

// Package libraries:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "DHWorkspace.h"
#include "RooFitHead.h"
#include "statistics.h"

class DHTestStat {
  
 public:
  
  DHTestStat(TString newConfigFile, TString newOptions,
	     RooWorkspace *newWorkspace);
  virtual ~DHTestStat() {};
  
  double accessValue(TString testStat, bool observed, int N);
  void calculateNewCL();
  void calculateNewP0();
  void clearData();
  void clearFitParamSettings();
  RooDataSet* createPseudoData(int seed, int valPoI, bool fixPoI);
  bool fitsAllConverged();
  double functionQMu(double x);
  double functionQMuTilde(double x, double asimovTestStat);
  double getCLFromCLs(double CLs);
  double getCLsFromCL(double CL);
  double getCLFromQMu(double qMu, double N);
  double getCLsFromQMu(double qMu, double N);
  double getFitNLL(TString datasetName, double muVal, bool fixMu,
		   double &profiledMu, bool resetParams=true);
  std::map<std::string,double> getGlobalObservables();
  std::vector<double> getNEventsToys();
  std::map<std::string,double> getNuisanceParameters();
  std::map<std::string,double> getParameters();
  double getP0FromQ0(double q0);
  double getPbFromN(double N);
  double getPbFromQMu(double qMu, double sigma, double mu);
  double getPMuFromQMu(double qMu);
  double getQ0FromNLL(double nllMu0, double nllMuHat, double muHat);
  double getQMuFromNLL(double nllMu, double nllMuHat, double muHat,
		       double muTest);
  double getQMuTildeFromNLL(double nllMu, double nllMu0, double nllMuHat,
			    double muHat, double muTest);
  double graphIntercept(TGraph *graph, double valueToIntercept);
  void loadStatsFromFile();
  TGraph *nllScanGraph();
  void resetParamsAfterFit(bool doResetParamsAfterFit);
  void saveSnapshots(bool doSaveSnapshot);
  std::map<int,double> scanNLL(TString scanName, TString datasetName,
			       TString varToScan,
			       std::vector<TString> varsToFix);
  void setPlotAxis(bool useLogScale, double yMin, double yMax, 
		   double GeVPerBin);
  void setPlotDirectory(TString directory);
  void setParam(TString paramName, double paramVal, bool doSetConstant);
  void setSubPlot(TString ratioOrSubtraction);
  RooWorkspace *theWorkspace();
  ModelConfig *theModelConfig();
  
 private:
  
  TString getKey(TString testStat, bool observed, int N);
  bool mapValueExists(TString mapKey);
  TString nameOfVar(TString varForm);
  TGraphErrors* plotDivision(TString dataName, TString pdfName, TString obsName,
			     double xMin, double xMax, double xBins);
  void plotFits(TString fitType, TString datasetName);
  TGraphErrors* plotSubtract(TString dataName, TString pdfName,
			     TString obsName, double xMin, double xMax,
			     double xBins);
  void printer(TString statement, bool isFatal);
  void printSet(TString setName, RooArgSet* set);
  void storeParams(RooArgSet *set, std::map<std::string,double>& map);

  // From the initialization:
  TString m_anaType;    // The analysis type ("Res", "NonRes").
  TString m_DHSignal;   // The specific signal under study.
  TString m_jobName;    // The name of the group of jobs (for I/O purposes).
  TString m_options;    // Job options.
  TString m_outputDir;  // The output directory for statistics.
  
  TString m_dataForObsQ0;  // The dataset for expected results.
  TString m_dataForObsQMu; // The dataset for observed results.
  TString m_dataForExpQ0;  // The dataset for expected results.
  TString m_dataForExpQMu; // The dataset for observed results.
  
  TString m_plotDir;    // The output directory for plots (not set by default).
  bool m_doSaveSnapshot;// Option to save snapshots of PoI and nuisance params.
  bool m_doPlot;        // Sets whether or not to plot fit results.
  
  // Pointer to the input file:
  TFile *inputFile;
  
  // The configuration of the analysis:
  Config *m_config;

  // Check whether all fits successful:
  bool m_allGoodFits;
  bool m_doResetParamsAfterFit;

  // The workspace for the fits:
  RooWorkspace *m_workspace;
  ModelConfig *m_mc;
  TGraph *m_graphNLL;
  
  // Store the calculated values:
  std::map<TString,double> m_calculatedValues;
  
  // Store fit parameters from NLL calculation:
  std::map<std::string,double> m_mapGlobs;
  std::map<std::string,double> m_mapNP;
  std::map<std::string,double> m_mapPars;
  
  // Store number of toy MC events generated:
  std::vector<double> m_numEventsPerCate;
  
  // In case special parameter settings are used for a fit:
  std::map<TString,double> m_paramValToSet;
  std::map<TString,bool> m_paramConstToSet;
  
  // Plot settings:
  bool m_useLogScale;
  double m_yMin;
  double m_yMax;
  double m_geVPerBin;
  TString m_ratioOrSubtraction;
};

#endif

