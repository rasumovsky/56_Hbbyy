////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParam.h                                                          //
//  Class: SigParam.cxx                                                       //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 01/11/2016                                                          //
//                                                                            //
//  Accessors access class data without modifying the member objects, while   //
//  mutators modify the state of the class (and also sometimes return data.   //
//                                                                            //
//  Private accessors and mutators are marked as such because they depend on  //
//  a number of internal settings that are hidden to the user. All of the     //
//  functionalities are available through proper use of the public methods.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef SigParam_h
#define SigParam_h

// C++ includes:
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

// ROOT includes:
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include "TRegexp.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

// ROOT math headers:
#include <Math/QuantFuncMathCore.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <Math/Factory.h>

// RooFit headers:
#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooAbsRealLValue.h>
#include <RooAddition.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooBernstein.h>
#include <RooCategory.h>
#include <RooCBShape.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooExponential.h>
#include <RooExtendPdf.h>
#include <RooFitResult.h>
#include <RooFormula.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
#include <RooLandau.h>
#include <RooMinimizer.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooNumIntConfig.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooProfileLL.h>
#include <RooRandom.h>
#include <RooRealSumPdf.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>
#include <RooVoigtian.h>

#include "RooStats/AsymptoticCalculator.h"

class SigParam {
  
 public:
  
  // Constructor and destructor:
  SigParam(TString signalType, TString directory);
  virtual ~SigParam() {};
  
  //----------Public Accessors----------//
  bool addSigToWS(RooWorkspace *&workspace, int cateIndex);
  bool addSigToWS(RooWorkspace *&workspace, double resonanceMass,
		  int cateIndex);
  double calculateStdDev(double resonanceMass, int cateIndex);
  TF1 *createTF1FromParameterization(TString varName, int cateIndex,
				     double xMin, double xMax);
  bool dataExists(double resonanceMass, int cateIndex);
  double extendedTerm();
  double generatedDataNorm();
  TString getKey(double resonanceMass, int cateIndex);
  double getMeanOrStdDev(TString value, double resonanceMass, int cateIndex);
  double getMeanOrStdDevInData(TString value, double resonanceMass,
			       int cateIndex);
  int getNParamsForVar(TString varName);
  double getParameterError(TString paramName, double resonanceMass,
			   int cateIndex);
  double getParameterError(TString paramName, int cateIndex);
  double getParameterValue(TString paramName, double resonanceMass, 
			   int cateIndex);
  double getParameterValue(TString paramName, int cateIndex);
  double getParameterizedValue(TString paramName, double resonanceMass,
			       int cateIndex);
  TString getParamState(TString paramName);
  RooAbsPdf* getResonance(int cateIndex);
  RooAbsPdf* getSingleResonance(double resonanceMass, int cateIndex);
  double getTestStat(TString statistic, int cateIndex);
  double getTestStat(TString statistic, double resonanceMass, int cateIndex);
  double getTestStatLatest(TString statistic);
  std::vector<TString> getVariableNames(double resonanceMass, int cateIndex);
  TString getVarParameterization(TString varName);
  RooWorkspace* getWorkspace();
  double getYieldErrorInCategory(double resonanceMass, int cateIndex);
  double getYieldErrorTotal(double resonanceMass);
  double getYieldInCategory(double resonanceMass, int cateIndex);
  double getYieldInWindow(double resonanceMass, int cateIndex, double obsMin,
			  double obsMax);
  double getYieldTotalInWindow(double resonanceMass, double obsMin, 
			       double obsMax);
  double getYieldTotal(double resonanceMass);
  std::vector<TString> listParamsForVar(TString varName);
  std::vector<double> massPointsForCategory(int cateIndex);
  std::vector<TString> variablesForFunction(TString function);
  
  //----------Public Mutators----------//
  void addMResSystematic(TString nameMResSys);
  void addMResSystematics(std::vector<TString> nameMResSys);
  void addMScaleSystematic(TString nameMScaleSys);
  void addMScaleSystematics(std::vector<TString> nameMScaleSys);
  void addDataSet(double resonanceMass, int cateIndex, RooDataSet* dataSet,
		  TString observableName);
  void addDataTree(double resonanceMass, int cateIndex, TTree *dataTree,
		   TString massBranchName, TString weightBranchName);
  void addMassPoint(double resonanceMass, int cateIndex, double diphotonMass,
		    double eventWeight);
  std::vector<double> doBiasTest(double resonanceMass, int cateIndex,
				 TString dataType, int seed);
  void doBinnedFit(bool doBinned, double geVPerBin = 1.0);
  bool generateAndFitData(double resonanceMass, int cateIndex, TString dataType,
			  int seed = 1);
  RooDataSet* generateData(double resonanceMass, int cateIndex,
			   TString dataType, int seed = 1);
  bool loadParameterization(TString directory, TString signalType);
  bool makeAllParameterizations(TString function);
  bool makeCategoryParameterization(int cateIndex, TString function);
  bool makeSingleResonance(double resonanceMass, int cateIndex,
			   TString function);
  void makeYieldParameterization(int cateIndex);
  void nameTheCategories(std::vector<TString> cateNames);
  void plotCategoryResonances(int cateIndex);
  void plotSingleResonance(double resonanceMass, int cateIndex, 
			   TString dataType = "");
  void printResTable(double resonanceMass);
  void plotYields(int cateIndex);
  void saveAll();
  void saveParameterization();
  void saveParameterList();
  void saveYieldList();
  void setDirectory(TString directory);
  void setLogYAxis(bool useLogYAxis);
  void setMassWindowSize(double fraction);
  void setMassWindowFixed(bool fixWindow, double windowMin, double windowMax);
  void setParamState(TString paramName, TString valueAndRange);
  void setPlotATLASLabel(TString atlasLabel);
  void setPlotFormat(TString fileFormat);
  void setPlotLuminosity(TString lumiLabel);
  void setPlotXAxisTitle(TString xAxisTitle);
  void setRatioPlot(bool doRatioPlot, double ratioMin, double ratioMax);
  void setResMassConstant(bool setConstant, double resonanceMass);
  void setResMassConstant(bool setConstant);
  void setSignalType(TString signalType);
  void setVarParameterization(TString varName, TString function);
  void useCommonCBGAMean(bool sameCBGAMean);
  void verbosity(bool beVerbose);
  
 private:
  
  //----------Private Accessors----------//
  bool equalMasses(double massValue1, double massValue2);
  bool functionIsDefined(TString function);
  double massIntToDouble(int massInteger);
  int massDoubleToInt(double resonanceMass);
  double normalizationError(RooAbsData *dataSet);

  //----------Private Mutators----------//
  void addVariable(TString paramName, int cateIndex);
  void binTheData(TString unbinnedDataName, double resonanceMass, 
		  int cateIndex);
  void binSingleDataSet(TString unbinnedName, TString binnedName,
			double resonanceMass, int cateIndex);
  RooFitResult* fitResult(int cateIndex, TString dataType = "", 
			  TString option = "");
  RooFitResult* fitResult(double resonanceMass, int cateIndex, 
			  TString dataType = "", TString option = "");
  int getNCategories();
  void parameterizeFunction(TString function, double mRegularized,
			    double mResonance, int cateIndex,
			    bool parameterized);
  void parameterizeVar(TString varName, double mRegularized, double mResonance,
		       int cateIndex, bool parameterized);
  RooDataSet* plotData(RooAbsData *data, RooRealVar *observable);
  TGraphErrors* plotSubtraction(RooAbsData *data, RooAbsPdf *pdf, 
				RooRealVar *observable, double xBins,
				double &chi2Prob);
  TGraphErrors* plotDivision(RooAbsData *data, RooAbsPdf *pdf, 
			     RooRealVar *observable, double xBins,
			     double &chi2Prob);
  double regularizedMass(double resonanceMass);
  void resonanceCreator(double resonanceMass, int cateIndex, TString function);
  void setParamsConstant(RooAbsPdf* pdf, bool isConstant);
  
  // Member variables:
  int m_nCategories;
  TString m_signalType;
  TString m_directory;
  TString m_fileFormat;

  // Objects for fitting:
  RooWorkspace *m_ws;
  std::map<int,RooCategory*> m_cat;
  
  // Store fit initial values and ranges and parameterizations:
  std::map<TString,TString> m_paramState;
  std::map<TString,TString> m_varParameterization;
  
  // Systematic: 
  TString m_listMRS;
  TString m_listMSS;
  
  // Yield data for plotting:
  std::map<int,TF1*> m_yieldFunc;
  std::map<int,TGraphErrors*> m_yieldGraph;
  
  // List to track imported datasets:
  std::vector<std::pair<double,int> > m_massCatePairs;
  
  // True iff fits should be binned.
  bool m_binned;
  //int m_nBinsPerGeV;
  int m_geVPerBin;
  double m_windowFraction;
  double m_windowMin;
  double m_windowMax;
  bool m_fixWindow;
  
  // Fit parameter options:
  bool m_sameCBGAMean;
  
  // Fit result information:
  double m_currChi2;
  double m_currNLL;
  double m_currExtendVal;
  double m_generatedDataNorm;
  std::map<TString,double> m_testStats;
  
  // Plot options:
  TString m_atlasLabel;
  std::vector<TString> m_cateNames;
  TString m_currFunction;
  bool m_doRatioPlot;
  TString m_lumiLabel;
  double m_ratioMin;
  double m_ratioMax;
  bool m_useLogYAxis;
  TString m_xAxisTitle;
  
  // A bool to control how spammy the tool is:
  bool m_verbose;
  
};

#endif
