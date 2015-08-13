////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParam.h                                                          //
//  Class: SigParam.cxx                                                       //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/06/2015                                                          //
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
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"

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
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooNumIntConfig.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooProfileLL.h>
#include <RooRealSumPdf.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>

class SigParam {
  
 public:
  
  // Constructor and destructor:
  SigParam(TString signalType, TString directory);
  virtual ~SigParam() {};
  
  //----------Public Accessors----------//
  bool addSigToWS(RooWorkspace *&workspace, int cateIndex);
  bool addSigToWS(RooWorkspace *&workspace, double resonanceMass,
		  int cateIndex);
  TString getKey(double resonanceMass, int cateIndex);
  double getParameterError(TString paramName, double resonanceMass,
			   int cateIndex);
  double getParameterError(TString paramName, int cateIndex);
  double getParameterValue(TString paramName, double resonanceMass, 
			   int cateIndex);
  double getParameterValue(TString paramName, int cateIndex);
  TString getParamState(TString paramName);
  RooAbsPdf* getResonance(int cateIndex);
  RooAbsPdf* getSingleResonance(double resonanceMass, int cateIndex);
  std::vector<TString> getVariableNames(double resonanceMass, int cateIndex);
  RooWorkspace* getWorkspace();
  double getYieldInCategory(double resonanceMass, int cateIndex);
  double getYieldTotal(double resonanceMass);
  
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
  void doBinnedFit(bool doBinned, double nBinsPerGeV = 1.0);
  bool loadParameterization(TString directory, TString signalType);
  bool makeAllParameterizations(TString function);
  bool makeCategoryParameterization(int cateIndex, TString function);
  bool makeSingleResonance(double resonanceMass, int cateIndex,
			   TString function);
  void makeYieldParameterization(int cateIndex);
  void plotCategoryResonances(int cateIndex);
  void plotSingleResonance(double resonanceMass, int cateIndex);
  TGraphErrors* plotSubtraction(RooAbsData *data, RooAbsPdf *pdf, double xMin,
				double xMax, double xBins);
  TGraphErrors* plotDivision(RooAbsData *data, RooAbsPdf *pdf, double xMin,
			     double xMax, double xBins);
  void plotYields(int cateIndex);
  void saveAll();
  void saveParameterization();
  void saveParameterList();
  void saveYieldList();
  void setDirectory(TString directory);
  void setLogYAxis(bool useLogYAxis);
  void setParamState(TString paramName, TString valueAndRange);
  void setPlotFormat(TString fileFormat);
  void setRatioPlot(bool doRatioPlot, double ratioMin, double ratioMax);
  void setSignalType(TString signalType);
  void useCommonCBGAMean(bool sameCBGAMean);
  
 private:
  
  //----------Private Accessors----------//
  std::vector<int> categoriesForMass(double resonanceMass);
  bool dataExists(double resonanceMass, int cateIndex);
  bool equalMasses(double massValue1, double massValue2);
  double massIntToDouble(int massInteger);
  int massDoubleToInt(double resonanceMass);
  std::vector<double> massPointsForCategory(int cateIndex);

  //----------Private Mutators----------//
  void addFitParameter(TString paramName, int cateIndex);
  void binTheData(TString unbinnedDataName, int cateIndex);
  void binSingleDataSet(TString unbinnedDataName, TString binnedDataName);
  RooFitResult* fitResult(int cateIndex);
  RooFitResult* fitResult(double resonanceMass, int cateIndex);
  int getNCategories();
  double regularizedMass(double resonanceMass);
  void resonanceCreator(double resonanceMass, int cateIndex, TString function);
  void setParamsConstant(RooAbsPdf* pdf, bool isConstant);
    
  // Member variables:
  int m_nCategories;
  TString m_signalType;
  TString m_directory;
  TString m_fileFormat;

  // Objects for fitting:
  RooRealVar *m_yy;
  RooRealVar *m_wt;
  RooRealVar *m_mResonance;
  RooWorkspace *m_ws;
  RooCategory *m_cat;
  
  // Store fit initial values and ranges:
  std::map<TString,TString> m_paramState;
  
  // Systematic: 
  TString m_listMRS;
  TString m_listMSS;
  
  // Yield data for plotting:
  std::map<int,TF1*> m_yieldFunc;
  std::map<int,TGraph*> m_yieldGraph;
  
  // List to track imported datasets:
  std::vector<std::pair<double,int> > m_massCatePairs;
  
  // True iff fits should be binned.
  bool m_binned;
  int m_nBinsPerGeV;
  
  // Fit parameter options:
  bool m_sameCBGAMean;

  // Plot options:
  bool m_useLogYAxis;
  bool m_doRatioPlot;
  double m_ratioMin;
  double m_ratioMax;
  TString m_currFunction;
};

#endif
