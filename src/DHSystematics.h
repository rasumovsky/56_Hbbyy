////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DHSystematics.h                                                           //
//  Class: DHSystematics.cxx                                                  //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 01/01/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DHSystematics_h
#define DHSystematics_h

// Package libraries:
#include "CommonHead.h"
//#include "CommonFunc.h"
#include "Config.h"

class DHSystematics {

 public:
  
  DHSystematics(TString newConfigFile, TString options);
  virtual ~DHSystematics() {};
  
  // Public accessor methods:
  double getInclSymmetricSystematic(TString sample, TString systematic);
  double getInclSystematic(TString sample, TString systematic, TString up_down);
  double getSymmetricSystematic(TString sample, TString systematic, 
				int category);
  double getSystematic(TString sample, TString systematic, int category,
		       TString up_down);
  std::vector<TString> listSampleNames();
  std::vector<TString> listSystematicNames();
  void printWorkspaceInput();
  
  // Public mutator methods:
  void addCategoryNames(std::vector<TString> cateNames);
  void addSystematic(TString sample, TString systematic, int category, 
		     TString up_down, double systematicValue);
  void addSystematicName(TString systematic);
  void groupSyst(TString sample, TString groupName,
		 std::vector<TString> sysComponents);
  void groupSystAllSamples(TString groupName, 
			   std::vector<TString> sysComponents);
  void loadSystematicsFile(TString fileName, TString sample);
  //void parameterizeSyst(std::map<TString,double> sampleToVar);
  void setSysToDefaults();
  void setConstrCenterTypeIncl(TString systematic, TString constraintType,
			       double centralValue, TString systematicType,
			       bool inclusive);
  
 private:
  
  // Private accessor methods:
  TString mKey(TString sample, TString systematic, int category,
	       TString up_down);
  void printer(TString statement, bool isFatal);

  
  // Private mutator methods:

  
  // Private member variables:
  TString m_outputDir;
  Config *m_config;
  
  int m_maxCategories;
  std::vector<TString> m_cateNames;
  std::vector<TString> m_sampleNames;
  std::map<TString,double> m_sysValues;
  std::vector<TString> m_sysNames;
  bool m_missingSys;
  
  // Systematic information:
  std::map<TString,TString> m_constraintType;
  std::map<TString,double> m_centralValue;
  std::map<TString,TString> m_sysType;
  std::map<TString,bool> m_inclusive;

};

#endif

