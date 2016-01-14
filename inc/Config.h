////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  A brief class that handles reading input from a configuration text file   //
//  based on root's TEnv format, i.e. key-value pairs, and shamelessly stolen //
//  from the HGamAnalysisTools package.                                       //
//                                                                            //
//  Author: Dag Gillberg                                                      //
//  Appropriator: Andrew Hard                                                 //
//  Email: ahard@cern.ch                                                      //
//  Date: 03/08/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef Config_h
#define Config_h

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include "TString.h"
#include "TEnv.h"
#include "THashList.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"

class Config {

 public:
  
  Config(TString fileName);
  virtual ~Config() {};
  
  Config &operator = (const Config &rhs); // assignment operator
  
  
  TString getStr(TString key, bool expand=true);
  TString getStr(TString key, TString dflt);
  std::vector<TString>  getStrV(TString key);
  int getInt(TString key);
  int getInt(TString key, int dflt);
  bool getBool(TString key);
  bool getBool(TString key, bool dflt);
  double getNum(TString key);
  double getNum(TString key, double dflt);
  std::vector<double> getNumV(TString key);
  
  bool isDefined(TString key);
  void addFile(TString fileName);
  void setValue(TString key, TString value);
  
  //! \brief accessor to the TEnv database
  inline const TEnv* getDB() { return &m_env; }
  void printDB();
  std::vector<TString> vectorize(TString str, TString sep);
  std::vector<double> vectorizeNum(TString str, TString sep);
  
  bool fileExist(TString fn);
  
 private:
  
  TEnv m_env; // TEnv objects holding the settings database
  
  void ensureDefined(TString key);
  inline void copyTable(const TEnv &env);
  
};



#endif
