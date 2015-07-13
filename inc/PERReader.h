////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PERReader.h                                                         //
//  Class: PERReader.cxx                                                      //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 15/04/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef PERReader_h
#define PERReader_h

#include <cmath>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
class PERReader {

 public:
  
  PERReader(TString inputFileName, int nCategories);
  ~PERReader();
  
  int sourceNameToIndex(TString name); 
  double getValue(TString name, int cateIndex);
  int getSign(TString name, int cateIndex);
  int getNumberOfSources();
  TString getNameOfSource(int indexPER);
  std::vector<TString> listSources();
  
 private:
  
  int nPERParams;
  
  // store the RES values [#systematics][#categories]
  double valuesPER[100][20];
  
  // store the RES names:
  std::vector<TString> nameListPER;
  
};

#endif
