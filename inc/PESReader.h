////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PESReader.h                                                         //
//  Class: PESReader.cxx                                                      //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 15/04/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef PESReader_h
#define PESReader_h

#include <cmath>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"

class PESReader {
 
 public:
  
  PESReader(TString inputFileName, int nCategories);
  ~PESReader();
  
  int sourceNameToIndex(TString name);
  double getValue(TString name, int cateIndex);
  int getSign(TString name, int cateIndex);
  int getNumberOfSources();
  TString getNameOfSource(int sourceIndex);
  std::vector<TString> listSources();
  
 private:
  
  int nPESParams; 
  
  // store the PES values [#systematics][#categories]
  double valuesPES[100][20];
  
  // store the PES names:
  std::vector<TString> nameListPES;
  
};

#endif
