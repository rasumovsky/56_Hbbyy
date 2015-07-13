////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PESReader.cxx                                                       //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 15/04/2015                                                          //
//                                                                            //
//  This code reads the PES table and provides the values for the workspace.  //
//                                                                            //
//  WARNING! The category indexing is a bit unclear.                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "PESReader.h"

/**
   -----------------------------------------------------------------------------
   Initialize the PESReader class.
   @param inputFileName - the file location containing PES values.
   @param nCategories - the number of analysis categories.
   @returns - void.
*/
PESReader::PESReader(TString inputFileName, int nCategories) {
  
  std::cout << "\nInitializing the PESReader class" << std::endl;
  
  nPESParams = 0;
  nameListPES.clear();
  
  TString tempSourceName = "";
  double tempValues[20] = {0.0};
  
  ifstream sysFile(inputFileName);
  while( !sysFile.eof() )
  {
    sysFile >> tempSourceName;
    for (int i_c = 0; i_c <= nCategories; i_c++) {
      sysFile >> tempValues[i_c];
      valuesPES[nPESParams][i_c] = tempValues[i_c];
    }
    nameListPES.push_back(tempSourceName);
    nPESParams++;
  }
  sysFile.close();
  std::cout << "Finished loading data for the PESReader class." << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Finds the index associated with a particular PES source.
   @param name - the name of the PES.
   @returns - the index of the PES.
*/
int PESReader::sourceNameToIndex(TString name) {
  for (int i_n = 0; i_n < (int)nameListPES.size(); i_n++) {
    TString current_name = nameListPES[i_n];
    if (name.Contains(current_name)) {
      return i_n;
    }
  }
  return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the value of a particular PES in a given category.
   @param name - the name of the PES.
   @param cateIndex - the index of the category.
   @returns - the value of the PES in the category.
*/
double PESReader::getValue(TString name, int cateIndex) {
  return fabs(valuesPES[sourceNameToIndex(name)][cateIndex]);
}

/**
   -----------------------------------------------------------------------------
   Returns the sign of the PES.
   @param name - the name of the PES.
   @param cateIndex - the index of the category.
   @returns - the sign of the PES in the category.
*/
int PESReader::getSign(TString name, int cateIndex) {
  double value = valuesPES[sourceNameToIndex(name)][cateIndex];
  // make it absolute value (sign comes from GetSign():
  double sign = value / fabs(value);
  int result = (sign >= 0) ? 1 : -1;
  return result;
}

/**
   -----------------------------------------------------------------------------
   Returns the number of PES sources that have been defined.
   @returns - the number of PES sources.
*/
int PESReader::getNumberOfSources() {
  return nPESParams;
}

/**
   -----------------------------------------------------------------------------
   Retrieve the name of a source from its index.
   @param sourceIndex - the index of the PES.
   @returns - the name of the PES.
*/
TString PESReader::getNameOfSource(int indexPES) {
  return nameListPES[indexPES];
}

/**
   -----------------------------------------------------------------------------
   Retrieve the vector of PES systematics names.
*/
std::vector<TString> PESReader::listSources() {
  return nameListPES;
}
