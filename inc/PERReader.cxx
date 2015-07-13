////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PERReader.cxx                                                       //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 15/04/2015                                                          //
//                                                                            //
//  This code reads the resolution systematic uncertainty table and provides  //
//  the values for the workspace code.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "PERReader.h"

/**
   -----------------------------------------------------------------------------
   Initialize the PERReader class.
   @param inputFileName - the file location containing PER sys values.
   @param nCategories - the number of analysis categories.
   @returns - void.
*/
PERReader::PERReader(TString inputFileName, int nCategories) {
  
  std::cout << "\nInitializing the PERReader class" << std::endl;
  
  nPERParams = 0;
  nameListPER.clear();
  
  TString tempSourceName = "";
  double tempValues[20] = {0.0};
  
  ifstream sysFile(inputFileName);
  while (!sysFile.eof()) {
    
    sysFile >> tempSourceName;
    for (int i_c = 0; i_c <= nCategories; i_c++) {
      sysFile >> tempValues[i_c];
      valuesPER[nPERParams][i_c] = tempValues[i_c];
    }
    nameListPER.push_back(tempSourceName);
    nPERParams++;
  }
  sysFile.close();
  std::cout << "Finished loading data for the PERReader class." << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Finds the index associated with a particular PER sys source.
   @param name - the name of the PER sys.
   @returns - the index of the PER sys.
*/
int PERReader::sourceNameToIndex(TString name) {
  for (int i_n = 0; i_n < (int)nameListPER.size(); i_n++) {
    if (name.Contains(nameListPER[i_n])) {
      return i_n;
    }
  }
  return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the value of a particular PER sys. in a given category.
   @param name - the name of the PER sys.
   @param cateIndex - the index of the category.
   @returns - the value of the PER sys. in the category.
*/
double PERReader::getValue(TString name, int cateIndex) {
  return fabs(valuesPER[sourceNameToIndex(name)][cateIndex]);
}

/**
   -----------------------------------------------------------------------------
   Returns the sign of the PER sys.
   @param name - the name of the PER sys.
   @param cateIndex - the index of the category.
   @returns - the sign of the PER sys. in the category.
*/
int PERReader::getSign(TString name, int cateIndex) {
  double value = valuesPER[sourceNameToIndex(name)][cateIndex];
  // make it absolute value (sign comes from GetSign():
  double sign = value / fabs(value);
  int result = ( sign >= 0 ) ? 1 : -1;
  return result;
}

/**
   -----------------------------------------------------------------------------
   Returns the number of PER sys. sources that have been defined.
   @returns - the number of PER sys. sources.
*/
int PERReader::getNumberOfSources() {
  return nPERParams;
}

/**
   -----------------------------------------------------------------------------
   Retrieve the name of a source from its index.
   @param sourceIndex - the index of the PER sys.
   @returns - the name of the PER sys.
*/
TString PERReader::getNameOfSource(int indexPER) {
  return nameListPER[indexPER];
}

/**
   -----------------------------------------------------------------------------
   Retrieve the vector of PES systematics names.
*/
std::vector<TString> PERReader::listSources() {
  return nameListPER;
}
