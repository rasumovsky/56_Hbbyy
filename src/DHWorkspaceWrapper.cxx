////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHWorkspaceWrapper.cxx                                              //
//                                                                            //
//  Creator: Andrew Hard,                                                     //
//  Email: ahard@cern.ch                                                      //
//  Date: 21/07/2015                                                          //
//                                                                            //
//  This class builds the workspace for the di-Higgs analysis fits.           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHWorkspace.h"

int main(int argc, char **argv) {

  // Check that arguments are provided.
  if (argc < 4) { 
    std::cout << "\nUsage: " << argv[0] 
	      << " <configFile> <analysisType> <options>" << std::endl;
    exit(0);
  }
  
  TString configFile = argv[1];
  TString analysisType = argv[2];
  TString options = argv[3];
  
  DHWorkspace *dmw = new DHWorkspace(configFile, analysisType, options);
  if (dmw->fitsAllConverged()) {
    std::cout << "DMWorkspaceWrapper: All OK!" << std::endl;
  }
  return 0;
}
