////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHWorkspaceWrapper.cxx                                              //
//                                                                            //
//  Creator: Andrew Hard,                                                     //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/11/2015                                                          //
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
  
  // Load job configuration and options:
  TString configFile = argv[1];
  TString options = argv[2];
  
  // Instantiate the workspace class:
  DHWorkspace *dmw = new DHWorkspace(configFile, options);
  if (dmw->fitsAllConverged()) {
    std::cout << "DMWorkspaceWrapper: All OK!" << std::endl;
  }
  return 0;
}
