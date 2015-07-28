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
	      << " <jobName> <cateScheme> <options>" << std::endl;
    exit(0);
  }
  
  TString jobName = argv[1];
  TString DHSignal = argv[2];
  TString cateScheme = argv[3];
  TString options = argv[4];
  
  DHTestStat *dhts = new DHTestStat(jobName, DHSignal, cateScheme,
				    "FromFile", m_combinedWS);

  DHWorkspace *dmw = new DHWorkspace(jobName, cateScheme, options);
  if (dmw->fitsAllConverged()) {
    std::cout << "DMWorkspaceWrapper: All OK!" << std::endl;
  }
  return 0;
}
