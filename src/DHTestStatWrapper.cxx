////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHTestStatWrapper.cxx                                               //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 06/07/2015                                                          //
//                                                                            //
//  Includes a main method for using the DHTestStat.cxx class. This is useful //
//  for grid jobs.                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHTestStat.h"

int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 5) {
    std::cout << "\nUsage: " << argv[0]
	      << " <jobName> <DHSignal> <cateScheme> <options>" << std::endl;
    exit(0);
  }
  
  TString jobName = argv[1];
  TString DHSignal = argv[2];
  TString cateScheme = argv[3];
  TString options = argv[4];
  
  // Define the input file, then make a local copy (for remote jobs):
  TString originFile = Form("%s/%s/workspaces/rootfiles/workspaceDH.root",
			    masterOutput.Data(), jobName.Data());
  TString copiedFile = "workspaceDH.root";
  system(Form("cp %s %s", originFile.Data(), copiedFile.Data()));
  
  // Load the RooWorkspace and ModelConfig:
  TFile inputFile(copiedFile, "read");
  RooWorkspace *workspace = (RooWorkspace*)inputFile.Get("combinedWS");
  
  DHTestStat *ts = new DMTestStat(jobName, DHSignal, cateScheme, "new",
				  workspace);
  ts->calculateNewCL();
  ts->calculateNewP0();
  if (ts->fitsAllConverged()) {
    std::cout << "DMTestStatWrapper: All OK!" << std::endl;
  }
  
  inputFile.Close();
  system(Form("rm %s", copiedFile.Data()));
  return 0;
}
