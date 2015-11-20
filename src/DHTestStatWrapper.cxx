////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHTestStatWrapper.cxx                                               //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/11/2015                                                          //
//                                                                            //
//  Includes a main method for using the DHTestStat.cxx class. This is useful //
//  for grid jobs.                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHTestStat.h"

int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 3) {
    std::cout << "\nUsage: " << argv[0] << " <config> <options>" << std::endl;
    exit(0);
  }
  
  TString configFile = argv[1];
  TString options = argv[2];
  
  // Load the analysis configuration file:
  Config *config = new Config(configFile);
  TString jobName = config->getStr("JobName");
  TString anaType = config->getStr("AnalysisType");
  
  // Define the input file, then make a local copy (for remote jobs):
  TString originFile = Form("%s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
			    (config->getStr("MasterOutput")).Data(), 
			    jobName.Data(), anaType.Data());
  
  TString copiedFile = "workspaceDH.root";
  system(Form("cp %s %s", originFile.Data(), copiedFile.Data()));
  
  // Load the RooWorkspace and ModelConfig:
  TFile inputFile(copiedFile, "read");
  RooWorkspace *workspace = (RooWorkspace*)inputFile.Get("combinedWS");
  
  // Instantiate test statistic class and calculate CL and p0:
  DHTestStat *ts = new DHTestStat(configFile, "new", workspace);
  ts->calculateNewCL();
  ts->calculateNewP0();
  if (ts->fitsAllConverged()) {
    std::cout << "DHTestStatWrapper: All OK!" << std::endl;
  }
  
  inputFile.Close();
  system(Form("rm %s", copiedFile.Data()));
  return 0;
}
