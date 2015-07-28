////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHMaster.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/06/2015                                                          //
//                                                                            //
//  This program is useful as an interface to the H->diphoton + DM analysis   //
//  tools. It centralizes the commands for creating inputs, plots, workspaces,//
//  and statistical results. Some of the commands will rely on accessing      //
//  classes (mass points, signal parameterization), while others will use     //
//  system commands to submit jobs to various clusters.                       //
//                                                                            //
//  MasterOption - Note: Each can be followed by the suffix "New"             //
//    - MassPoints                                                            //
//    - SigParam                                                              //
//    - Workspace                                                             //
//    - ResubmitWorkspace                                                     //
//    - TossPseudoExp                                                         //
//    - PlotPseudoExp                                                         //
//    - TestStat                                                              //
//    - ResubmitTestStat                                                      //
//    - MuLimit                                                               //
//                                                                            //
//  Need to rethink the SigParam handling of the RooDataSet. Maybe we         //
//  should just hand it a RooDataSet?                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHMaster.h"


/**
   -----------------------------------------------------------------------------
   Submits the workspace jobs to the lxbatch server. 
   @param exeJobName - the job name.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
   @param exeCateScheme - the categorization scheme for the executable.
*/
void submitWSViaBsub(TString exeJobName, TString exeOption, TString exeSignal,
		     TString exeCateScheme) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_Workspace", DHAnalysis::clusterFileLocation.Data(),
		     exeJobName.Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  if (isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", DHAnalysis::exeWorkspace.Data()));
    system(Form("chmod +x %s/%s", DHAnalysis::packageLocation.Data(), 
		DHAnalysis::jobScriptWorkspace.Data()));
    system(Form("cp -f %s/%s %s/jobFileWorkspace.sh", 
		DHAnalysis::packageLocation.Data(), 
		DHAnalysis::jobScriptWorkspace.Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(), exeJobName.Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), exeJobName.Data(), 
			     exeSignal.Data());
  
  // Define the arguments for the job script:
  TString nameJobScript = Form("%s/jobFileWorkspace.sh %s %s %s %s %s %s",
			       exe.Data(), exeJobName.Data(), inputFile.Data(),
			       exeOption.Data(), 
			       DHAnalysis::exeWorkspace.Data(),
			       exeSignal.Data(), exeCateScheme.Data());
  // Submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJobScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   Submits the test statistics jobs to the lxbatch server. 
   @param exeJobName - the job name.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
   @param exeCateScheme - the categorization scheme for the executable.
*/
void submitTSViaBsub(TString exeJobName, TString exeOption, TString exeSignal,
		     TString exeCateScheme) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_TestStat", DHAnalysis::clusterFileLocation.Data(),
		     exeJobName.Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  if (isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", DHAnalysis::exeTestStat.Data()));
    system(Form("chmod +x %s/%s", DHAnalysis::packageLocation.Data(), 
		DHAnalysis::jobScriptTestStat.Data()));
    system(Form("chmod +x %s/%s/Workspace/rootfiles/workspace_%s.root",
		masterOutput.Data(), exeJobName.Data(), exeSignal.Data()));
    system(Form("cp -f %s/%s %s/jobFileTestStat.sh",
		DHAnalysis::packageLocation.Data(), 
		DHAnalysis::jobScriptTestStat.Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(),
			     exeJobName.Data(), exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), exeJobName.Data(),
			     exeSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJobScript = Form("%s/jobFileTestStat.sh %s %s %s %s %s %s", 
			       exe.Data(), exeJobName.Data(), inputFile.Data(),
			       DHAnalysis::exeTestStat.Data(), exeSignal.Data(),
			       exeCateScheme.Data(), exeOption.Data());
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(),
	      nameErrFile.Data(), nameJobScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   Submits the mu limit jobs to the lxbatch server. 
   @param exeJobName - the job name.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
*/
void SubmitMuLimitViaBsub(TString exeJobName, TString exeOption,
			  TString exeSignal) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_MuLimit", DHAnalysis::clusterFileLocation.Data(),
		     exeJobName.Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  if (isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", DHAnalysis::exeMuLimit.Data()));
    system(Form("chmod +x %s", DHAnalysis::jobScriptMuLimit.Data()));
    system(Form("chmod +x %s/%s/Workspace/rootfiles/workspace_%s.root", 
		masterOutput.Data(), exeJobName.Data(), exeSignal.Data()));
    system(Form("cp -f %s/%s %s/jobFileMuLimit.sh", packageLocation.Data(), 
		DHAnalysis::jobScriptMuLimit.Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(), exeJobName.Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), exeJobName.Data(),
			     exeSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJobScript = Form("%s/jobFileMuLimit.sh %s %s %s %s %s", 
			       exe.Data(), exeJobName.Data(), inputFile.Data(),
			       DHAnalysis::exeMuLimit.Data(), exeSignal.Data(),
			       exeOption.Data());
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJobScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   Submits the mu limit jobs to the lxbatch server. 
   @param exeJobName - the job name.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
   @param exeCateScheme - the categorization scheme for the executable.
   @param int exeSeed - the seed for the randomized dataset generation.
   @param int exeToysPerJob - the number of toy datasets to create per job.
*/
void submitPEViaBsub(TString exeJobName, TString exeOption, TString exeSignal,
		     TString exeCateScheme, int exeSeed, int exeToysPerJob) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_PseudoExp", DHAnalysis::clusterFileLocation.Data(),
		     exeJobName.Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  if (isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", DHAnalysis::exePseudoExp.Data()));
    system(Form("chmod +x %s", DHAnalysis::jobScriptPseudoExp.Data()));
    system(Form("chmod +x %s/%s/Workspace/rootfiles/workspace_%s.root", 
		masterOutput.Data(), exeJobName.Data(), exeSignal.Data()));
    system(Form("cp -f %s/%s %s/jobFilePseudoExp.sh", packageLocation.Data(), 
		DHAnalysis::jobScriptPseudoExp.Data(), exe.Data()));
    system(Form("chmod +x %s/jobFilePseudoExp.sh", exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s_%d.out", dir.Data(),
			     exeJobName.Data(), exeSignal.Data(), exeSeed);
  TString nameErrFile = Form("%s/err/%s_%s_%d.err", dir.Data(),
			     exeJobName.Data(), exeSignal.Data(), exeSeed);
 
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFilePseudoExp.sh %s %s %s %s %s %s %d %d", 
			     exe.Data(), exeJobName.Data(), inputFile.Data(),
			     DHAnalysis::exePseudoExp.Data(), exeSignal.Data(),
			     exeCateScheme.Data(), exeOption.Data(),
			     exeSeed, exeToysPerJob);
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   This is the main DHMaster method:
*/
int main (int argc, char **argv) {
  // Check arguments:
  if (argc < 4) {
    printf("\nUsage: %s <jobName> <option> <cateScheme>\n\n", argv[0]);
    exit(0);
  }
  
  // The job name and options (which analysis steps to perform):
  TString masterJobName = argv[1];
  TString masterOption = argv[2];
  TString masterCateScheme = argv[3];
  
  // Submit jobs to bsub or grid, etc.:
  bool runInParallel = false;
  isFirstJob = true;
  
  // Options for each step:
  TString workspaceOptions = "New_nosys";//"FromFile_nosys";
  TString pseudoExpOptions = "FixMu";
  TString toyPlotOptions   = "null";
  TString testStatOptions  = "New";//"FromFile";
  TString muLimitOptions   = "null";
  
  //--------------------------------------//
  // Step 4.1: Create the workspace for fitting:
  if (masterOption.Contains("Workspace") && 
      !masterOption.Contains("ResubmitWorkspace")) {
    std::cout << "DHMaster: Step 4.1 - Making the workspaces." << std::endl;
    
    int jobCounterWS = 0;
    for (int i_s = 0; i_s < DHAnalysis::nDHModes; i_s++) {
      TString currSignal = DHAnalysis::sigDHModes[i_s];
      if (runInParallel) {
	submitWSViaBsub(masterJobName, workspaceOptions, currSignal,
			masterCateScheme);
	jobCounterWS++;
	isFirstJob = false;
      }
      else {
	DHWorkspace *ws = new DHWorkspace(masterJobName, currSignal,
					  masterCateScheme, workspaceOptions);
	if (dhw->fitsAllConverged()) {
	  jobCounterWS++;
	}
	else {
	  std::cout << "DHMaster: Problem with workspace fit!" << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Submitted/completed " << jobCounterWS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 4.2: Resubmit any failed workspace jobs:
  if (masterOption.Contains("ResubmitWorkspace")) {
    std::cout << "DHMaster: Step 4.2 - Resubmit failed workspace." << std::endl;
    
    int jobCounterWS = 0;
    // Get the points to resubmit:
    CheckJobs *cj = new CheckJobs(masterJobName);
    vector<TString> resubmitSignals = cj->getResubmitList("Workspace");
    cj->printResubmitList("Workspace");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitSignals.size()
	      << " workspace jobs." << std::endl;
    for (int i_s = 0; i_s < (int)resubmitSignals.size(); i_s++) {
      TString currSignal = resubmitSignals[i_s];
      
      if (runInParallel) {
	submitWSViaBsub(exeWorkspace, masterJobName, workspaceOptions, 
			currSignal);
	jobCounterWS++;
	isFirstJob = false;
      }
      else {
	DHWorkspace *dhw = new DHWorkspace(masterJobName, currSignal,
					   masterCateScheme, workspaceOptions);
	if (dhw->fitsAllConverged()) {
	  jobCounterWS++;
	}
	else {
	  std::cout << "DHMaster: Problem with workspace fit!" << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Resubmitted " << jobCounterWS << " jobs" << std::endl;
  }
  /*
  //--------------------------------------//
  // Step 5.1: Create pseudoexperiment ensemble:
  TString currSignal = sigModes[2];
  if (masterOption.Contains("TossPseudoExp")) {
    cout << "DHMaster: Step 5.1 - Creating pseudoexperiments for signal "
	 << currSignal << std::endl;
    
    int toySeed = 1987;
    int nToysTotal = 10000;
    int nToysPerJob = 50;
    int increment = nToysPerJob;
    int highestSeed = toySeed + nToysTotal;
    
    for (int i_s = toySeed; i_s < highestSeed; i_s += increment) {
      submitPEViaBsub(masterJobName, pseudoExpOptions, currSignal,
		      masterCateScheme, i_s, nToysPerJob);
      isFirstJob = false;
    }
    std::cout << "DHMaster: Submitted " << (int)(nToysTotal/nToysPerJob) 
	      << " total pseudo-experiments." << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.2: Plot pseudo-experiment ensemble results:
  if (masterOption.Contains("PlotPseudoExp")) {
    std::cout << "DHMaster: Step 5.2 - Plot pseudoexperiment results for "
	      << currSignal << std::endl;    
    ToyAnalysis *ta = new ToyAnalysis(masterJobName, currSignal,
				      masterCateScheme, toyPlotOptions);
  }
  
  //--------------------------------------//
  // Step 6.1: Calculate the test statistics:
  if (masterOption.Contains("TestStat") && 
      !masterOption.Contains("ResubmitTestStat")) {
    std::cout << "DHMaster: Step 6.1 - Calculating CL and p0." << std::endl;

    int jobCounterTS = 0;
    for (int i_s = 0; i_s < DHAnalysis::nDHModes; i_s++) {
      TString currSignal = DHAnalysis::sigDHModes[i_s];
      
      if (runInParallel) {
	submitTSViaBsub(masterJobName, testStatOptions, currSignal,
			masterCateScheme);
	jobCounterTS++;
	isFirstJob = false;
      }
      else {
	TestStat *ts = new TestStat(masterJobName, currSignal, 
				    masterCateScheme, testStatOptions, NULL);
	ts->calculateNewCL();
	ts->calculateNewP0();
	if (ts->fitsAllConverged()) jobCounterTS++;
	else {
	  std::cout << "DHMaster: Problem with test-stat fit!" << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Submitted/completed " << jobCounterTS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 6.2: Resubmit any failed test statistics jobs:
  if (masterOption.Contains("ResubmitTestStat")) {
    std::cout << "DHMaster: Step 6.2 - Resubmit failed test stat." << std::endl;
    
    int jobCounterTS = 0;
    // Get the points to resubmit:
    CheckJobs *cj = new CheckJobs(masterJobName);
    vector<TString> resubmitSignals = cj->getResubmitList("TestStat");
    cj->printResubmitList("TestStat");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitSignals.size()
	      << " workspace jobs." << std::endl;
    for (int i_s = 0; i_s < (int)resubmitSignals.size(); i_s++) {
      TString currSignal = resubmitSignals[i_s];
      
      if (runInParallel) {
	submitTSViaBsub(masterJobName, testStatOptions, currSignal, 
			masterCateScheme);
      	jobCounterTS++;
	isFirstJob = false;
      }
      else {
	TestStat *ts = new TestStat(masterJobName, currSignal,
				    masterCateScheme, testStatOptions, NULL);
	ts->calculateNewCL();
	ts->calculateNewP0();
	if (ts->fitsAllConverged()) jobCounterTS++;
	else {
	  std::cout << "DHMaster: Problem with test-stat fit!" << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Resubmitted " << jobCounterTS << " jobs" << std::endl;
  }
  */
  //--------------------------------------//
  // Step 7.1: Calculate the limits on the dark matter signal strength.
  if (masterOption.Contains("MuLimit") &&
      !masterOption.Contains("ResubmitMuLimit")) {
    std::cout << "DHMaster: Step 7.1 - Calculate 95%CL mu value." << std::endl;

    int jobCounterML = 0;
    for (int i_s = 0; i_s < DHAnalysis::nDHModes; i_s++) {
      TString currSignal = DHAnalysis::sigDHModes[i_s];
      
      if (runInParallel) {
	submitMLViaBsub(masterJobName, muLimitOptions, currSignal);
	isFirstJob = false;
      }
      else {
	TString muCommand = Form(".%s/bin/%s %s %s %s", packageLocation.Data(), 
				 exeMuLimit.Data(), masterJobName.Data(),
				 currSignal.Data(), muLimitOptions.Data());
	std::cout << "Executing following system command: \n\t"
		  << muCommand << std::endl;
	system(muCommand);
      }
      jobCounterML++;
    }
    std::cout << "Submitted/completed " << jobCounterML << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 7.2: Resubmit any failed mu limit jobs:
  if (masterOption.Contains("ResubmitMuLimit")) {
    std::cout << "DHMaster: Step 7.2 - Resubmit failed mu limits." << std::endl;
    
    int jobCounterML = 0;
    // Get the points to resubmit:
    CheckJobs *cj = new CheckJobs(masterJobName);
    vector<TString> resubmitSignals = cj->getResubmitList("MuLimit");
    cj->printResubmitList("MuLimit");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitSignals.size()
	      << " mu limit jobs." << std::endl;
    for (int i_s = 0; i_s < (int)resubmitSignals.size(); i_s++) {
      TString currSignal = resubmitSignals[i_s];
      
      if (runInParallel) {
	submitMLViaBsub(masterJobName, muLimitOptions, currSignal);
	isFirstJob = false;
      }
      else {
	system(Form(".%s/bin/%s %s %s %s", packageLocation.Data(), 
		    exeMuLimit.Data(), masterJobName.Data(),
		    currSignal.Data(), muLimitOptions.Data()));
      }
      jobCounterML++;
    }
    std::cout << "Resubmitted " << jobCounterML << " jobs" << std::endl;
  }
  
  return 0;
}
