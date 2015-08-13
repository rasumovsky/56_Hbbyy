////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHMaster.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/08/2015                                                          //
//                                                                            //
//  This program is useful as an interface to the di-Higgs (bb+yy) analysis   //
//  tools. It centralizes the commands for creating inputs, plots, workspaces,//
//  and statistical results. Some of the commands will rely on accessing      //
//  classes (mass points, signal parameterization), while others will use     //
//  system commands to submit jobs to various clusters.                       //
//                                                                            //
//  MasterOption - Note: Each can be followed by the suffix "New"             //
//    - MassPoints                                                            //
//    - SigParam                                                              //
//    - Workspace                                                             //
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
   Submits the test statistics jobs to the lxbatch server. 
   @param exeConfigFile - the config file.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
 */
void submitTSViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal) {  
  // Make directories for job info:
  TString dir = Form("%s/%s_DHTestStat", 
		     (m_config->getStr("clusterFileLocation")).Data(),
		     (m_config->getStr("jobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  TString exeAna = DHAnalysis::getAnalysisType(m_config, exeSignal);
  
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exeTestStat")).Data()));
    system(Form("chmod +x %s/%s", (m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptTestStat")).Data()));
    system(Form("chmod +x %s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
		(m_config->getStr("masterOutput")).Data(), 
		(m_config->getStr("jobName")).Data(), exeAna.Data()));
    system(Form("cp -f %s/%s %s/jobFileTestStat.sh",
		(m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptTestStat")).Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(),
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), 
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFileTestStat.sh %s %s %s %s %s %s", 
			     exe.Data(), (m_config->getStr("jobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exeTestStat")).Data(),
			     exeSignal.Data(), exeOption.Data());
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(),
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   Submits the mu limit jobs to the lxbatch server. 
   @param exeConfigFile - the config file.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
*/
void SubmitMuLimitViaBsub(TString exeConfigFile, TString exeOption,
			  TString exeSignal) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_DHMuLimit", 
		     (m_config->getStr("clusterFileLocation")).Data(),
		     (m_config->getStr("jobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  TString exeAna = DHAnalysis::getAnalysisType(m_config, exeSignal);
  
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exeMuLimit")).Data()));
    system(Form("chmod +x %s", (m_config->getStr("jobScriptMuLimit")).Data()));
    system(Form("chmod +x %s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root", 
		(m_config->getStr("masterOutput")).Data(), 
		(m_config->getStr("jobName")).Data(), exeAna.Data()));
    system(Form("cp -f %s/%s %s/jobFileMuLimit.sh", 
		(m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptMuLimit")).Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(), 
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), 
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFileMuLimit.sh %s %s %s %s %s %s", 
			     exe.Data(), (m_config->getStr("jobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exeMuLimit")).Data(),
			     exeSignal.Data(), exeOption.Data());
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   Submits the mu limit jobs to the lxbatch server. 
   @param exeConfigFile - the config file.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
   @param int exeSeed - the seed for the randomized dataset generation.
   @param int exeToysPerJob - the number of toy datasets to create per job.
*/
void submitPEViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal, int exeSeed, int exeToysPerJob) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_PseudoExp",
		     (m_config->getStr("clusterFileLocation")).Data(),
		     (m_config->getStr("jobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  TString exeAna = DHAnalysis::getAnalysisType(m_config, exeSignal);
  
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exePseudoExp")).Data()));
    system(Form("chmod +x %s",(m_config->getStr("jobScriptPseudoExp")).Data()));
    system(Form("chmod +x %s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root", 
		(m_config->getStr("masterOutput")).Data(), 
		(m_config->getStr("jobName")).Data(), exeAna.Data()));
    system(Form("cp -f %s/%s %s/jobFilePseudoExp.sh", 
		(m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptPseudoExp")).Data(), exe.Data()));
    system(Form("chmod +x %s/jobFilePseudoExp.sh", exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s_%d.out", dir.Data(),
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data(), exeSeed);
  TString nameErrFile = Form("%s/err/%s_%s_%d.err", dir.Data(),
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data(), exeSeed);
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFilePseudoExp.sh %s %s %s %s %s %s %d %d", 
			     exe.Data(), (m_config->getStr("jobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exePseudoExp")).Data(),
			     exeSignal.Data(), exeOption.Data(), exeSeed,
			     exeToysPerJob);
  
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
  if (argc < 3) {
    printf("\nUsage: %s <jobName> <option> <cateScheme>\n\n", argv[0]);
    exit(0);
  }
  
  // The job name and options (which analysis steps to perform):
  TString masterOption = argv[1];
  TString configFileName = argv[2];
  
  // Submit jobs to bsub or grid, etc.:
  bool runInParallel = false;
  m_isFirstJob = true;
  
  // Load the config class and file:
  std::cout << "DHMaster: Loading the global config file." << std::endl;
  m_config = new Config(configFileName);
  m_config->printDB();
  TString fullConfigPath = Form("%s/%s",
				(m_config->getStr("packageLocation")).Data(),
				configFileName.Data());
  
  // Options for each step:
  TString workspaceOptions = m_config->getStr("workspaceOptions");
  TString pseudoExpOptions = m_config->getStr("pseudoExpOptions");
  TString toyPlotOptions   = m_config->getStr("toyPlotOptions");
  TString testStatOptions  = m_config->getStr("testStatOptions");
  TString muLimitOptions   = m_config->getStr("muLimitOptions");
  
  //--------------------------------------//
  // Step 4.1: Create the workspace for fitting:
  if (masterOption.Contains("Workspace") && 
      !masterOption.Contains("ResubmitWorkspace")) {
    std::cout << "DHMaster: Step 4.1 - Making the workspaces." << std::endl;
    std::vector<TString> analysisTypes = m_config->getStrV("analysisTypes");
    for (int i_a = 0; i_a < (int)analysisTypes.size(); i_a++) {
      
      if (workspaceOptions.Contains(Form("No%s",analysisTypes[i_a].Data()))) {
	continue;
      }
      std::cout << "DHMaster: Creating ws for " << analysisTypes[i_a] 
		<< " analysis." << std::endl;
      
      DHWorkspace *dhws = new DHWorkspace(configFileName, analysisTypes[i_a],
					  workspaceOptions);
      if (!dhws->fitsAllConverged()) {
	std::cout << "DHMaster: Problem with workspace fit!" << std::endl;
	exit(0);
      }
    }
    std::cout << "DHMaster: Constructed " << analysisTypes.size() 
	      << " workspaces." << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.1: Create pseudoexperiment ensemble:
  TString currToySignal = m_config->getStr("exampleSignal");
  if (masterOption.Contains("TossPseudoExp")) {
    cout << "DHMaster: Step 5.1 - Creating pseudoexperiments for signal "
	 << currToySignal << std::endl;
    
    int toySeed = m_config->getInt("toySeed");
    int nToysTotal = m_config->getInt("nToysTotal");
    int nToysPerJob = m_config->getInt("nToysPerJob");
    int increment = nToysPerJob;
    int highestSeed = toySeed + nToysTotal;
    
    for (int i_s = toySeed; i_s < highestSeed; i_s += increment) {
      submitPEViaBsub(fullConfigPath, pseudoExpOptions, currToySignal, i_s,
		      nToysPerJob);
      m_isFirstJob = false;
    }
    std::cout << "DHMaster: Submitted " << (int)(nToysTotal/nToysPerJob) 
	      << " total pseudo-experiments." << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.2: Plot pseudo-experiment ensemble results:
  if (masterOption.Contains("PlotPseudoExp")) {
    std::cout << "DHMaster: Step 5.2 - Plot pseudoexperiment results for "
	      << currToySignal << std::endl;
    DHToyAnalysis *dhta = new DHToyAnalysis(configFileName, currToySignal);
  }
  
  //--------------------------------------//
  // Step 6.1: Calculate the test statistics:
  if (masterOption.Contains("TestStat") && 
      !masterOption.Contains("ResubmitTestStat")) {
    std::cout << "DHMaster: Step 6.1 - Calculating CL and p0." << std::endl;

    int jobCounterTS = 0;
    std::vector<TString> sigDHModes = m_config->getStrV("sigDHModes");
    for (int i_s = 0; i_s < (int)sigDHModes.size(); i_s++) {
      TString currSignal = sigDHModes[i_s];
      
      if (runInParallel) {
	submitTSViaBsub(fullConfigPath, testStatOptions, currSignal);
	jobCounterTS++;
	m_isFirstJob = false;
      }
      else {
	DHTestStat *ts = new DHTestStat(configFileName, currSignal,
					testStatOptions, NULL);
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
    DHCheckJobs *dhc = new DHCheckJobs(configFileName);
    vector<TString> resubmitSignals = dhc->getResubmitList("DHTestStat");
    dhc->printResubmitList("DHTestStat");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitSignals.size()
	      << " test statistic jobs." << std::endl;
    for (int i_s = 0; i_s < (int)resubmitSignals.size(); i_s++) {
      TString currSignal = resubmitSignals[i_s];
      
      if (runInParallel) {
	submitTSViaBsub(fullConfigPath, testStatOptions, currSignal);
      	jobCounterTS++;
	m_isFirstJob = false;
      }
      else {
	DHTestStat *dhts = new DHTestStat(configFileName, currSignal,
					  testStatOptions, NULL);
	dhts->calculateNewCL();
	dhts->calculateNewP0();
	if (dhts->fitsAllConverged()) jobCounterTS++;
	else {
	  std::cout << "DHMaster: Problem with test-stat fit!" << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Resubmitted " << jobCounterTS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 7.1: Calculate the limits on the dark matter signal strength.
  if (masterOption.Contains("MuLimit") &&
      !masterOption.Contains("ResubmitMuLimit")) {
    std::cout << "DHMaster: Step 7.1 - Calculate 95%CL mu value." << std::endl;

    int jobCounterML = 0;
    std::vector<TString> sigDHModes = m_config->getStrV("sigDHModes");
    for (int i_s = 0; i_s < (int)sigDHModes.size(); i_s++) {
      TString currSignal = sigDHModes[i_s];
      
      if (runInParallel) {
	submitMLViaBsub(fullConfigPath, muLimitOptions, currSignal);
	m_isFirstJob = false;
      }
      else {
	TString muCommand = Form(".%s/bin/%s %s %s %s", 
				 (m_config->getStr("packageLocation")).Data(), 
				 (m_config->getStr("exeMuLimit")).Data(),
				 fullConfigPath.Data(), currSignal.Data(),
				 muLimitOptions.Data());
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
    DHCheckJobs *dhc = new DHCheckJobs(configFileName);
    vector<TString> resubmitSignals = dhc->getResubmitList("DHMuLimit");
    dhc->printResubmitList("DHMuLimit");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitSignals.size()
	      << " mu limit jobs." << std::endl;
    for (int i_s = 0; i_s < (int)resubmitSignals.size(); i_s++) {
      TString currSignal = resubmitSignals[i_s];
      
      if (runInParallel) {
	submitMLViaBsub(fullConfigPath, muLimitOptions, currSignal);
	m_isFirstJob = false;
      }
      else {
	system(Form(".%s/bin/%s %s %s %s", 
		    (m_config->getStr("packageLocation")).Data(), 
		    (m_config->getStr("exeMuLimit")).Data(),
		    fullConfigPath.Data(), currSignal.Data(),
		    muLimitOptions.Data()));
      }
      jobCounterML++;
    }
    std::cout << "Resubmitted " << jobCounterML << " jobs" << std::endl;
  }
  
  return 0;
}
