////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHMaster.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 01/01/2016                                                          //
//                                                                            //
//  This program is useful as an interface to the di-Higgs (bb+yy) analysis   //
//  tools. It centralizes the commands for creating inputs, plots, workspaces,//
//  and statistical results. Some of the commands will rely on accessing      //
//  classes (mass points, signal parameterization), while others will use     //
//  system commands to submit jobs to various clusters.                       //
//                                                                            //
//  MasterOption - Note: Each can be followed by the suffix "New"             //
//    - Systematics                                                           //
//    - Workspace                                                             //
//    - TossPseudoExp                                                         //
//    - PlotPseudoExp                                                         //
//    - TestStat                                                              //
//    - ResubmitTestStat                                                      //
//    - MuLimit                                                               //
//    - CLScanAnalysis                                                        //
//    - CLScanSubmitToys                                                      //
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
		     (m_config->getStr("ClusterFileLocation")).Data(),
		     (m_config->getStr("JobName")).Data());
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
    system(Form("chmod +x %s/%s", (m_config->getStr("PackageLocation")).Data(), 
		(m_config->getStr("jobScriptTestStat")).Data()));
    system(Form("chmod +x %s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
		(m_config->getStr("MasterOutput")).Data(), 
		(m_config->getStr("JobName")).Data(), exeAna.Data()));
    system(Form("cp -f %s/%s %s/jobFileTestStat.sh",
		(m_config->getStr("PackageLocation")).Data(), 
		(m_config->getStr("jobScriptTestStat")).Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(),
			     (m_config->getStr("JobName")).Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), 
			     (m_config->getStr("JobName")).Data(),
			     exeSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFileTestStat.sh %s %s %s %s %s %s", 
			     exe.Data(), (m_config->getStr("JobName")).Data(),
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
		     (m_config->getStr("ClusterFileLocation")).Data(),
		     (m_config->getStr("JobName")).Data());
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
		(m_config->getStr("MasterOutput")).Data(), 
		(m_config->getStr("JobName")).Data(), exeAna.Data()));
    system(Form("cp -f %s/%s %s/jobFileMuLimit.sh", 
		(m_config->getStr("PackageLocation")).Data(), 
		(m_config->getStr("jobScriptMuLimit")).Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(), 
			     (m_config->getStr("JobName")).Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), 
			     (m_config->getStr("JobName")).Data(),
			     exeSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFileMuLimit.sh %s %s %s %s %s %s", 
			     exe.Data(), (m_config->getStr("JobName")).Data(),
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
   @param int exeSeed - the seed for the randomized dataset generation.
   @param int exeToysPerJob - the number of toy datasets to create per job.
*/
void submitPEViaBsub(TString exeConfigFile, TString exeOption, int exeSeed, 
		     int exeToysPerJob) {
  //@param exeSignal - the signal to process in the executable.

  // Make directories for job info:
  TString dir = Form("%s/%s_PseudoExp",
		     (m_config->getStr("ClusterFileLocation")).Data(),
		     (m_config->getStr("JobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  TString exeAna = m_config->getStr("AnalysisType");
  
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exePseudoExp")).Data()));
    system(Form("chmod +x %s",(m_config->getStr("jobScriptPseudoExp")).Data()));
    system(Form("chmod +x %s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root", 
		(m_config->getStr("MasterOutput")).Data(), 
		(m_config->getStr("JobName")).Data(), exeAna.Data()));
    system(Form("cp -f %s/%s %s/jobFilePseudoExp.sh", 
		(m_config->getStr("PackageLocation")).Data(), 
		(m_config->getStr("jobScriptPseudoExp")).Data(), exe.Data()));
    system(Form("chmod +x %s/jobFilePseudoExp.sh", exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%d.out", dir.Data(),
			     (m_config->getStr("JobName")).Data(), exeSeed);
  TString nameErrFile = Form("%s/err/%s_%d.err", dir.Data(),
			     (m_config->getStr("JobName")).Data(), exeSeed);
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFilePseudoExp.sh %s %s %s %s %s %d %d", 
			     exe.Data(), (m_config->getStr("JobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exePseudoExp")).Data(),
			     exeOption.Data(), exeSeed, exeToysPerJob);
  
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
    printf("\nUsage: %s <option> <configFileName>\n\n", argv[0]);
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
  TString fullConfigPath
    = Form("%s/%s", (m_config->getStr("PackageLocation")).Data(),
	   configFileName.Data());
  
  //--------------------------------------//
  // Step 3: Calculate systematic uncertainties and output a file for .cfg
  if (masterOption.Contains("Systematics")) {
    DHSystematics *sys = new DHSystematics(configFileName, "");
    sys->addCategoryNames(m_config->getStrV("CateNames"));
    std::vector<TString> sysComponents = m_config->getStrV("SysComponents");
    for (int i_s = 0; i_s < (int)sysComponents.size(); i_s++) {
      TString fileName = Form("%s/%s_allSys.txt",
			      (m_config->getStr("SysDirectory")).Data(),
			      (sysComponents[i_s]).Data());
      sys->loadSystematicsFile(fileName, sysComponents[i_s]);
    }
    sys->setSysToDefaults();
    
    // B flavor grouping:
    std::vector<TString> vec_FT_EFF_Eigen_B; vec_FT_EFF_Eigen_B.clear(); 
    vec_FT_EFF_Eigen_B.push_back("FT_EFF_Eigen_B_0");
    vec_FT_EFF_Eigen_B.push_back("FT_EFF_Eigen_B_1");
    vec_FT_EFF_Eigen_B.push_back("FT_EFF_Eigen_B_2");
    vec_FT_EFF_Eigen_B.push_back("FT_EFF_Eigen_B_3");
    vec_FT_EFF_Eigen_B.push_back("FT_EFF_Eigen_B_4");
    vec_FT_EFF_Eigen_B.push_back("FT_EFF_Eigen_B_5");
    sys->groupSystAllSamples("FT_EFF_Eigen_B", vec_FT_EFF_Eigen_B);
    
    // C flavor grouping:
    std::vector<TString> vec_FT_EFF_Eigen_C; vec_FT_EFF_Eigen_C.clear(); 
    vec_FT_EFF_Eigen_C.push_back("FT_EFF_Eigen_C_0");
    vec_FT_EFF_Eigen_C.push_back("FT_EFF_Eigen_C_1");
    vec_FT_EFF_Eigen_C.push_back("FT_EFF_Eigen_C_2");
    vec_FT_EFF_Eigen_C.push_back("FT_EFF_Eigen_C_3");
    sys->groupSystAllSamples("FT_EFF_Eigen_C", vec_FT_EFF_Eigen_C);
    
    // Light flavor grouping:
    std::vector<TString> vec_FT_EFF_Eigen_Light;
    vec_FT_EFF_Eigen_Light.clear(); 
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_0");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_1");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_2");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_3");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_4");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_5");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_6");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_7");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_8");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_9");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_10");
    vec_FT_EFF_Eigen_Light.push_back("FT_EFF_Eigen_Light_11");
    sys->groupSystAllSamples("FT_EFF_Eigen_Light", vec_FT_EFF_Eigen_Light);
    
    // JET uncertainty group:
    std::vector<TString> vec_JET_GroupedNP; vec_JET_GroupedNP.clear(); 
    vec_JET_GroupedNP.push_back("JET_GroupedNP_1");
    vec_JET_GroupedNP.push_back("JET_GroupedNP_2");
    vec_JET_GroupedNP.push_back("JET_GroupedNP_3");
    sys->groupSystAllSamples("JET_GroupedNP", vec_JET_GroupedNP);
    
    // Set the types for the nuisance parameters:
    sys->setConstrCenterTypeIncl("FT_EFF_extrapolation_from_charm",
				 "gaus", 1.0, "yield", false);
    sys->setConstrCenterTypeIncl("FT_EFF_extrapolation",
				 "gaus", 1.0, "yield", false);
    sys->setConstrCenterTypeIncl("JET_JER_SINGLE_NP",
				 "gaus", 1.0, "yield", false);
    sys->setConstrCenterTypeIncl("PH_EFF_ID_Uncertainty",
				 "gaus", 1.0, "yield", true);
    sys->setConstrCenterTypeIncl("PH_EFF_TRKISO_Uncertainty",
				 "gaus", 1.0, "yield", true);
    sys->setConstrCenterTypeIncl("PRW_DATASF",
				 "gaus", 1.0, "yield", false);
    sys->setConstrCenterTypeIncl("FT_EFF_Eigen_B",
				 "gaus", 1.0, "yield", false);
    sys->setConstrCenterTypeIncl("FT_EFF_Eigen_C",
				 "gaus", 1.0, "yield", false);
    sys->setConstrCenterTypeIncl("FT_EFF_Eigen_Light",
				 "asym", 1.0, "yield", false);
    sys->setConstrCenterTypeIncl("JET_GroupedNP",
				 "asym", 1.0, "yield", false);
    
    // Print the systematic uncertainty input:
    sys->printWorkspaceInput();
    delete sys;
  }
  
  //--------------------------------------//
  // Step 4.1: Create the workspace for fitting:
  if (masterOption.Contains("Workspace") && 
      !masterOption.Contains("ResubmitWorkspace")) {
    std::cout << "DHMaster: Step 4.1 - Making the workspaces." << std::endl;
    DHWorkspace *dhws = new DHWorkspace(configFileName, 
					m_config->getStr("WorkspaceOptions"));
    if (!dhws->fitsAllConverged()) {
      std::cout << "DHMaster: Problem with workspace fit!" << std::endl;
      exit(0);
    }
    std::cout << "DHMaster: Successfully built workspace!" << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.1: Create pseudoexperiment ensemble:
  if (masterOption.Contains("TossPseudoExp")) {
    std::cout << "DHMaster: Step 5.1 - Create pseudoexperiments." << std::endl;
    
    int toySeed = m_config->getInt("toySeed");
    int nToysTotal = m_config->getInt("nToysTotal");
    int nToysPerJob = m_config->getInt("nToysPerJob");
    int increment = nToysPerJob;
    int highestSeed = toySeed + nToysTotal;
    
    for (int i_s = toySeed; i_s < highestSeed; i_s += increment) {
      submitPEViaBsub(fullConfigPath, m_config->getStr("PseudoExpOptions"),
		      i_s, nToysPerJob);
      m_isFirstJob = false;
    }
    std::cout << "DHMaster: Submitted " << (int)(nToysTotal/nToysPerJob) 
	      << " total pseudo-experiments." << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.2: Plot pseudo-experiment ensemble results:
  if (masterOption.Contains("PlotPseudoExp")) {
    std::cout << "DHMaster: Step 5.2 - Plot pseudoexperiment results."
	      << std::endl;
    DHToyAnalysis *dhta = new DHToyAnalysis(configFileName, "NONE");
    //DHToyAnalysis *dhta = new DHToyAnalysis(configFileName, "ForcePlotCLScan3");
  }
  
  //--------------------------------------//
  // Step 6.1: Calculate the test statistics:
  /*
  if (masterOption.Contains("TestStat") && 
      !masterOption.Contains("ResubmitTestStat")) {
    std::cout << "DHMaster: Step 6.1 - Calculating CL and p0." << std::endl;

    int jobCounterTS = 0;
    std::vector<TString> sigDHModes = m_config->getStrV("sigDHModes");
    for (int i_s = 0; i_s < (int)sigDHModes.size(); i_s++) {
      TString currSignal = sigDHModes[i_s];
      
      if (runInParallel) {
	submitTSViaBsub(fullConfigPath, m_config->getStr("TestStatOptions"), currSignal);
	jobCounterTS++;
	m_isFirstJob = false;
      }
      else {
	DHTestStat *ts = new DHTestStat(configFileName, currSignal,
					m_config->getStr("TestStatOptions"), NULL);
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
	submitTSViaBsub(fullConfigPath, m_config->getStr("TestStatOptions"), currSignal);
      	jobCounterTS++;
	m_isFirstJob = false;
      }
      else {
	DHTestStat *dhts = new DHTestStat(configFileName, currSignal,
					  m_config->getStr("TestStatOptions"), NULL);
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
  */
  
  //--------------------------------------//
  // Step 7.1: Calculate the limits on the dark matter signal strength.
  if (masterOption.Contains("MuLimit") &&
      !masterOption.Contains("ResubmitMuLimit")) {
    std::cout << "DHMaster: Step 7.1 - Calculate 95% CL mu value." << std::endl;

    //int jobCounterML = 0;
    
    //if (runInParallel) {
    //submitMLViaBsub(fullConfigPath, m_config->getStr("MuLimitOptions"), currSignal);
    //m_isFirstJob = false;
    //}
    //else {
    TString muCommand 
      = Form("./bin/%s %s %s", (m_config->getStr("exeMuLimit")).Data(),
	     fullConfigPath.Data(), m_config->getStr("MuLimitOptions").Data());
    std::cout << "Executing following system command: \n\t"
	      << muCommand << std::endl;
    system(muCommand);
    //}
    //jobCounterML++;
    //}
    //std::cout << "Submitted/completed " << jobCounterML << " jobs" << std::endl;
  }
  
  /*
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
	submitMLViaBsub(fullConfigPath, m_config->getStr("MuLimitOptions"), currSignal);
	m_isFirstJob = false;
      }
      else {
	system(Form(".%s/bin/%s %s %s %s", 
		    (m_config->getStr("packageLocation")).Data(), 
		    (m_config->getStr("exeMuLimit")).Data(),
		    fullConfigPath.Data(), currSignal.Data(),
		    m_config->getStr("MuLimitOptions").Data()));
      }
      jobCounterML++;
    }
    std::cout << "Resubmitted " << jobCounterML << " jobs" << std::endl;
  }
  */
  
  //--------------------------------------//
  // Step 8.1: Toss toys for a CL scan:
  if (masterOption.Contains("CLScanSubmitToys")) {
    std::cout << "DHMaster: Step 5.1 - Create pseudoexperiments." << std::endl;
    
    int toySeed = m_config->getInt("toySeed");
    int nToysTotal = m_config->getInt("nToysTotal");
    int nToyJobs = ((m_config->getNum("CLScanMax") - 
		     m_config->getNum("CLScanMin")) / 
		    m_config->getNum("CLScanStep"));
    for (int currIndex = 0; currIndex < nToyJobs; currIndex++) {
      TString currOptions = m_config->getStr("CLScanToyOptions");
      currOptions.Append(Form("%d",currIndex));
      submitPEViaBsub(fullConfigPath, currOptions, toySeed, nToysTotal);
      m_isFirstJob = false;
      toySeed++;
    }
    std::cout << "DHMaster: Submitted " << nToyJobs
	      << " total pseudo-experiments for CL scan." << std::endl;
  }
  
  //--------------------------------------//
  // Step 8.2: Create the CL scan plot:
  if (masterOption.Contains("CLScanAnalysis")) {
    std::cout << "DHMaster: Step 8 - Get ." << std::endl;
    system(Form("./bin/DHCLScan %s %s", fullConfigPath.Data(),
		(m_config->getStr("CLScanOptions")).Data()));
  }
  
  return 0;
}
