////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHMaster.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 08/02/2016                                                          //
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
//    - CLScanAnalysis                                                        //
//    - CLScanSubmitToys                                                      //
//    - PlotCLvsMX                                                            //
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
  
  TString exeAna = m_config->getStr("AnalysisType");
  
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
  
  TString exeAna = m_config->getStr("AnalysisType");
  
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
   @param exeSeed - the seed for the randomized dataset generation.
   @param exeToysPerJob - the number of toy datasets to create per job.
   @param resonanceMass - The mass of the resonance (resonant analysis only)
*/
void submitPEViaBsub(TString exeConfigFile, TString exeOption, int exeSeed, 
		     int exeToysPerJob, int resonanceMass) {
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
  TString nameJScript = Form("%s/jobFilePseudoExp.sh %s %s %s %s %s %d %d %d", 
			     exe.Data(), (m_config->getStr("JobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exePseudoExp")).Data(),
			     exeOption.Data(), exeSeed, exeToysPerJob,
			     resonanceMass);
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   This is the main DHMaster method. Runs every analysis component.
   @param masterOption - The task to run (see document header or README).
   @param configFileName - The name of the config file. 
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
    sys->addCategoryNames(m_config->getStrV("SysCateNames"));
    std::vector<TString> sysComponents = m_config->getStrV("SysComponents");
    for (int i_s = 0; i_s < (int)sysComponents.size(); i_s++) {
      TString fileName = "";
      if (m_config->getStr("AnalysisType").EqualTo("Resonant")) {
	fileName = Form("%s/%s.txt", (m_config->getStr("SysDirectory")).Data(),
			(sysComponents[i_s]).Data());
      }
      // Non-resonant files:
      else {
	fileName = Form("%s/%s_allSys.txt",
			(m_config->getStr("SysDirectory")).Data(),
			(sysComponents[i_s]).Data());
      }
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
    
    // Also create a list of systematic uncertainties to ignore:
    std::vector<TString> sysToIgnore; sysToIgnore.clear();
    sysToIgnore.push_back("FT_EFF_Eigen_B_0");
    sysToIgnore.push_back("FT_EFF_Eigen_B_1");
    sysToIgnore.push_back("FT_EFF_Eigen_B_2");
    sysToIgnore.push_back("FT_EFF_Eigen_B_3");
    sysToIgnore.push_back("FT_EFF_Eigen_B_4");
    sysToIgnore.push_back("FT_EFF_Eigen_B_5");
    sysToIgnore.push_back("FT_EFF_Eigen_C_0");
    sysToIgnore.push_back("FT_EFF_Eigen_C_1");
    sysToIgnore.push_back("FT_EFF_Eigen_C_2");
    sysToIgnore.push_back("FT_EFF_Eigen_C_3");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_0");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_1");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_2");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_3");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_4");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_5");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_6");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_7");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_8");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_9");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_10");
    sysToIgnore.push_back("FT_EFF_Eigen_Light_11");
    sysToIgnore.push_back("JET_GroupedNP_1");
    sysToIgnore.push_back("JET_GroupedNP_2");
    sysToIgnore.push_back("JET_GroupedNP_3");
    sysToIgnore.push_back("MUON");
    sys->ignoreSystematics(sysToIgnore);

    // Then parameterize the resonant signal uncertainties:
    if (m_config->getStr("AnalysisType").EqualTo("Resonant")) {
      std::map<TString,double> sampleMap; sampleMap.clear();
      sampleMap["X275-Window"] = 275.0;
      sampleMap["X300-Window"] = 300.0;
      sampleMap["X325-Window"] = 325.0;
      sampleMap["X350-Window"] = 350.0;
      sampleMap["X400-Window"] = 400.0;
      sys->parameterizeSyst("SigBSM2H", sampleMap);
      sys->parameterizeSyst("SigSM", sampleMap);
    }
    
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
		      i_s, nToysPerJob, 0);
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
    std::vector<double> scanMXValues = m_config->getNumV("MXScanValues");
    DHToyAnalysis *dhta
      = new DHToyAnalysis(configFileName,"NONE", (int)(scanMXValues[0]));
    //DHToyAnalysis *dhta = new DHToyAnalysis(configFileName, "ForcePlotCLScan3", (int)(scanMXValues[0]));
  }
  
  //--------------------------------------//
  // Step 6.0: Calculate CL and p0 statistics with asymptotic formulae:
  if (masterOption.Contains("TestStat")) {
    std::cout << "DHMaster: Step 6.0 - Asymptotic CL and p0 calculation."
	      << std::endl;
    // Instantiate the statistics class:
    DHTestStat *dhts = new DHTestStat(configFileName, 
				      m_config->getStr("TestStatOptions"),
				      NULL);
    dhts->setParam(m_config->getStr("CrossSectionVar"),
		   m_config->getNum("CrossSectionValue"), true);
    
    // Non-resonant analysis: specify cross-section:
    
    // Resonant analysis: specify cross-section and resonance mass:
    if ((m_config->getStr("AnalysisType")).EqualTo("Resonant")) {
      dhts->setParam(m_config->getStr("ResonanceMassVar"),
		     m_config->getNum("ResonanceMassValue"), true);
    }
    
    // Calculate CL and then p0:
    dhts->calculateNewCL();
    dhts->calculateNewP0();
    
    // Check that fits succeeded:
    if (!dhts->fitsAllConverged()) {
      std::cout << "\nDMMaster: ERROR! Statistics fits failed.\n" << std::endl;
    }
  }
  
  //--------------------------------------//
  // Step 8.1: Toss toys for a CL scan:
  if (masterOption.Contains("CLScanSubmitToys")) {
    std::cout << "DHMaster: Step 8.1 - Create toys for CL Scan." << std::endl;
    
    int toySeed = m_config->getInt("toySeed");
    int nToysTotal = m_config->getInt("nToysTotal");
    int nToysPerJob = m_config->getInt("nToysPerJob");
    int highestSeed = toySeed + nToysTotal;
    
    // Loop over scan points:
    std::vector<double> scanCLValues = m_config->getNumV("CLScanValues");
    for (int currIndex = 0; currIndex < (int)scanCLValues.size(); currIndex++) {
      TString currOptions = m_config->getStr("CLScanToyOptions");
      currOptions.Append(Form("%d",currIndex));
      
      // Then split the toy jobs for each scan point:
      for (int i_s = toySeed; i_s < highestSeed; i_s += nToysPerJob) {
	// Loop over resonance masses for resonant analysis:
	if ((m_config->getStr("AnalysisType")).EqualTo("Resonant")) {
	  std::vector<double> scanMXValues = m_config->getNumV("MXScanValues");
	  for (int i_m = 0; i_m < (int)scanMXValues.size(); i_m++) {
	    submitPEViaBsub(fullConfigPath, currOptions, i_s, nToysPerJob, 
			    (int)(scanMXValues[i_m]));
	  }
	}
	else {
	  submitPEViaBsub(fullConfigPath, currOptions, i_s, nToysPerJob, 0);
	}
	m_isFirstJob = false;
      }
    }
    int nJobSubmits = ((int)scanCLValues.size()*(int)(nToysTotal/nToysPerJob));
    if ((m_config->getStr("AnalysisType")).EqualTo("Resonant")) {
      std::vector<double> scanMXValues = m_config->getNumV("MXScanValues");
      nJobSubmits = (((int)scanCLValues.size()*(int)(nToysTotal/nToysPerJob)) * 
		     (int)(scanMXValues.size()));
    }
    std::cout << "DHMaster: Submitted " << nJobSubmits << " total CL scan toys."
	      << std::endl;
  }
  
  //--------------------------------------//
  // Step 8.2: Create the CL scan plot:
  if (masterOption.Contains("CLScanAnalysis")) {
    std::cout << "DHMaster: Step 8.2 - Get CL Scan Plot" << std::endl;
    
    if ((m_config->getStr("AnalysisType")).EqualTo("Resonant")) {
      
      std::vector<double> scanMXValues = m_config->getNumV("MXScanValues");
      for (int i_m = 0; i_m < (int)scanMXValues.size(); i_m++) {
	system(Form("./bin/DHCLScan %s %s %d", fullConfigPath.Data(),
		    (m_config->getStr("CLScanOptions")).Data(),
		    (int)(scanMXValues[i_m])));
      }
    }
    else {
      system(Form("./bin/DHCLScan %s %s 0", fullConfigPath.Data(),
		  (m_config->getStr("CLScanOptions")).Data()));
    }
  }
  
  //--------------------------------------//
  // Step 8.2: Create the CL scan plot:
  if (masterOption.Contains("PlotCLvsMX")) {
    std::cout << "DHMaster: Step 8.2 - Get Plot of CL vs. MX" << std::endl;
    system(Form("./bin/DHPlotCLvsMX %s %s", fullConfigPath.Data(),
		(m_config->getStr("CLMXPlotOptions")).Data()));
  }
  
  return 0;
}
