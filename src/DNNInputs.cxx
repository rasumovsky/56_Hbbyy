////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DNNInputs.cxx                                                             //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 30/08/2016                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This program converts TTree data to text file data that can more easily   //
//  be used for the neural network analysis.                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "MiniNtup.h"

// Globally scoped variables:
TString m_outputDir;
TString m_options;
MiniNtup *m_tree;

/**
   -----------------------------------------------------------------------------
   Copy files from a slow resource (e.g. EOS) to the local disk for faster
   processing.
   @param fileNames - The original file names.
   @return - An updated list of file names.
*/
std::vector<TString> makeLocalFileCopies(std::vector<TString> fileNames) {
  std::cout << "DNNInputs: Making local copies of inputs."
	    << std::endl;
  std::vector<TString> result; result.clear();
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    TString newName = Form("tempFile%d.root", i_f);
    if (fileNames[i_f].Contains("root://eosatlas/")) {
      system(Form("xrdcp %s %s", fileNames[i_f].Data(), newName.Data()));
    }
    else if (fileNames[i_f].Contains("/eos/atlas/")) {
      system(Form("eos cp %s %s", fileNames[i_f].Data(), newName.Data()));
    }
    else {
      system(Form("cp %s %s", fileNames[i_f].Data(), newName.Data()));
    }
    result.push_back(newName);
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Remove any files that were copied over for speed.
   @param fileNames - The original file names.
*/
void removeLocalFileCopies(std::vector<TString> fileNames) {
  std::cout << "DNNInputs: Removing local copies of inputs."
	    << std::endl;
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    system(Form("rm %s", fileNames[i_f].Data()));
  }
}

/**
   -----------------------------------------------------------------------------
   Prints a progress bar to screen to provide elapsed time and remaining time
   information to the user. This is useful when processing large datasets. 
   @param index - The current event index.
   @param total - The total number of events.
*/
void printProgressBar(int index, int total) {
  if (index%100000 == 0) {
    TString print_bar = " [";
    for (int bar = 0; bar < 20; bar++) {
      double current_fraction = double(bar) / 20.0;
      if (double(index)/double(total) > current_fraction) print_bar.Append("/");
      else print_bar.Append(".");
    }
    print_bar.Append("] ");
    double percent = 100.0 * (double(index) / double(total));
    TString text = Form("%s %2.2f ", print_bar.Data(), percent);
    std::cout << text << "%\r" << std::flush; 
  }
}

/**
   -----------------------------------------------------------------------------
   The main method for this utility.
   @param configFile - The configuration file for this code (usually in data/)
   @param options - The job options for execution.
*/
int main(int argc, char *argv[])
{
  // Check that arguments are provided.
  if (argc < 2) {
    std::cout << "\nUsage: " << argv[0] << " <configFile> <options>"
	      << std::endl;
    exit(0);
  }
  
  Config *config = new Config(TString(argv[1]));
  m_options = argv[2];

  // Check that output directory exists:
  m_outputDir = Form("%s/%s/DNNInputs", (config->getStr("MasterOutput")).Data(),
		     (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Define output file with variables:
  std::ofstream outputDNN(Form("%s/variablesForDNN_target%d.txt",
			       m_outputDir.Data(),
			       config->getInt("TrainingLabel")));
  
  // Prepare for loop over input MxAOD/TTree:
  std::vector<TString> fileNames = config->getStrV("TTreesForData");
  // Make local copies of files if requested, to improve speed:
  if (config->getBool("MakeLocalMxAODCopies")) {
    fileNames = makeLocalFileCopies(fileNames);
  }
  // Create TChain of input files:
  TChain *chain = new TChain(config->getStr("TreeName"));
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    chain->AddFile(fileNames[i_f]);
  }
  m_tree = new MiniNtup(chain);
  
  //--------------------------------------//
  // Loop over events:
  int nEvents = m_tree->fChain->GetEntries();
  std::cout << "There are " << nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {

    // Load event from MxAOD:
    m_tree->fChain->GetEntry(index);
    printProgressBar(index, nEvents);
    
    //---------- Select the event ----------//
    
    // The section below takes low-level physics variables and converts them
    // into a number between -0.5 and 0.5 with a mean ~0. 
    
    // Loop over photons:
    for (int i_p = 0; i_p < 2; i_p++) {
      outputDNN << (m_tree->photon_eta[i_p] / 4.0) << " ";
      outputDNN << (m_tree->photon_phi[i_p] / (2.0*TMath::Pi())) << " ";
      outputDNN << (TMath::Log10(m_tree->photon_pt[i_p]+1.0)-1.9)/2.0 << " ";
      outputDNN << ((double)(m_tree->photon_isTight[i_p]) - 0.5) << " ";
      outputDNN << ((TMath::Log10(0.001*m_tree->photon_ptcone20[i_p]+1.0) - 3.1)
		    / 2.0) << " ";
      outputDNN << (TMath::Log10(m_tree->photon_topoEtcone40[i_p] + 2000) - 3.3)
		<< " ";
    }
    
    // Loop over jets:
    outputDNN << ((double)(m_tree->jet_n) - 4.0) / 5.0 << " ";
    for (int i_j = 0; i_j < 5; i_j++) {
      if (i_j <= m_tree->jet_n) {
	outputDNN << TMath::Log10(m_tree->jet_pt[i_j] + 1.0) / 4.0 << " ";
	outputDNN << (m_tree->jet_eta[i_j] / 8.0) << " ";
	outputDNN << (m_tree->jet_phi[i_j] / (2.0*TMath::Pi())) << " ";
	outputDNN << TMath::Log10(m_tree->jet_m[i_j] + 1.0) / 2.0 << " ";
	outputDNN << TMath::Log10(1.5 - m_tree->jet_Jvt[i_j]) << " ";
	outputDNN << ((double)(m_tree->jet_MV2c10_FixedCutBEff_70[i_j]) - 0.5)
		  << " ";
      }
      else {
	outputDNN << -1.0 << " ";
	outputDNN << -1.0 << " ";
	outputDNN << -1.0 << " ";
	outputDNN << -1.0 << " ";
	outputDNN << -1.0 << " ";
	outputDNN << -1.0 << " ";
      }
    }
    
    outputDNN << (TMath::Log10(m_tree->met+1.0)-1.5) / 2.0 << " ";
    outputDNN << (m_tree->met_phi / (2.0*TMath::Pi())) << " ";
    outputDNN << (TMath::Log10(m_tree->met_sumet+1.0) - 2.7) / 3.0 << " ";
    outputDNN << config->getInt("TrainingLabel") << std::endl;
    
  }// End of loop over events
  outputDNN.close();
  std::cout << "End of event loop" << std::endl;

  //----------------------------------------//
  // Remove files that were copied:
  if (config->getBool("MakeLocalMxAODCopies")) removeLocalFileCopies(fileNames);
  
  delete m_tree;
  delete chain;
  delete config;
  return 0;
}
