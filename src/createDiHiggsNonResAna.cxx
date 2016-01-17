////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  createDiHiggsNonResAna.cxx                                                //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 03/12/2015                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This main method provides a tool for performing individual fits to the    //
//  single Standard Model Higgs background in the non-resonant di-Higgs       //
//  search.                                                                   //
//                                                                            //
//  Associated config file: data/signalParam/diHiggsNonResAna.cfg             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamTools/HggTwoSidedCBPdf.h"
#include "HGamTools/SigParam.h"
#include "HGamTools/AtlasStyle.h"

#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"

/**
   -----------------------------------------------------------------------------
   Retrieve the total number of events in the file for event normalization:
*/
double getNTotEvtFromHist(TFile *f) {
  // Find the cutflow histograms from the file based on limited name info:
  TIter next(f->GetListOfKeys());
  TObject *currObj;
  while ((currObj = (TObject*)next())) {
    TString currName = currObj->GetName();
    if (currName.Contains("CutFlow") && currName.Contains("weighted")
    	&& currName.Contains("noDalitz")) {
      return (((TH1F*)f->Get(currName))->GetBinContent(3));
    }
  }
  std::cout << "createDiHiggsNonResAna: ERROR! MxAOD doesn't have cutflow"
	    << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Copy files from a slow resource (e.g. EOS) to the local disk for faster
   processing.
   @param fileNames - The original file names.
   @returns - An updated list of file names.
*/
std::vector<TString> makeLocalFileCopies(std::vector<TString> fileNames) {
  std::cout << "createDiHiggsNonResAna: Making local copies of inputs."
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
  std::cout << "createDiHiggsNonResAna: Removing local copies of inputs."
	    << std::endl;
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    system(Form("rm %s", fileNames[i_f].Data()));
  }
}

/**
   -----------------------------------------------------------------------------
   Returns the maximum entry in a vector of doubles.
   @param currList - The vector of doubles.
   @returns - The maximum entry in the vector.
*/
double maxEntry(std::vector<double> currList) {
  double maximum = 0.0;
  for (std::vector<double>::iterator mIter = currList.begin(); 
       mIter != currList.end(); mIter++) {
    if (*mIter > maximum) maximum = *mIter;
  }
  return maximum;
}

/**
   -----------------------------------------------------------------------------
   Returns the minimum entry in a vector of doubles.
   @param currList - The vector of doubles.
   @returns - The minimum entry in the vector.
*/
double minEntry(std::vector<double> currList) {
  double minimum = 100000.0;
  for (std::vector<double>::iterator mIter = currList.begin(); 
       mIter != currList.end(); mIter++) {
    if (*mIter < minimum) minimum = *mIter;
  }
  return minimum;
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
   The main method for this utility. Provide 1 argument - the location of the 
   config (.cfg) file, which should be stored in the data/ directory. The main()
   method runs over the samples provided, performs the fits requests, and gives
   comparisons of parameterized and non-parameterized fits. 
*/
int main(int argc, char *argv[])
{
  // Check that the config file location is provided.
  if (argc < 2) HG::fatal("No arguemnts provided");
  HG::Config *settings = new HG::Config(TString(argv[1]));

  // Print configuration for benefit of user:
  std::cout << "createDiHiggsNonResAna will run with parameters:"
	    << std::endl;
  settings->printDB();
  
  // Set the function type:
  TString function = settings->getStr("SignalFunctionalForm");
  
  // Check that output directory exists:
  TString outputDir = settings->getStr("OutputDir");
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Set the ATLAS Style for plots:
  SetAtlasStyle();
  
  // Get the type of analysis ("Resonant" or "NonResonant"):
  TString analysisType = settings->getStr("AnalysisType");
  
  // Instantiate SigParam class for individual & parameterized fits:
  SigParam *sps
    = new SigParam(settings->getStr("SampleName"), outputDir+"/Individual");
  sps->verbosity(settings->getBool("Verbose"));// Sets the print level.
  sps->doBinnedFit(false, settings->getNum("GeVPerPlotBin"));
  
  // Set initial values and ranges if desired:
  if (settings->getBool("UseDefinedParams")) {
    std::vector<TString> params = settings->getStrV("DefinedParams");
    for (int i_p = 0; i_p < (int)params.size(); i_p++) {
      sps->setParamState(params[i_p], settings->getStr("Param_"+params[i_p]));  
    }
  }
  
  // Prepare for loop over input MxAOD/TTree:
  std::vector<TString> fileNames = settings->getStrV("InputFile");
  // Make local copies of files if requested, to improve speed:
  if (settings->getBool("MakeLocalCopies")) {
    fileNames = makeLocalFileCopies(fileNames);
  }
  // Create TChain of input files:
  TChain *chain = new TChain(settings->getStr("TreeName"));
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    chain->AddFile(fileNames[i_f]);
  }
  
  // Get the luminosity and reserve a variable for total event norm:
  double luminosity = settings->getNum("Luminosity") / 1000.0;
  //double nTotEvt = 1000.0;
  TString currFileName = "";
  
  // Assign a dummy resonance mass:
  double resonanceMass = settings->getNum("SMHiggsMass");
  
  // Assign the MxAOD/TTree branches to variables:
  float v_myyjj;
  float v_weight;
  float v_xsbreff;
  int v_cutFlow;
  int v_currCate;
  float v_myy;
  chain->SetBranchAddress(settings->getStr("MassBranchName"), &v_myyjj);
  chain->SetBranchAddress(settings->getStr("WeightBranchName"), &v_weight);
  chain->SetBranchAddress(settings->getStr("XSBREffBranchName"), &v_xsbreff);
  chain->SetBranchAddress(settings->getStr("CutFlowBranchName"), &v_cutFlow);
  chain->SetBranchAddress(settings->getStr("CategoryBranchName"), &v_currCate);
  chain->SetBranchAddress(settings->getStr("MyyBranchName"), &v_myy);
  int nEvents = chain->GetEntries();
  int nCategories = 0;
  std::map<TString,double> yieldWeighted; yieldWeighted.clear();
  std::map<TString,int> yieldUnweighted; yieldUnweighted.clear();

  std::map<int,double> yieldInWindow; yieldInWindow.clear();
  std::map<int,double> yieldTotal; yieldTotal.clear();
  
  //--------------------------------------//
  // Loop over events to build dataset for signal parameterization:
  std::cout << "There are " << nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {
    chain->GetEntry(index);
    printProgressBar(index, nEvents);

    // Change the nTotEvt normalization factor for each new file:
    if (!currFileName.EqualTo(chain->GetFile()->GetName())) {
      currFileName = chain->GetFile()->GetName();
      //nTotEvt = getNTotEvtFromHist(chain->GetFile());
      yieldWeighted[Form("%2.2f",resonanceMass)] = 0.0;
      yieldUnweighted[Form("%2.2f",resonanceMass)] = 0.0;
    }
    
    // Only use events passing the selection:
    if (v_cutFlow < settings->getInt("CutFlowIndex")) continue;
    
    // The category index fed into the SigParam tool should start at 0.
    int currCate = v_currCate;
    if (currCate < 0) continue;
    if (currCate >= nCategories) nCategories = currCate + 1;
    
    // The observed mass fed into the SigParam tool (should always be in GeV):
    double massToUse = -100.0;// update depending on analysis type below.
        
    // For the resonant analysis, also place a cut on diphoton mass:
    if (analysisType.EqualTo("Resonant")) {
      massToUse = (double)v_myyjj / 1000.0;
      if (fabs((v_myy/1000.0) - settings->getNum("SMHiggsMass")) >
	  settings->getNum("ResonantMyyWindowWidth")) {
	continue;
      }
    }
    else if (analysisType.EqualTo("NonResonant")) {
      massToUse = (double)v_myy/1000.0;
    }
    else {
      std::cout << "createDiHiggsParam: ERROR! " << analysisType
		<< " not recognized as analysis type..." << std::endl;
      exit(0);
    }
    
    // Don't use events with negative mass:
    if (massToUse < 0.0) continue;
    
    // Calculate the weight to use:
    //double weightToUse = luminosity * v_xsbreff * v_weight / nTotEvt;
    double weightToUse = luminosity * v_weight;
    
    // Add the mass and weight values to the datasets for fitting:
    sps->addMassPoint(resonanceMass, currCate, massToUse, weightToUse);
    
    // Add to the yield:
    yieldWeighted[Form("%2.2f",resonanceMass)] += weightToUse;
    yieldUnweighted[Form("%2.2f",resonanceMass)]++;
    
    // Also count events inside myy window:
    if (yieldTotal.count(0) == 0) {
      yieldTotal[currCate] = 0.0;
      yieldInWindow[currCate] = 0.0;
    }
    
    yieldTotal[currCate]++;
    if (fabs((v_myy/1000.0) - settings->getNum("SMHiggsMass")) <
	settings->getNum("ResonantMyyWindowWidth")) {
      yieldInWindow[currCate]++;
    }
  }
  
  //--------------------------------------//
  // Now fit and plot the resonance shapes!
  std::cout << "createDiHiggsNonResAna: Start fitting and plotting!" 
	    << std::endl;
  
  // Check which fits failed or succeeded:
  std::vector<TString> fitFailures; fitFailures.clear();
  std::vector<TString> fitSuccesses; fitSuccesses.clear();
  
  // Set the category names for output plots and tables:
  std::vector<TString> categoryNames = settings->getStrV("CategoryNames");
  sps->nameTheCategories(categoryNames);
  
  // Loop over the analysis categories:
  std::cout << "createDiHiggsNonResAna: Loop over the categories."
	    << std::endl;
  for (int i_c = 0; i_c < nCategories; i_c++) {
    
    // The fit method returns true iff successful, else false:
    if (sps->makeSingleResonance(resonanceMass, i_c, function)) {
      fitSuccesses.push_back(Form("mass=%2.2f GeV in category %d", 
				  resonanceMass, i_c));
      // Plot the fit result for the single resonance:
      sps->plotSingleResonance(resonanceMass, i_c);

    }
    else {
      std::cout << "createDiHiggsNonResAna: Fit at mRes=" << resonanceMass 
		<< ", cate=" << i_c << " did not converge :(" << std::endl;
      fitFailures.push_back(Form("mass=%2.2f GeV in category %d",
				 resonanceMass, i_c));
    }
    
    // Save the individual resonances to file:
    sps->saveAll();
  }
  
  // Remove any files copied locally for speed, then exit:
  if (settings->getBool("MakeLocalCopies")) removeLocalFileCopies(fileNames);
  
  // Print out a list of good and bad fits:
  std::cout << "createSingleSignal: Printing Summary" << std::endl;
  std::cout << "\tFits that succeeded (" << (int)fitSuccesses.size() 
	    << " total):" << std::endl;
  for (int i_s = 0; i_s < (int)fitSuccesses.size(); i_s++) {
    std::cout << "\t\t" << fitSuccesses[i_s] << std::endl;
  }
  std::cout << "\tFits that failed (" << (int)fitFailures.size() << " total):" 
	    << std::endl;
  for (int i_f = 0; i_f < (int)fitFailures.size(); i_f++) {
    std::cout << "\t\t" << fitFailures[i_f] << std::endl;
  }
  
  // Print yield information:
  std::cout << "createDiHiggsNonResAna: Printing inclusive yields from counting"
	    << std::endl;
  std::cout << "\tmass = " << resonanceMass << ", yield = " 
	    << yieldWeighted[Form("%2.2f",resonanceMass)] << " ("
	    << yieldUnweighted[Form("%2.2f",resonanceMass)] << ")" << std::endl;
  
  // Get yields in each category from the signal parameterization:
  std::cout << "createDiHiggsNonResAna: Print yields in each category from tool"
	    << std::endl;
  for (int i_c = 0; i_c < (int)categoryNames.size(); i_c++) {
    std::cout << "\t" << categoryNames[i_c] << " \t"
	      << sps->getYieldInCategory(resonanceMass,i_c) << "\n" 
	      << std::endl;
  }
  
  // Finally, print the resonance variable values:
  std::cout << "createDiHiggsParam: Printing parameterized resonance at mH="
	    << resonanceMass << " GeV" << std::endl;
  // List the parameters for the resonance function:
  std::vector<TString> parameters = sps->variablesForFunction(function);
  for (int i_c = 0; i_c < nCategories; i_c++) {
    std::cout << "\tCategory = " << i_c << std::endl;
    for (int i_p = 0; i_p < (int)parameters.size(); i_p++) {
      std::cout << "\t\t" << parameters[i_p] << "\t" 
		<< sps->getParameterValue(parameters[i_p], resonanceMass, i_c)
		<< std::endl;
    }
    // Then calculate the standard deviation:  
    std::cout << "\t\t68% interval is " 
	      << sps->getMeanOrStdDev("StdDev", resonanceMass, i_c)
	      << std::endl;
  }
  
  // Then calculate the fraction of events in the myy window:
  std::cout << "createDiHiggsNonResAna: Calculating myy efficiency."
	    << std::endl;
  for (int i_c = 0; i_c < nCategories; i_c++) {
    std::cout << "\tcategory " << i_c << std::endl;
    std::cout << "\t  MC = " << (yieldInWindow[i_c] / yieldTotal[i_c])
	      << std::endl;
    
    double width = settings->getNum("ResonantMyyWindowWidth");
    double numFit = sps->getYieldInWindow(resonanceMass, i_c,
					  resonanceMass-width, 
					  resonanceMass+width);
    double denomFit = sps->getYieldInCategory(resonanceMass, i_c);
    std::cout << "\t  fit = " << (numFit / denomFit) << std::endl;
  }
  return 0;
}
