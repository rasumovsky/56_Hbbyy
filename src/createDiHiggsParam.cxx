////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  createDiHiggsParam.cxx                                                    //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 22/11/2015                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This main method provides a tool for performing individual and parameter- //
//  ized fits to the resonant di-Higgs Monte Carlo. Settings for the utility  //
//  are provided in data/signalParam/diHiggsParam.cfg                         //
//                                                                            //
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
   Checks whether the given mass value is contained in the vector, and adds it
   if it is not found.
   @param currList - The current vector of mass points.
   @param newMass - The mass point to check for membership in the list.
   @returns - A list of unique resonance mass points.
*/
std::vector<double> checkMResList(std::vector<double> currList, double newMass){
  for (int i_m = 0; i_m < (int)currList.size(); i_m++) {
    if (fabs(currList[i_m] - newMass) <= 0.1) return currList;
  }
  currList.push_back(newMass);
  return currList;
}

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
  std::cout << "createDiHiggsParam: ERROR! MxAOD doesn't have cutflow"
	    << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Get the mass of a resonance X->hh using the name of the input file.
   @param fileName - The name of the MxAOD file.
   @return - The resonance mass.
*/
double getResMassFromName(TString fileName) {
  for (int i_m = 100; i_m < 1000; i_m += 5) {
    if (fileName.Contains(Form("X%dtohh.root",i_m))) {
      return (double)i_m;
    }
  }
  std::cout << "createDiHiggsParam: ERROR getting X mass from file name."
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
  std::cout << "createDiHiggsParam: Making local copies of inputs."
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
  std::cout << "createDiHiggsParam: Removing local copies of inputs."
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
  std::cout << "createDiHiggsParam will run with parameters:"
	    << std::endl;
  settings->printDB();
  
  // Set the function type:
  TString function = settings->getStr("SignalFunctionalForm");
  
  // Check that output directory exists:
  TString outputDir = settings->getStr("OutputDir");
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Set the ATLAS Style for plots:
  SetAtlasStyle();
  
  // Instantiate SigParam class for individual & parameterized fits:
  SigParam *sps
    = new SigParam(settings->getStr("SampleName"), outputDir+"/Individual");
  SigParam *spp
    = new SigParam(settings->getStr("SampleName"), outputDir+"/Parameterized");
  // Set the print level of the signal parameterization tools:
  sps->verbosity(settings->getBool("Verbose"));// Sets the print level.
  spp->verbosity(settings->getBool("Verbose"));
  
  if (settings->getBool("UseDefinedParams")) {
    std::vector<TString> params = settings->getStrV("DefinedParams");
    for (int i_p = 0; i_p < (int)params.size(); i_p++) {
      sps->setParamState(params[i_p], settings->getStr("Param_"+params[i_p]));  
      spp->setParamState(params[i_p], settings->getStr("Param_"+params[i_p]));  
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
  double nTotEvt = 1000.0;
  TString currFileName = "";
  
  // Assign the MxAOD/TTree branches to variables:
  float v_mass;
  float v_weight;
  float v_xsbreff;
  double resMassToUse = 0.0;
  int v_cutFlow;
  int currCate;
  float v_myy;
  chain->SetBranchAddress(settings->getStr("MassBranchName"), &v_mass);
  chain->SetBranchAddress(settings->getStr("WeightBranchName"), &v_weight);
  chain->SetBranchAddress(settings->getStr("XSBREffBranchName"), &v_xsbreff);
  chain->SetBranchAddress(settings->getStr("CutFlowBranchName"), &v_cutFlow);
  chain->SetBranchAddress(settings->getStr("CategoryBranchName"), &currCate);
  chain->SetBranchAddress(settings->getStr("MyyBranchName"), &v_myy);
  int nEvents = chain->GetEntries();
  int nCategories = 0;
  std::vector<double> mResList; mResList.clear();
  std::map<TString,double> yieldWeighted; yieldWeighted.clear();
  std::map<TString,int> yieldUnweighted; yieldUnweighted.clear();
  
  //--------------------------------------//
  // Loop over events to build dataset for signal parameterization:
  std::cout << "There are " << nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {
    chain->GetEntry(index);
    printProgressBar(index, nEvents);

    // Change the nTotEvt normalization factor for each new file:
    if (!currFileName.EqualTo(chain->GetFile()->GetName())) {
      currFileName = chain->GetFile()->GetName();
      nTotEvt = getNTotEvtFromHist(chain->GetFile());
      resMassToUse = getResMassFromName(currFileName);
      yieldWeighted[Form("%2.2f",resMassToUse)] = 0.0;
      yieldUnweighted[Form("%2.2f",resMassToUse)] = 0.0;
    }
    
    // Only use events passing the selection:
    if (v_cutFlow < settings->getInt("CutFlowIndex")) continue;
    
    // The category index fed into the SigParam tool should start at 0.
    if (currCate < 0) continue;
    if (currCate >= nCategories) nCategories = currCate + 1;
    
    // The observed mass fed into the SigParam tool should always be in GeV.
    double massToUse = (double)v_mass / 1000.0;
    if (massToUse < 0.0) continue;
    
    // For the resonant analysis, also place a cut on diphoton mass:
    if (fabs((v_myy/1000.0) - settings->getNum("SMHiggsMass")) >
	settings->getNum("ResonantMyyWindowWidth")) {
      continue;
    }
    
    // Calculate the weight to use:
    //double weightToUse = luminosity * v_xsbreff * v_weight / nTotEvt;
    double weightToUse = luminosity * v_weight;
    
    // Add the mass and weight values to the datasets for fitting:
    sps->addMassPoint(resMassToUse, currCate, massToUse, weightToUse);
    spp->addMassPoint(resMassToUse, currCate, massToUse, weightToUse);

    // Add to the yield:
    yieldWeighted[Form("%2.2f",resMassToUse)] += weightToUse;
    yieldUnweighted[Form("%2.2f",resMassToUse)]++;
    
    // Add the mass point to the list of points:
    mResList = checkMResList(mResList, resMassToUse);
  }
  
  //--------------------------------------//
  // Now fit and plot the resonance shapes!
  std::cout << "createDiHiggsParam: Start fitting and plotting!" 
	    << std::endl;
  
  // Check which fits failed or succeeded:
  std::vector<TString> fitFailures; fitFailures.clear();
  std::vector<TString> fitSuccesses; fitSuccesses.clear();
  
  // Set the category names for output plots and tables:
  sps->nameTheCategories(settings->getStrV("CategoryNames"));
  spp->nameTheCategories(settings->getStrV("CategoryNames"));

  // Find the resonance mass range that was fitted:
  double minMass = minEntry(mResList);
  double maxMass = maxEntry(mResList);
  
  // List the parameters for the resonance function:
  std::vector<TString> parameters = spp->variablesForFunction(function);
  
  // Loop over the analysis categories:
  std::cout << "createDiHiggsParam: Loop over the categories."
	    << std::endl;
  for (int i_c = 0; i_c < nCategories; i_c++) {
    std::vector<TF1*> fFits; fFits.clear();
    std::vector<TGraphErrors*> gParams; gParams.clear();
    std::vector<TGraphErrors*> gRatio; gRatio.clear();
    
    for (int i_p = 0; i_p < (int)parameters.size(); i_p++) {
      TGraphErrors *currGParams	= new TGraphErrors();
      currGParams->SetNameTitle(Form("graph_%s",parameters[i_p].Data()),
				Form("graph_%s",parameters[i_p].Data()));
      TGraphErrors *currGRatio = new TGraphErrors();
      currGRatio->SetNameTitle(Form("ratio_%s",parameters[i_p].Data()),
			       Form("ratio_%s",parameters[i_p].Data()));
      gParams.push_back(currGParams);
      gRatio.push_back(currGRatio);
    }
    
    // Fit the individual mass points:
    std::cout << "createDiHiggsParam: # of masses: " 
	      << mResList.size() << std::endl;
    for (int i_m = 0; i_m < (int)mResList.size(); i_m++) {
      // The fit method returns true iff successful, else false:
      if (sps->makeSingleResonance(mResList[i_m], i_c, function)) {
	fitSuccesses
	  .push_back(Form("mass=%2.2f GeV in category %d",mResList[i_m],i_c));
	// Plot the fit result for the single resonance:
	sps->plotSingleResonance(mResList[i_m], i_c);
	// Save the parameter values in the graph:
	for (int i_p = 0; i_p < (int)parameters.size(); i_p++) {
	  double currVal
	    = sps->getParameterValue(parameters[i_p], mResList[i_m], i_c);
	  double currErr
	    = sps->getParameterError(parameters[i_p], mResList[i_m], i_c);
	  gParams[i_p]->SetPoint(i_m, mResList[i_m], currVal);
	  gParams[i_p]->SetPointError(i_m, 0.0, currErr);
	}
      }
      else {
	std::cout << "createDiHiggsParam: Fit at mRes="
		  << mResList[i_m] << ", cate=" << i_c << " did not converge :("
		  << std::endl;
	fitFailures
	  .push_back(Form("mass=%2.2f GeV in category %d",mResList[i_m],i_c));
      }
    }
    // Save the individual resonances to file:
    sps->saveAll();
    
    // Do a parameterized fit for the resonance:
    if (spp->makeCategoryParameterization(i_c, function)) {
      fitSuccesses.push_back(Form("parameterization in category %d", i_c));
      std::cout << "createDiHiggsParam: Simultaneous fit converged!" 
		<< std::endl;
      // Plot the parameterized resonances and yields, then save to file:
      spp->plotCategoryResonances(i_c);
      spp->plotYields(i_c);
      spp->saveAll();
      
      // Get the parameterization functions for each variable in TF1 format:
      for (int i_p = 0; i_p < (int)parameters.size(); i_p++) {
	TF1 *currTF1 = spp
	  ->createTF1FromParameterization(parameters[i_p],i_c,minMass,maxMass);
	fFits.push_back(currTF1);
      }
      
      // Start plotting the results:
      TCanvas *can = new TCanvas("can", "can", 800, 800);
      // Loop over the fit parameters:
      for (int i_p = 0; i_p < (int)parameters.size(); i_p++) {
	can->cd();
	can->Clear();
	// Format the pads:
	TPad *pad1 = new TPad("pad1", "pad1", 0.00, 0.33, 1.00, 1.00);
	TPad *pad2 = new TPad("pad2", "pad2", 0.00, 0.00, 1.00, 0.33);
	pad1->SetBottomMargin(0.00001);
	pad1->SetBorderMode(0);
	pad2->SetTopMargin(0.00001);
	pad2->SetBottomMargin(0.4);
	pad2->SetBorderMode(0);
	can->cd();
	pad1->Draw();
	pad2->Draw();
	pad1->cd();
	// Format the plot axes and draw options:
	fFits[i_p]->GetXaxis()->SetTitle("M_{Resonance} [GeV]");
	gParams[i_p]->GetYaxis()
	  ->SetTitle(Form("parameter %s", parameters[i_p].Data()));
	gParams[i_p]->GetYaxis()->SetTitleOffset(1.2);
	fFits[i_p]->SetLineWidth(2);
	fFits[i_p]->SetLineColor(kRed);
	gParams[i_p]->Draw("AEP");
	fFits[i_p]->Draw("LSAME");
	// Write some text:
	TLatex text; text.SetNDC(); text.SetTextColor(1);
	text.DrawLatex(0.2, 0.88, Form("category %d", i_c));
	double chi2Prob = TMath::Prob(gParams[i_p]->Chisquare(fFits[i_p]),
				      mResList.size());
	text.DrawLatex(0.2, 0.82, Form("#chi^{2} prob. = %2.2f",chi2Prob));
	
	// Create a ratio plot for individual fits and parameterization:
	pad2->cd();
	for (int i_n = 0; i_n < gParams[i_p]->GetN(); i_n++) {
	  double xVal; double yVal;
	  gParams[i_p]->GetPoint(i_n, xVal, yVal);
	  double fVal = fFits[i_p]->Eval(xVal);
	  gRatio[i_p]->SetPoint(i_n, xVal, (yVal/fVal));
	  double yErr = gParams[i_p]->GetErrorY(i_n);
	  double ratioErr = (yErr / yVal);
	  gRatio[i_p]->SetPointError(i_n, 0.0, ratioErr);
	}
	// Format the ratio plot:
	gRatio[i_p]->GetXaxis()->SetTitle("M_{Resonance} [GeV]");
	gRatio[i_p]->GetYaxis()->SetTitle("Ratio");
	gRatio[i_p]->GetXaxis()->SetTitleOffset(1.5);
	gRatio[i_p]->GetYaxis()->SetTitleOffset(0.6);
	gRatio[i_p]->GetXaxis()->SetTitleSize(0.1);
	gRatio[i_p]->GetYaxis()->SetTitleSize(0.1);
	gRatio[i_p]->GetXaxis()->SetLabelSize(0.1);
	gRatio[i_p]->GetYaxis()->SetLabelSize(0.1);
	gRatio[i_p]->GetXaxis()
	  ->SetRangeUser(gParams[i_p]->GetXaxis()->GetXmin(), 
			 gParams[i_p]->GetXaxis()->GetXmax());
	gRatio[i_p]->GetYaxis()->SetNdivisions(4);
  	gRatio[i_p]->Draw("AEP");
	// Draw a line at a ratio of 1.0:
	TLine *line = new TLine();
	line->SetLineStyle(2);
	line->SetLineWidth(1);
	line->SetLineColor(1);
	line->DrawLine(gRatio[i_p]->GetXaxis()->GetXmin(), 1,
		       gRatio[i_p]->GetXaxis()->GetXmax(), 1);
	
	can->Print(Form("%s/plotParam_%s_cate%d.eps", outputDir.Data(),
			  parameters[i_p].Data(), i_c));
	can->Clear();
      }
      delete can;
    }
    else {
      std::cout << "createDiHiggsParam: Parameterization failed cate="
		<< i_c << std::endl;
      fitFailures.push_back(Form("parameterization in category %d", i_c));
    }
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
  std::cout << "createDiHiggsParam: Printing inclusive yields:" << std::endl;
  for (int i_m = 0; i_m < (int)mResList.size(); i_m++) {
    std::cout << "\tmass = " << mResList[i_m] << ", yield = " 
	      << yieldWeighted[Form("%2.2f",mResList[i_m])] << "("
	      << yieldUnweighted[Form("%2.2f",mResList[i_m])] << ")" 
	      << std::endl;
  }
  
  return 0;
}

