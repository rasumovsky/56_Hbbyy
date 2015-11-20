////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  resDiHiggsParameterization.cxx                                          //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 08/11/2015                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This main method provides a tool for performing individual and parameter- //
//  ized fits to the resonance Monte Carlo. Settings for the utility are      //
//  provided in signalParamExample.cfg.                                       //
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
  std::cout << "resDiHiggsParameterization will run with parameters:"
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
  
  /*
  std::vector<TString> params = settings->getStrV("DefinedParams");
  for (int i_p = 0; i_p < (int)params.size(); i_p++) {
    spp->setParamState(params[i_p], settings->getStr("Param_"+params[i_p]));  
  }
  */
  
  std::vector<TString> categoryNames = settings->getStrV("CategoryNames");
  int nCategories = categoryNames.size();
  std::vector<double> mResList = settings->getNumV("MassPoints");
  
  //--------------------------------------//
  // Loop over mass points to get files, and then categories.
  for (std::vector<double>::iterator iterMass = mResList.begin(); 
       iterMass != mResList.end(); iterMass++) {
    TString fileNameForm = settings->getStr("InputFileForm");
    fileNameForm.ReplaceAll("MASS",Form("%d",((int)*iterMass)));
    TFile *tempFile = new TFile(Form(fileNameForm));
    
    // Loop over categories:
    for (int i_c = 0; i_c < (int)categoryNames.size(); i_c++) {
      TString histName = Form("Hist%s",(categoryNames[i_c]).Data());
      TH1F *currHist = (TH1F*)tempFile->Get(settings->getStr(histName));
      // Loop over bins, treating each as a mass point:
      for (int i_b = 1; i_b <= currHist->GetXaxis()->GetNbins(); i_b++) {
	sps->addMassPoint(*iterMass, i_c, currHist->GetBinCenter(i_b), currHist->GetBinContent(i_b));
	spp->addMassPoint(*iterMass, i_c, currHist->GetBinCenter(i_b), currHist->GetBinContent(i_b));
      }
      delete currHist;
    }
    tempFile->Close();
    delete tempFile;
  }
    
  //--------------------------------------//
  // Now fit and plot the resonance shapes!
  std::cout << "resDiHiggsParameterization: Start fitting and plotting!" 
	    << std::endl;

  // Check which fits failed or succeeded:
  std::vector<TString> fitFailures; fitFailures.clear();
  std::vector<TString> fitSuccesses; fitSuccesses.clear();
  
  // Set the category names for output plots and tables:
  sps->nameTheCategories(settings->getStrV("CategoryNames"));
  spp->nameTheCategories(settings->getStrV("CategoryNames"));

  // Find the resonance mass range that was fitted:
  double minMass = mResList[0];
  double maxMass = mResList[(int)mResList.size()-1];
  
  // List the parameters for the resonance function:
  std::vector<TString> parameters = spp->variablesForFunction(function);
  
  // Loop over the analysis categories:
  std::cout << "resDiHiggsParameterization: Loop over the categories."
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
    std::cout << "resDiHiggsParameterization: # of masses: " 
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
	std::cout << "resDiHiggsParameterization: Fit at mRes="
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
      std::cout << "resDiHiggsParameterization: Simultaneous fit converged!" 
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
      std::cout << "resDiHiggsParameterization: Parameterization failed cate="
		<< i_c << std::endl;
      fitFailures.push_back(Form("parameterization in category %d", i_c));
    }
  }
  
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
  
  return 0;
}

