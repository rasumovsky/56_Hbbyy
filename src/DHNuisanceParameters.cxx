////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHNuisanceParameters.cxx                                            //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 23/02/2016                                                          //
//                                                                            //
//  Obtains the nuisance parameter pull and ranking plot.                     //
//                                                                            //
//  Macro options:                                                            //
//  - "New"        Calculate everything from scratch.                         //
//  - "FromFile"   Load CL values from file.                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "DHTestStat.h"
#include "DHToyAnalysis.h"

/**
   -----------------------------------------------------------------------------
   The main method scans the 95% CL for various signal cross-sections.
   @param configFile - The analysis configuration file.
   @param options - Job options: "New","FromFile","toy","asymptotic","NEvents"
   @param resMass - The resonance mass.
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 4) {
    std::cout << "\nUsage: " << argv[0]
	      << " <configFile> <options> <resMass>" << std::endl;
    exit(0);
  }
  
  TString configFile = argv[1];
  TString options = argv[2];
  int resonanceMass = atoi(argv[3]);
  
  // Load the analysis configuration file:
  Config *config = new Config(configFile);
  TString jobName = config->getStr("JobName");
  TString anaType = config->getStr("AnalysisType");
  TString outputDir = Form("%s/%s/DHNuisanceParameters", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Create a nametag for results:
  TString tag = (anaType.EqualTo("Resonant")) ? 
    Form("ResonantMX%d", resonanceMass) : "NonResonant";
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();

  // Store (unranked) names of nuisance parameters:
  std::vector<TString> paramNames; paramNames.clear();
  
  // Store impact on mu:
  double impactRank[100] = {0};
  double impactRankErr[100] = {0.5};
  double points[100] = {0.0};
  double muLo[100] = {0.0};
  double muHi[100] = {0.0};
  
  // Store the NP info for the mu=0 fit:
  double valMu0[100] = {0.0};
  double errLoMu0[100] = {0.0};
  double errHiMu0[100] = {0.0};
  
  // Store the NP info for the mu=1 fit:
  double valMu1[100] = {0.0};
  double errLoMu1[100] = {0.0};
  double errHiMu1[100] = {0.0};
  
  // Store the NP info for the mu=free fit:
  double valMuFree[100] = {0.0};
  double errLoMuFree[100] = {0.0};
  double errHiMuFree[100] = {0.0};
  
  TString fileNameNP = Form("%s/nuisance_parameters_%s.txt", 
			    outputDir.Data(), tag.Data());
  
  // Keep track of whether all fits work:
  bool allGoodFits = true;
  
  //----------------------------------------//
  // Open nuisance parameter values from file:
  if (options.Contains("FromFile")) {
    
    // Open the saved CL values from asymptotics:
    
    ifstream inputFile(fileNameNP);
    if (inputFile.is_open()) {
      int i_n = 0; TString currName = "";
      while (inputFile >> currName >> points[i_n] >> muLo[i_n] >> muHi[i_n]
	     >> valMu0[i_n] >> errLoMu0[i_n] >> errHiMu0[i_n] >> valMu1[i_n]
	     >> errLoMu1[i_n] >> errHiMu1[i_n] >> valMuFree[i_n] 
	     >> errLoMuFree[i_n] >> errHiMuFree[i_n]) {
	paramNames.push_back(currName);
	i_n++;
      }
      inputFile.close();
    }
    else {
      std::cout << "DHNuisanceParameters: ERROR opnening file " 
		<< fileNameNP << std::endl;
      exit(0);
    }
  }
  
  //----------------------------------------//
  // Calculate CL values from scratch and save them:
  else {
    // Open the workspace:
    TString originFile = Form("%s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
			      (config->getStr("MasterOutput")).Data(), 
			      jobName.Data(), anaType.Data());
    TFile wsFile(originFile, "read");
    RooWorkspace *workspace = (RooWorkspace*)wsFile.Get("combinedWS");
    
    // Instantiate the test statistic class for calculations and plots:
    DHTestStat *dhts = new DHTestStat(configFile, "new", workspace);
    
    // Get the name of the dataset for fitting:
    TString dataset = config->getStr("RankNPData");
    
    // Set the value of the variable to scan:
    dhts->setParam("xsec_SigBSM2H", config->getNum("RankNPCrossSection"), true);
    
    // If doing the resonant analysis, also set the res mass value:
    if (anaType.EqualTo("Resonant")) {
      dhts->setParam("mResonance", config->getNum("RankNPResonanceMass"), true);
    }
    
    // Keep nuisance parameter values after fitting (so that they can be saved):
    dhts->resetParamsAfterFit(false);
    
    // Create an (ordered) vector of nuisance parameter names:
    const RooArgSet *nuisParameters 
      = dhts->theModelConfig()->GetNuisanceParameters();
    TIterator *iterNuis = nuisParameters->createIterator();
    RooRealVar *currNuis = NULL;
    while ((currNuis = (RooRealVar*)iterNuis->Next())) {
      paramNames.push_back((TString)(currNuis->GetName()));
    }
    
    //----------------------------------------//
    // First get NP values from profiling data:
    std::cout << "DHNuisanceParameters: get NP values from profiling."
	      << std::endl;
    double profiledPOIVal = -999.0;
    // Mu = 0 fits:
    double nllMu0 = dhts->getFitNLL(dataset, 0, true, profiledPOIVal);
    for (int i_n = 0; i_n < (int)paramNames.size(); i_n++) {
      if (dhts->theWorkspace()->var(paramNames[i_n])) {
	valMu0[i_n] = dhts->theWorkspace()->var(paramNames[i_n])->getVal();
	errLoMu0[i_n] =dhts->theWorkspace()->var(paramNames[i_n])->getError();
	errHiMu0[i_n] =dhts->theWorkspace()->var(paramNames[i_n])->getError();
      }
      else {
	std::cout << "DHNuisanceParameters: parameter " << paramNames[i_n]
		  << " not found." << std::endl;
      }
    }
    
    // Mu = 1 fits:
    double nllMu1 = dhts->getFitNLL(dataset, 1, true, profiledPOIVal);
    for (int i_n = 0; i_n < (int)paramNames.size(); i_n++) {
      if (dhts->theWorkspace()->var(paramNames[i_n])) {
	valMu1[i_n] = dhts->theWorkspace()->var(paramNames[i_n])->getVal();
	errLoMu1[i_n] =dhts->theWorkspace()->var(paramNames[i_n])->getError();
	errHiMu1[i_n] =dhts->theWorkspace()->var(paramNames[i_n])->getError();
      }
      else {
	std::cout << "DHNuisanceParameters: parameter " << paramNames[i_n]
		  << " not found." << std::endl;
      }
    }
    
    // Mu free fits:
    double nllMuFree = dhts->getFitNLL(dataset, 1, false, profiledPOIVal);
    for (int i_n = 0; i_n < (int)paramNames.size(); i_n++) {
      if (dhts->theWorkspace()->var(paramNames[i_n])) {
	valMu1[i_n] = dhts->theWorkspace()->var(paramNames[i_n])->getVal();
	errLoMu1[i_n] =dhts->theWorkspace()->var(paramNames[i_n])->getError();
	errHiMu1[i_n] =dhts->theWorkspace()->var(paramNames[i_n])->getError();
      }
      else {
	std::cout << "DHNuisanceParameters: parameter " << paramNames[i_n]
		  << " not found." << std::endl;
      }
    }
    
    //----------------------------------------//
    // Then get the impact of each NP on the signal strength:
    std::cout << "DHNuisanceParameters: get NP impact on signal strength."
	      << std::endl;
    // First get the full uncertainty on mu:
    std::vector<TString> nullSet; nullSet.clear();
    std::map<int,double> mlScan
      = dhts->scanNLL("mlScan", dataset, config->getStr("RankNPPoI"), nullSet);
    
    // Loop over the NP:
    for (int i_n = 0; i_n < (int)paramNames.size(); i_n++) {
      
      // Set the current parameter constant in the fit:
      std::vector<TString> currParamsToFix;
      currParamsToFix.clear();
      currParamsToFix.push_back(paramNames[i_n]);
      
      // Get the uncertainty on mu from this fit:
      std::map<int,double> currScan
	= dhts->scanNLL("currScan", dataset, config->getStr("RankNPPoI"),
			currParamsToFix);
      
      // Take the difference with usual error as uncertainty on mu from this NP:
      muLo[i_n] = -1.0 * fabs(mlScan[-1] - currScan[-1]);
      muHi[i_n] = fabs(currScan[1] - mlScan[1]);
      points[i_n] = i_n;
    }
    allGoodFits = dhts->fitsAllConverged();

    //----------------------------------------//
    // Write nuisance parameter values to file:
    ofstream outFile(fileNameNP);
    for (int i_n = 0; i_n < (int)paramNames.size(); i_n++) {
      outFile << paramNames[i_n] << " " << points[i_n] << " " << muLo[i_n] 
	      << " " << muHi[i_n] << " " << valMu0[i_n] << " " << errLoMu0[i_n] 
	      << " " << errHiMu0[i_n] << " " << valMu1[i_n] << " " 
	      << errLoMu1[i_n] << " " << errHiMu1[i_n] << " " << valMuFree[i_n] 
	      << " " << errLoMuFree[i_n] << " " << errHiMuFree[i_n] 
	      << std::endl;
    }
    outFile.close();
    
    delete dhts;
    delete workspace;
    wsFile.Close();
  }

  //----------------------------------------//
  // Rank the nuisance parameters by impact on mu:
  std::cout << "DHNuisanceParameters: Rank NP according to mu impact."
	      << std::endl;
  std::vector<TString> sortedParamNames; sortedParamNames.clear();
  std::vector<double> sortedParamValues; sortedParamValues.clear();
  for (int i_n = 0; i_n < (int)paramNames.size(); i_n++) {
    double currMaxErr = TMath::Max(fabs(muLo[i_n]), fabs(muHi[i_n]));
    if (i_n == 0) {
      sortedParamNames.push_back(paramNames[i_n]);
      sortedParamValues.push_back(currMaxErr);
    }
    else {
      std::vector<TString>::iterator iterName = sortedParamNames.begin();
      std::vector<double>::iterator iterVal= sortedParamValues.begin();
      bool inserted = false;
      while (iterName != sortedParamNames.end()) {
	if (currMaxErr >= *iterVal) {
	  sortedParamNames.insert(iterName, paramNames[i_n]);
	  sortedParamValues.insert(iterVal, currMaxErr);
	  inserted = true;
	  break;
	}
	else {
	  iterName++;
	  iterVal++;
	}
      }
      if (!inserted) {
	sortedParamNames.push_back(paramNames[i_n]);
	sortedParamValues.push_back(currMaxErr);
      }
    }
  }
  
  // Finally, assign the rank of each NP:
  for (int i_n = 0; i_n < (int)paramNames.size(); i_n++) {
    for (int i_r = 0; i_r < (int)sortedParamNames.size(); i_n++) {
      if ((sortedParamNames[i_r]).EqualTo(paramNames[i_n])) {
	impactRank[i_n] = (int)(sortedParamNames.size()) - i_r;
      }
    }
  }
  
  //----------------------------------------//
  // Plot the results:
  std::cout << "DHNuisanceParameters: plot NP ranking results."
	    << std::endl;
  
  // Start plotting:
  TCanvas *can = new TCanvas("can","can", 500, 800);
  gStyle->SetPadLeftMargin(0.2);
  can->cd();
  
  // Stupid histogram:
  TH1F *hTemp = new TH1F("temp", "temp", (int)sortedParamNames.size()+5, 0,
			 (int)sortedParamNames.size()+5);
  for (int i_h = 0; i_h < (int)sortedParamNames.size()+5; i_h++) {
    hTemp->SetBinContent(i_h+1, 0);
    hTemp->GetXaxis()->SetBinLabel(i_h+1, sortedParamNames[i_h]);
  }
  hTemp->GetYaxis()->SetTitle("#delta#hat{#mu}");
  hTemp->GetYaxis()->SetRangeUser(-0.25, 2.5);
  hTemp->Draw("hbar");
  
  // Fill graphs with fit data:
  double zeros[100] = {0.0};
  TGraphAsymmErrors *gMuRank
    = new TGraphAsymmErrors((int)paramNames.size(), zeros, impactRank, 
			    muLo, muHi, impactRankErr, impactRankErr);
  TGraphAsymmErrors *gFitMu0
    = new TGraphAsymmErrors((int)paramNames.size(), valMu0, impactRank,
			    errLoMu0, errHiMu0, zeros, zeros);
  TGraphAsymmErrors *gFitMu1
    = new TGraphAsymmErrors((int)paramNames.size(), valMu1, impactRank,
			    errLoMu1, errHiMu1, zeros, zeros);
  TGraphAsymmErrors *gFitMuFree
    = new TGraphAsymmErrors((int)paramNames.size(), valMuFree, impactRank,
			    errLoMuFree, errHiMuFree, zeros, zeros);
  
  //
  gMuRank->SetFillColor(kOrange+1);
  gMuRank->Draw("a2");
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.44, 0.48, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.56, 0.48, config->getStr("ATLASLabel"));
  t.DrawLatex(0.44, 0.42, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			       (config->getNum("AnalysisLuminosity")/1000.0)));
  
  // Legend:
  TLegend leg(0.7,0.85,0.9,0.92);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(gMuRank, "Post-fit #hat{#mu} impact", "F");
  leg.AddEntry(gFitMu0, "Pull, #mu=0", "LEP");
  leg.AddEntry(gFitMu1, "Pull, #mu=1", "LEP");
  leg.AddEntry(gFitMuFree, "#Pull, #hat{#mu}", "LEP");
  leg.Draw("SAME");
  
  // 95% CL Line:
  /*
  TLine *line = new TLine();
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  if (options.Contains("toy") || options.Contains("both")) {
    line->DrawLine(gCLExp_toy->GetXaxis()->GetXmin(), 0.95,
		   gCLExp_toy->GetXaxis()->GetXmax(), 0.95);
  }
  else {
    line->DrawLine(gCLExp_asym->GetXaxis()->GetXmin(), 0.95,
		   gCLExp_asym->GetXaxis()->GetXmax(), 0.95);
  }
  */
    
  // Print the canvas:
  can->Print(Form("%s/plot_parameters_%s.eps", outputDir.Data(), tag.Data()));
  can->Print(Form("%s/plot_parameters_%s.C", outputDir.Data(), tag.Data()));
  
  // Save graphs to a ROOT file:
  TFile rootFileOutput(Form("%s/nuisance_parameter_graphs_%s.root",
			    outputDir.Data(), tag.Data()));
  gMuRank->Write();
  gFitMu0->Write();
  gFitMu1->Write();
  gFitMuFree->Write();
  rootFileOutput.Close();
  
  // Print a warning if some fits failed:
  if (!allGoodFits) {
    std::cout << "\nDHNuisanceParameters: ERROR! Fits failed somewhere..."
	      << "\n\t...queue Charlie Brown music" << std::endl;
  }
  else std::cout << "DHNuisanceParameters: Finished!" << std::endl;
  
  //delete line;
  delete can;
  delete config;
  return 0;
}
