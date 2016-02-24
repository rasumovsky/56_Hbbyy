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
  
  // Define the input file, then make a local copy (for remote jobs):
  TString originFile = Form("%s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
			    (config->getStr("MasterOutput")).Data(), 
			    jobName.Data(), anaType.Data());
  
  // Create a nametag for results:
  TString tag = (anaType.EqualTo("Resonant")) ? 
    Form("ResonantMX%d", resonanceMass) : "NonResonant";
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Open the workspace:
  TFile wsFile(originFile, "read");
  RooWorkspace *workspace = (RooWorkspace*)wsFile.Get("combinedWS");
  
  // Instantiate the test statistic class for calculations and plots:
  DHTestStat *dhts = new DHTestStat(configFile, "new", workspace);
    
  //----------------------------------------//
  // Open nuisance parameter values from file:
  if (options.Contains("FromFile")) {
    
    // Open the saved CL values from asymptotics:
    TString inputNPName = Form("%s/nuisanceParameterValues_%s.txt",
				 outputDir.Data(), tag.Data());
    ifstream inputFile(inputNPName);
    if (inputFile.is_open()) {
      while (inputFile >> varValues_asym[n_asym] >> CLObs_asym[n_asym] 
	     >> CLExp_asym[n_asym] >> CLExp_asym_p2[n_asym] 
	       >> CLExp_asym_p1[n_asym] >> CLExp_asym_n1[n_asym] 
	       >> CLExp_asym_n2[n_asym]) {
	  n_asym++;
	}
      }
      else {
	std::cout << "DHNuisanceParameters: ERROR opnening file " 
		  << inputNPName << std::endl;
	exit(0);
      }
      
      inputFile_asym.close();
    }
  }
  
  //----------------------------------------//
  // Calculate CL values from scratch and save them:
  else {
    // Save values for plotting again!
    ofstream outFile_asym.open(Form("%s/scan_CLvalues_asym_%s.txt", 
				    outputDir.Data(), tag.Data()));
          
    // Set the value of the variable to scan:
    dhts->setParam(config->getStr("CLScanVar"), crossSection, true);
    
    // If doing the resonant analysis, also set the res mass value:
    if (anaType.EqualTo("Resonant")) {
      dhts->setParam("mResonance", resonanceMass, true);
    }
  
    // Keep nuisance parameter values after fitting (so that they can be saved):
    dhts->resetParamsAfterFit(true);
    
    //----------------------------------------//
    // Storage for nuisance parameter data:
    std::vector<TString> paramNames; paramNames.clear();
    RooArgSet *nuisParameters = dhts->theModelConfig()->GetNuisanceParameters();
    TIterator *iterNuis = nuisParameters->createIterator();
    RooRealVar *currNuis = NULL;
    while ((currNuis = (RooRealVar*)iterNuis->Next())) {
      paramNames.push_back((TString)(currNuis->GetName()));
    }
    
    double points[100] = {0.0};
    double muLo[100] = {0.0};
    double muHi[100] = {0.0};
    
    double valMu0[100] = {0.0};
    double errLoMu0[100] = {0.0};
    double errHiMu0[100] = {0.0};
    
    double valMu1[100] = {0.0};
    double errLoMu1[100] = {0.0};
    double errHiMu1[100] = {0.0};
    
    double valMuFree[100] = {0.0};
    double errLoMuFree[100] = {0.0};
    double errHiMuFree[100] = {0.0};
    
    //----------------------------------------//
    // First get NP values from profiling data:
    
    double profiledPOIVal = -999.0;
    // Mu = 0 fits:
    double nllMu0 = dhts->getFitNLL(m_dataToPlot, 0, true, profiledPOIVal);
    if (!dhts->fitsAllConverged()) m_allGoodFits = false;
    else {
      RooArgSet *setMu0 = dhts->theModelConfig()->GetNuisanceParameters();
    }
    
    // Mu = 1 fits:
    double nllMu1 = dhts->getFitNLL(m_dataToPlot, 1, true, profiledPOIVal);
    if (!dhts->fitsAllConverged()) m_allGoodFits = false;
    else {
      RooArgSet *setMu1 = dhts->theModelConfig()->GetNuisanceParameters();
    }

    // Mu free fits:
    double nllMuFree = dhts->getFitNLL(m_dataToPlot, 1, false, profiledPOIVal);
    if (!dhts->fitsAllConverged()) m_allGoodFits = false;
    else {
      RooArgSet *setMuFree = dhts->theModelConfig()->GetNuisanceParameters();
    }
        
    
    
    //----------------------------------------//
    // Then get the impact of each NP on the signal strength:
    
    // First get the full uncertainty on mu:
    std::vector<TString> nullSet; nullSet.clear();
    std::map<int,double> mlScan
      = dhts->scanNLL("mlScan", config->getStr("DataToScan"),
		      config->getStr("ParamToScan"), nullSet);
    
    // Loop over the NP:
    for (int i_n = 0; i_n < (int)paramNames.size(); i_n++) {
      
      // Set the current parameter constant in the fit:
      std::vector<TString> currParamsToFix;
      currParamToFix.clear();
      currParamsToFix.push_back(paramNames[i_n]);
      
      // Get the uncertainty on mu from this fit:
      std::map<int,double> currScan
	= dhts->scanNLL("currScan", config->getStr("DataToScan"),
			config->getStr("ParamToScan"), currParamsToFix);
      
      // Take the difference with usual error as uncertainty on mu from this NP:
      muLo[i_n] = -1.0 * fabs(mlScan[-1] - currScan[-1]);
      muHi[i_n] = fabs(currScan[1] - mlScan[1]);
      points[i_n] = i_n;
    }

    
        
    // Close the files that save CL data:
    outFile_asym.close();
  }
  
  //----------------------------------------//
  // Plot the results:
    
  // Start plotting:
  TCanvas *can = new TCanvas("can","can");
  can->cd();
  
  // Legend:
  TLegend leg(0.54,0.20,0.79,0.38);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(gCLExp_toy,"Exp. limit","l");
  leg.AddEntry(gCLExp_toy_1s,"Exp. limit #pm1#sigma_{exp}","F");
  leg.AddEntry(gCLExp_toy_2s,"Exp. limit #pm2#sigma_{exp}","F");
  leg.Draw("SAME");
  
  // 95% CL Line
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
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.54, 0.48, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.66, 0.48, config->getStr("ATLASLabel"));
  t.DrawLatex(0.54, 0.42, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			       (config->getNum("AnalysisLuminosity")/1000.0)));
  
  // Print the canvas:
  can->Print(Form("%s/scan95CL_asymptotic_%s.eps",outputDir.Data(),tag.Data()));
  can->Print(Form("%s/scan95CL_asymptotic_%s.C",outputDir.Data(),tag.Data()));
  
  // Delete pointers, close files, return:
  std::cout << "DHNuisanceParameters: Finished!" << std::endl;
  delete line;
  delete can;
  delete gCLObs_asym;
  delete gCLExp_asym;
  delete gCLObs_toy;
  delete gCLExp_toy;
  delete gCLExp_asym_2s;
  delete gCLExp_asym_1s;
  delete gCLExp_toy_2s;
  delete gCLExp_toy_1s;
  delete workspace;
  delete config;
  wsFile.Close();
  return 0;
}
