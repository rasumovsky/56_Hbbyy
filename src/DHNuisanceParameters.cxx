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
#include "DHNP.h"

/**
   -----------------------------------------------------------------------------
   Get the intersection point for the graph.
*/
double getIntercept(TGraph *graph, double valueToIntercept) {
  
  // Loop over points in the graph to get search range:
  double rangeMin = 0.0; double rangeMax = 0.0;
  for (int i_p = 0; i_p < graph->GetN(); i_p++) {
    double xCurr; double yCurr;
    graph->GetPoint(i_p, xCurr, yCurr);
    if (i_p == 0) rangeMin = xCurr;
    if (i_p == (graph->GetN()-1)) rangeMax = xCurr;
  }
  
  // Bisection method to search for intercept:
  double precision = 0.0001;
  int nIterations = 0;
  int maxIterations = 30;
  double stepSize = (rangeMax - rangeMin) / 2.0;
  double currXValue = (rangeMax + rangeMin) / 2.0;
  double currYValue = graph->Eval(currXValue);
  while ((fabs(currYValue - valueToIntercept)/valueToIntercept) > precision && 
	 nIterations <= maxIterations) {
    
    currYValue = graph->Eval(currXValue);
    
    nIterations++;
    stepSize = 0.5 * stepSize;

    if (currYValue > valueToIntercept) currXValue -= stepSize;
    else currXValue += stepSize;
  }
  
  // Print error message and return bad value if convergence not achieved:
  if (nIterations == maxIterations) {
    std::cout << "DHNuisanceParameters: ERROR! Intercept not found."
	      << std::cout;
    return -999;
  }
  
  return currXValue;
}

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
    std::map<TString,DHNP*> parameterData; parameterData.clear();
    RooArgSet *nuisParameters = dhts->theModelConfig()->GetNuisanceParameters();
    TIterator *iterNuis = nuisParameters->createIterator();
    RooRealVar *currNuis = NULL;
    while ((currNuis = (RooRealVar*)iterNuis->Next())) {
      TString currNuisName = (TString)(currNuis->GetName());
      parameterData[currNuisName] = new DHNP(currNuisName);
    }
    
    //----------------------------------------//
    // First get NP values from profiling data:
    
    double profiledPOIVal = -999.0;
    // Mu = 0 fits:
    double nllMu0 = dhts->getFitNLL(m_dataToPlot, 0, true, profiledPOIVal);
    if (!dhts->fitsAllConverged()) m_allGoodFits = false;
    else {
      saveNuisData("Mu0", parameterData, 
		   dhts->theModelConfig()->GetNuisanceParameters());
    }
    
    // Mu = 1 fits:
    double nllMu1 = dhts->getFitNLL(m_dataToPlot, 1, true, profiledPOIVal);
    if (!dhts->fitsAllConverged()) m_allGoodFits = false;
    else {
      saveNuisData("Mu1", parameterData, 
		   dhts->theModelConfig()->GetNuisanceParameters());
    }

    // Mu free fits:
    double nllMuFree = dhts->getFitNLL(m_dataToPlot, 1, false, profiledPOIVal);
    if (!dhts->fitsAllConverged()) m_allGoodFits = false;
    else {
      saveNuisData("MuFree", parameterData, 
		   dhts->theModelConfig()->GetNuisanceParameters());
      
      // Also get total uncertainty on mu here.
    }
    
    

    
    
    
    
    
    //----------------------------------------//
    // Then get the impact of each NP on the signal strength:
    
    // Loop over the NP:
    
    // fit the NP to the best-fit value from mu-free fit:
    
    // float other NP:
    
    // Get the uncertainty on mu from this fit
    
    // Take the difference as the uncertainty on mu from this NP
    


    
    // Calculate the 95% CL value:
    dhts->calculateNewCL();
    
    // Check that the fit converges:
    if (dhts->fitsAllConverged()) {
      
      dhts->clearData();
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
