////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHCLScan.cxx                                                        //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 27/01/2016                                                          //
//                                                                            //
//  Performs a scan of the 95% CL for the non-resonant Di-Higgs analysis as a //
//  function of signal cross-section.                                         //
//                                                                            //
//  Macro options:                                                            //
//  - "New"        Calculate everything from scratch.                         //
//  - "FromFile"   Load CL values from file.                                  //
//  - "toy"        Get limits from toy MC jobs.                               //
//  - "asymptotic" Get asymptotic limits.                                     //
//  - "both"       Compares toy and asymptotic limits.                        //
//  - "NEvents"    Calculate limits in terms of number of events.             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "DHTestStat.h"
#include "DHToyAnalysis.h"

/**
   -----------------------------------------------------------------------------
   Get the intersection point for the graph.
   @param graph - The graph for which intercept points will be found.
   @param valueToIntercept - The y-value to intercept.
   @return - The x-value of the intercept.
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
    std::cout << "DHCLScan: ERROR! Intercept not found." << std::cout;
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
  TString outputDir = Form("%s/%s/DHCLScan", 
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
  if (options.Contains("NEvents")) tag.Append("_NEvent");
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Open the workspace:
  TFile wsFile(originFile, "read");
  RooWorkspace *workspace = (RooWorkspace*)wsFile.Get("combinedWS");
  
  // Instantiate the test statistic class for calculations and plots:
  DHTestStat *dhts = new DHTestStat(configFile, "new", workspace);
  
  // Arrays to store band information:
  double varValues_asym[100] = {0};
  double CLObs_asym[100] = {0};
  double CLExp_asym_p2[100] = {0};
  double CLExp_asym_p1[100] = {0};
  double CLExp_asym[100] = {0};
  double CLExp_asym_n1[100] = {0};
  double CLExp_asym_n2[100] = {0};
  double varValues_toy[100] = {0};  
  double CLObs_toy[100] = {0};
  double CLExp_toy_p2[100] = {0};  
  double CLExp_toy_p1[100] = {0};
  double CLExp_toy[100] = {0};
  double CLExp_toy_n1[100] = {0};
  double CLExp_toy_n2[100] = {0};
  
  double qMuObs_toy[100] = {0.0};
  double qMuExp_toy[100] = {0.0};
  double qMuExp_toy_p2[100] = {0.0};
  double qMuExp_toy_p1[100] = {0.0};
  double qMuExp_toy_n1[100] = {0.0};
  double qMuExp_toy_n2[100] = {0.0};
  
  int n_asym = 0;
  int n_toy = 0;
  
  // Scan information:
  std::vector<double> scanValues = config->getNumV("CLScanValues");
  
  //----------------------------------------//
  // Open CL values from file:
  if (options.Contains("FromFile")) {

    // Open the saved CL values from asymptotics:
    TString inputAsymName = Form("%s/scan_CLvalues_asym_%s.txt",
				 outputDir.Data(), tag.Data());
    if (options.Contains("asymptotic") || options.Contains("both")) {
      ifstream inputFile_asym(inputAsymName);
      if (inputFile_asym.is_open()) {
	while (inputFile_asym >> varValues_asym[n_asym] >> CLObs_asym[n_asym] 
	       >> CLExp_asym[n_asym] >> CLExp_asym_p2[n_asym] 
	       >> CLExp_asym_p1[n_asym] >> CLExp_asym_n1[n_asym] 
	       >> CLExp_asym_n2[n_asym]) {
	  n_asym++;
	}
      }
      else {
	std::cout << "DHCLScan: ERROR opnening asymptotic file " 
		  << inputAsymName << std::endl;
	exit(0);
      }
      
      inputFile_asym.close();
    }
    //bool skipped = false;
    // Open the saved CL values from toys:
    TString inputToyName = Form("%s/scan_CLvalues_toy_%s.txt",
				outputDir.Data(), tag.Data());
    if (options.Contains("toy") || options.Contains("both")) {
      ifstream inputFile_toy(inputToyName);
      if (inputFile_toy.is_open()) {
	while (inputFile_toy >> varValues_toy[n_toy] >> CLObs_toy[n_toy] 
	       >> CLExp_toy[n_toy] >> CLExp_toy_p2[n_toy] 
	       >> CLExp_toy_p1[n_toy] >> CLExp_toy_n1[n_toy] 
	       >> CLExp_toy_n2[n_toy]) {
	  n_toy++;
	}
      }
      else {
	std::cout << "DHCLScan: ERROR opnening toy file " << inputToyName
		  << std::endl;
	exit(0);
      }
      
      inputFile_toy.close();
    }
  }
  
  //----------------------------------------//
  // Calculate CL values from scratch and save them:
  else {
    // Save values for plotting again!
    ofstream outFile_asym;
    ofstream outFile_toy;
    if (options.Contains("asymptotic") || options.Contains("both")) {
      outFile_asym.open(Form("%s/scan_CLvalues_asym_%s.txt", 
			     outputDir.Data(), tag.Data()));
    }
    if (options.Contains("toy") || options.Contains("both")) {
      outFile_toy.open(Form("%s/scan_CLvalues_toy_%s.txt",
			    outputDir.Data(), tag.Data()));
    }
    
    // Loop over cross-section:
    for (int i_x = 0; i_x < (int)scanValues.size(); i_x++) {
      double crossSection = scanValues[i_x];
      std::cout << "DHCLScan: cross-section = " << crossSection << std::endl;
      
      // Set the value of the variable to scan:
      dhts->setParam(config->getStr("CLScanVar"), crossSection, true);
      
      // If doing the resonant analysis, also set the res mass value:
      if (anaType.EqualTo("Resonant")) {
	dhts->setParam("mResonance", resonanceMass, true);
      }
      
      // Asymptotic calculation of CL:
      if (options.Contains("asymptotic") || options.Contains("both")) {
	std::cout << "DHCLScan: calculating asymptotic CL" << std::endl;
	
	// Calculate the 95% CL value:
	dhts->calculateNewCL();
	
	// Check that the fit converges:
	if (dhts->fitsAllConverged()) {
	  if (options.Contains("NEvents")) {
	    varValues_asym[n_asym]
	      = (dhts->theWorkspace())->function("nRaw_SigBSM2H_bb")->getVal();
	  }
	  else varValues_asym[n_asym] = crossSection / 1000.0;
	  
	  CLObs_asym[n_asym] = dhts->accessValue("CL", true, 0);
	  CLExp_asym[n_asym] = dhts->accessValue("CL", false, 0);
	  CLExp_asym_p2[n_asym] = dhts->accessValue("CL", false, 2);
	  CLExp_asym_p1[n_asym] = dhts->accessValue("CL", false, 1);
	  CLExp_asym_n1[n_asym] = dhts->accessValue("CL", false, -1);
	  CLExp_asym_n2[n_asym] = dhts->accessValue("CL", false, -2);
	  
	  // Write CL values to file:
	  outFile_asym << varValues_asym[n_asym] << " " 
		       << CLObs_asym[n_asym] << " " 
		       << CLExp_asym[n_asym] << " " 
		       << CLExp_asym_p2[n_asym] << " " 
		       << CLExp_asym_p1[n_asym] << " " 
		       << CLExp_asym_n1[n_asym] << " " 
		       << CLExp_asym_n2[n_asym] << std::endl;
	}
	n_asym++;
	dhts->clearData();
      }
      
      // Pseudo-experiment calculation of CL:
      if (options.Contains("toy") || options.Contains("both")) {
	std::cout << "DHCLScan: calculating toy CL" << std::endl;
	
	// Calculate the observed qMu values:
	double obsPoI = 0.0;
	double nllObsMu1
	  = dhts->getFitNLL("obsData", 1, true, obsPoI, false);
	double nllObsMuFree
	  = dhts->getFitNLL("obsData", 1, false, obsPoI, false);
	qMuObs_toy[n_toy]
	  = dhts->getQMuFromNLL(nllObsMu1, nllObsMuFree, obsPoI, 1);
	std::cout << "nllObsMu1=" << nllObsMu1 
		  << ", nllObsMuFree=" << nllObsMuFree
		  << ", qMuObs=" << qMuObs_toy[n_toy] << std::endl;
	if (options.Contains("NEvents")) {
	  varValues_toy[n_toy]
	    = (dhts->theWorkspace())->function("nRaw_SigBSM2H_bb")->getVal();
	}
	else varValues_toy[n_toy] = crossSection / 1000.0;
	
	n_toy++;
      }
    }
    delete dhts;
    
    // Then specify file in DHToyAnalysis to load:
    if (options.Contains("toy") || options.Contains("both")) {
      for (int i_t = 0; i_t < n_toy; i_t++) {
	
	// Load the tool to analyze toys.
	// NOTE: this was moved outside the loop above because of interference
	// between the DHToyAnalysis class and DHTestStat, which is called
	// in DHToyAnalysis...
	TString toyScanOption = Form("CLScan%d",i_t);
	if (i_t == 1) toyScanOption.Append("_ForcePlot");
	DHToyAnalysis *dhta
	  = new DHToyAnalysis(configFile, toyScanOption, resonanceMass);
	if (!(dhta->areInputFilesOK())) {
	  std::cout << "DHCLScan: ERROR with toy scan option " << toyScanOption
		    << std::endl;
	  continue;
	}
	
	// Calculate the expected qmu:
	qMuExp_toy[i_t] = dhta->calculateBkgQMuForN(0);
	qMuExp_toy_p2[i_t] = dhta->calculateBkgQMuForN(2.0);
	qMuExp_toy_p1[i_t] = dhta->calculateBkgQMuForN(1.0);
	qMuExp_toy_n1[i_t] = dhta->calculateBkgQMuForN(-1.0);
	qMuExp_toy_n2[i_t] = dhta->calculateBkgQMuForN(-2.0);
	
	CLObs_toy[i_t] = dhta->calculateCLFromToy(qMuObs_toy[i_t]);
	CLExp_toy[i_t] = dhta->calculateCLFromToy(qMuExp_toy[i_t]);
	CLExp_toy_p2[i_t] = dhta->calculateCLFromToy(qMuExp_toy_p2[i_t]);
	CLExp_toy_p1[i_t] = dhta->calculateCLFromToy(qMuExp_toy_p1[i_t]);
	CLExp_toy_n1[i_t] = dhta->calculateCLFromToy(qMuExp_toy_n1[i_t]);
	CLExp_toy_n2[i_t] = dhta->calculateCLFromToy(qMuExp_toy_n2[i_t]);
	
	// Write CL values to file:
	outFile_toy << varValues_toy[i_t] << " " 
		    << CLObs_toy[i_t] << " " 
		    << CLExp_toy[i_t] << " "
		    << CLExp_toy_p2[i_t] << " " 
		    << CLExp_toy_p1[i_t] << " " 
		    << CLExp_toy_n1[i_t] << " " 
		    << CLExp_toy_n2[i_t] << std::endl;
	
	delete dhta;
      }
    }
    
    // Close the files that save CL data:
    outFile_asym.close();
    outFile_toy.close();
  }
  
  //----------------------------------------//
  // Plot the results:
  
  double errExp_asym_p2[100] = {0};
  double errExp_asym_p1[100] = {0};
  double errExp_asym_n1[100] = {0};
  double errExp_asym_n2[100] = {0};
  for (int i_a = 0; i_a < n_asym; i_a++) {
    errExp_asym_p2[i_a] = fabs(CLExp_asym_p2[i_a] - CLExp_asym[i_a]);
    errExp_asym_p1[i_a] = fabs(CLExp_asym_p1[i_a] - CLExp_asym[i_a]);
    errExp_asym_n1[i_a] = fabs(CLExp_asym_n1[i_a] - CLExp_asym[i_a]);
    errExp_asym_n2[i_a] = fabs(CLExp_asym_n2[i_a] - CLExp_asym[i_a]);
  }
  
  double errExp_toy_p2[100] = {0};  
  double errExp_toy_p1[100] = {0};
  double errExp_toy_n1[100] = {0};
  double errExp_toy_n2[100] = {0};
  for (int i_t = 0; i_t < n_toy; i_t++) {
    errExp_toy_p2[i_t] = fabs(CLExp_toy_p2[i_t] - CLExp_toy[i_t]);
    errExp_toy_p1[i_t] = fabs(CLExp_toy_p1[i_t] - CLExp_toy[i_t]);
    errExp_toy_n1[i_t] = fabs(CLExp_toy_n1[i_t] - CLExp_toy[i_t]);
    errExp_toy_n2[i_t] = fabs(CLExp_toy_n2[i_t] - CLExp_toy[i_t]);
  }
  
  // Median expected and observed results:
  TGraph *gCLObs_asym = new TGraph(n_asym, varValues_asym, CLObs_asym);
  TGraph *gCLExp_asym = new TGraph(n_asym, varValues_asym, CLExp_asym);
  TGraph *gCLExp_asym_p1 = new TGraph(n_asym, varValues_asym, CLExp_asym_p1);
  TGraph *gCLExp_asym_p2 = new TGraph(n_asym, varValues_asym, CLExp_asym_p2);
  TGraph *gCLExp_asym_n1 = new TGraph(n_asym, varValues_asym, CLExp_asym_n1);
  TGraph *gCLExp_asym_n2 = new TGraph(n_asym, varValues_asym, CLExp_asym_n2);
  
  TGraph *gCLObs_toy = new TGraph(n_toy, varValues_toy, CLObs_toy);
  TGraph *gCLExp_toy = new TGraph(n_toy, varValues_toy, CLExp_toy);
  TGraph *gCLExp_toy_p1 = new TGraph(n_toy, varValues_toy, CLExp_toy_p1);
  TGraph *gCLExp_toy_p2 = new TGraph(n_toy, varValues_toy, CLExp_toy_p2);
  TGraph *gCLExp_toy_n1 = new TGraph(n_toy, varValues_toy, CLExp_toy_n1);
  TGraph *gCLExp_toy_n2 = new TGraph(n_toy, varValues_toy, CLExp_toy_n2);
  
  // Also plot the bands:
  TGraphAsymmErrors *gCLExp_asym_2s 
    = new TGraphAsymmErrors(n_asym, varValues_asym, CLExp_asym, 0, 0, 
			    errExp_asym_n2, errExp_asym_p2);
  TGraphAsymmErrors *gCLExp_asym_1s
    = new TGraphAsymmErrors(n_asym, varValues_asym, CLExp_asym, 0, 0, 
			    errExp_asym_n1, errExp_asym_p1);
  TGraphAsymmErrors *gCLExp_toy_2s 
    = new TGraphAsymmErrors(n_toy, varValues_toy, CLExp_toy, 0, 0, 
			    errExp_toy_n2, errExp_toy_p2);
  TGraphAsymmErrors *gCLExp_toy_1s
    = new TGraphAsymmErrors(n_toy, varValues_toy, CLExp_toy, 0, 0, 
			    errExp_toy_n1, errExp_toy_p1);
  
  // Start plotting:
  TCanvas *can = new TCanvas("can","can");
  can->cd();
  
  // Toy graph formatting:
  if (options.Contains("NEvents")) {
    gCLExp_toy->GetXaxis()->SetTitle("X#rightarrowhh Events");
    gCLObs_toy->GetXaxis()->SetTitle("X#rightarrowhh Events");
  }
  else {
    gCLExp_toy->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
    gCLObs_toy->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
  }
  gCLExp_toy->GetYaxis()->SetTitle("CL_{s} value");
  gCLObs_toy->GetYaxis()->SetTitle("CL_{s} value");
  gCLExp_toy->SetLineColor(kBlack);
  gCLObs_toy->SetLineColor(kBlack);
  gCLExp_toy->SetLineStyle(2);
  gCLObs_toy->SetLineStyle(1);
  gCLExp_toy->SetLineWidth(2);
  gCLObs_toy->SetLineWidth(2);
  gCLExp_toy->GetYaxis()->SetRangeUser(0.0, 1.0);
  gCLExp_toy_2s->SetFillColor(kYellow);
  gCLExp_toy_1s->SetFillColor(kGreen);

  // Asymptotic graph formatting:
  if (options.Contains("NEvents")) {
    gCLExp_asym->GetXaxis()->SetTitle("X#rightarrowhh Events");
    gCLObs_asym->GetXaxis()->SetTitle("X#rightarrowhh Events");
  }
  else {
    gCLExp_asym->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
    gCLObs_asym->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
  }
  gCLExp_asym->GetYaxis()->SetTitle("CL_{s} value");
  gCLObs_asym->GetYaxis()->SetTitle("CL_{s} value");
  gCLExp_asym->SetLineColor(kBlack);
  gCLObs_asym->SetLineColor(kBlack);
  gCLExp_asym->SetLineStyle(2);
  gCLObs_asym->SetLineStyle(1);
  gCLExp_asym->SetLineWidth(2);
  gCLObs_asym->SetLineWidth(2);
  gCLExp_asym->GetYaxis()->SetRangeUser(0.0, 1.0);
  gCLExp_asym_2s->SetFillColor(kYellow);
  gCLExp_asym_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.64,0.20,0.89,0.38);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!config->getBool("DoBlind")) leg.AddEntry(gCLObs_toy,"Obs. CL_{s}","l");
  leg.AddEntry(gCLExp_toy,"Exp. CL_{s}","l");
  leg.AddEntry(gCLExp_toy_1s,"Exp. CL_{s} #pm1#sigma_{exp}","F");
  leg.AddEntry(gCLExp_toy_2s,"Exp. CL_{s} #pm2#sigma_{exp}","F");
  
  // Plotting options:
  if (options.Contains("toy")) {
    gCLExp_toy->Draw("AL");
    gCLExp_toy_2s->Draw("3same");
    gCLExp_toy_1s->Draw("3same");
    gCLExp_toy->Draw("LSAME");
    if (!config->getBool("DoBlind")) gCLObs_toy->Draw("LSAME");
  }
  else if (options.Contains("asymptotic")) {
    gCLExp_asym->Draw("AL");
    gCLExp_asym_2s->Draw("3same");
    gCLExp_asym_1s->Draw("3same");
    gCLExp_asym->Draw("LSAME");
    if (!config->getBool("DoBlind")) gCLObs_asym->Draw("LSAME");
  }
  else if (options.Contains("both")) {
    gCLExp_toy->Draw("AL");
    gCLExp_toy_2s->Draw("3same");
    gCLExp_toy_1s->Draw("3same");
    gCLExp_toy->Draw("LSAME");
    if (!config->getBool("DoBlind")) gCLObs_toy->Draw("LSAME");
    
    gCLExp_asym->SetLineColor(kRed+1);
    gCLObs_asym->SetLineColor(kRed+1);
    gCLExp_asym->SetLineStyle(2);
    gCLObs_asym->SetLineStyle(2);
    
    gCLExp_asym->Draw("LSAME");
    if (!config->getBool("DoBlind")) {
      gCLObs_asym->Draw("LSAME");
      leg.AddEntry(gCLObs_asym,"Obs. limit (asymptotic)","l");
    }
    leg.AddEntry(gCLExp_asym,"Exp. limit (asymptotic)","l");
  }
  gPad->RedrawAxis();
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
  t.DrawLatex(0.64, 0.48, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.76, 0.48, config->getStr("ATLASLabel"));
  t.DrawLatex(0.64, 0.42, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			       (config->getNum("AnalysisLuminosity")/1000.0)));
  
  // Print the canvas:
  if (options.Contains("asymptotic")) {
    can->Print(Form("%s/scan95CL_asymptotic_%s.eps",
		    outputDir.Data(), tag.Data()));
  }
  else if (options.Contains("toy")) {
    can->Print(Form("%s/scan95CL_toy_%s.eps",outputDir.Data(), tag.Data()));
  }
  else {
    can->Print(Form("%s/scan95CL_comparison_%s.eps",
		    outputDir.Data(), tag.Data()));
  }
  
  // Then save the intercepts to file:
  if (options.Contains("asymptotic")) {
    std::ofstream limitOutput(Form("%s/limits_asymptotic_%s.txt",
				   outputDir.Data(), tag.Data()));
    // Get the actual limit values:
    double observedCL_asym = getIntercept(gCLObs_asym, 0.95);
    double expectedCL_asym = getIntercept(gCLExp_asym, 0.95);
    double expectedCL_p1_asym = getIntercept(gCLExp_asym_p1, 0.95);
    double expectedCL_p2_asym = getIntercept(gCLExp_asym_p2, 0.95);
    double expectedCL_n1_asym = getIntercept(gCLExp_asym_n1, 0.95);
    double expectedCL_n2_asym = getIntercept(gCLExp_asym_n2, 0.95);  
    limitOutput << "observedCL " << observedCL_asym << std::endl;
    limitOutput << "expectedCL " << expectedCL_asym << std::endl;
    limitOutput << "expectedCL_p1 " << expectedCL_p1_asym << std::endl;
    limitOutput << "expectedCL_p2 " << expectedCL_p2_asym << std::endl;
    limitOutput << "expectedCL_n1 " << expectedCL_n1_asym << std::endl;
    limitOutput << "expectedCL_n2 " << expectedCL_n2_asym << std::endl;
    limitOutput.close();

    // Then print to screen:
    std::cout << "\nDHCLScan: Results" << std::endl;
    std::cout << "observedCL " << observedCL_asym << std::endl;
    std::cout << "expectedCL " << expectedCL_asym << std::endl;
    std::cout << "expectedCL_p1 " << expectedCL_p1_asym << std::endl;
    std::cout << "expectedCL_p2 " << expectedCL_p2_asym << std::endl;
    std::cout << "expectedCL_n1 " << expectedCL_n1_asym << std::endl;
    std::cout << "expectedCL_n2 " << expectedCL_n2_asym << std::endl;
  }
  else if (options.Contains("toy")) {
    std::ofstream limitOutput(Form("%s/limits_toy_%s.txt",
				   outputDir.Data(), tag.Data()));
    // Get the actual limit values:
    double observedCL_toy = getIntercept(gCLObs_toy, 0.95);
    double expectedCL_toy = getIntercept(gCLExp_toy, 0.95);
    double expectedCL_p1_toy = getIntercept(gCLExp_toy_p1, 0.95);
    double expectedCL_p2_toy = getIntercept(gCLExp_toy_p2, 0.95);
    double expectedCL_n1_toy = getIntercept(gCLExp_toy_n1, 0.95);
    double expectedCL_n2_toy = getIntercept(gCLExp_toy_n2, 0.95);  
    limitOutput << "observedCL " << observedCL_toy << std::endl;
    limitOutput << "expectedCL " << expectedCL_toy << std::endl;
    limitOutput << "expectedCL_p1 " << expectedCL_p1_toy << std::endl;
    limitOutput << "expectedCL_p2 " << expectedCL_p2_toy << std::endl;
    limitOutput << "expectedCL_n1 " << expectedCL_n1_toy << std::endl;
    limitOutput << "expectedCL_n2 " << expectedCL_n2_toy << std::endl;
    limitOutput.close();
    
    // Then print to screen:
    std::cout << "\nDHCLScan: Results" << std::endl;
    std::cout << "observedCL " << observedCL_toy << std::endl;
    std::cout << "expectedCL " << expectedCL_toy << std::endl;
    std::cout << "expectedCL_p1 " << expectedCL_p1_toy << std::endl;
    std::cout << "expectedCL_p2 " << expectedCL_p2_toy << std::endl;
    std::cout << "expectedCL_n1 " << expectedCL_n1_toy << std::endl;
    std::cout << "expectedCL_n2 " << expectedCL_n2_toy << std::endl;
  }
  
  // Delete pointers, close files, return:
  std::cout << "DHCLScan: Finished!" << std::endl;
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
