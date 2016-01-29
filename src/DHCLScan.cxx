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
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "DHTestStat.h"
#include "DHToyAnalysis.h"

/**
   -----------------------------------------------------------------------------
   Get the intersection point for the graph.

double getIntercept(TGraph *graph, double valueToIntercept) {
  
  double xLow = graph->GetMinimum();
  double xHigh = graph->GetMaximum();
  
  // Check that it is possible to find the intercept:
  if (valueToIntercept > xLow || valueToIntercept < xHigh) {
    std::cout << "DHCLScan: ERROR! Intercept not in graph range." << std::cout;
    exit(0);
  }
  
  // Loop over points in the graph
  for (int i_p = 0; i_p < graph->GetN(); i_p++) {
    double xCurr; double yCurr;
    graph->GetPoint(i_p, xCurr, yCurr);
    
    if (yCurr >= valueToIntercept && yCurr < graph->Eval(xHigh)) {
      xHigh = xCurr;
    }
    else if (yCurr <= valueToIntercept && yCurr > graph->Eval(xLow)) {
      xLow = xCurr;
    }
  }
  
  // Then use bisection method:
  
  
}
*/

/**
   -----------------------------------------------------------------------------
   The main method scans the 95% CL for various signal cross-sections.
   @param configFile - The analysis configuration file.
   @param options - Job options.
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 3) {
    std::cout << "\nUsage: " << argv[0]
	      << " <configFile> <options>" << std::endl;
    exit(0);
  }
  
  TString configFile = argv[1];
  TString options = argv[2];
  
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
  
  int n_asym = 0;
  int n_toy = 0;
  
  // Scan information:
  double xsMin = config->getNum("CLScanMin");
  double xsMax = config->getNum("CLScanMax");
  double step = config->getNum("CLScanStep");
  
  //----------------------------------------//
  // Open CL values from file:
  if (options.Contains("FromFile")) {
    
    // Open the saved CL values from asymptotics:
    if (options.Contains("asymptotic") || options.Contains("both")) {
      ifstream inputFile_asym(Form("%s/scan_CLvalues_asym.txt",
				   outputDir.Data()));
      if (inputFile_asym.is_open()) {
	while (!inputFile_asym.eof()) {
	  inputFile_asym >> varValues_asym[n_asym] >> CLObs_asym[n_asym] 
			 >> CLExp_asym[n_asym] >> CLExp_asym_p2[n_asym] 
			 >> CLExp_asym_p1[n_asym] >> CLExp_asym_n1[n_asym] 
			 >> CLExp_asym_n2[n_asym];
	  n_asym++;
	}
      }
      inputFile_asym.close();
    }
    
    // Open the saved CL values from toys:
    if (options.Contains("toy") || options.Contains("both")) {
      ifstream inputFile_toy(Form("%s/scan_CLvalues_toy.txt",
				  outputDir.Data()));
      if (inputFile_toy.is_open()) {
	while (!inputFile_toy.eof()) {
	  inputFile_toy >> varValues_toy[n_toy] >> CLObs_toy[n_toy] 
			>> CLExp_toy[n_toy] >> CLExp_toy_p2[n_toy] 
			>> CLExp_toy_p1[n_toy] >> CLExp_toy_n1[n_toy] 
			>> CLExp_toy_n2[n_toy];
	  n_toy++;
	}
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
      outFile_asym.open(Form("%s/scan_CLvalues_asym.txt", outputDir.Data()));
    }
    if (options.Contains("toy") || options.Contains("both")) {
      outFile_toy.open(Form("%s/scan_CLvalues_toy.txt", outputDir.Data()));
    }
    
    // Loop over number of accepted events:
    for (double crossSection = xsMin; crossSection < xsMax;
	 crossSection += step) {
      std::cout << "DHCLScan: cross-section = " << crossSection << std::endl;
      
      // Set output directory for plots:
      //dhts->setPlotDirectory(outputDir.Data());
      // Set the value of the variable to scan:
      dhts->setParam(config->getStr("CLScanVar"), crossSection, true);
      
      // Asymptotic calculation of CL:
      if (options.Contains("asymptotic") || options.Contains("both")) {
	std::cout << "DHCLScan: calculating asymptotic CL" << std::endl;
	
	// Calculate the 95% CL value:
	dhts->calculateNewCL();
	
	// Check that the fit converges:
	if (dhts->fitsAllConverged()) {
	  varValues_asym[n_asym] = crossSection;
	  CLObs_asym[n_asym] = dhts->accessValue("CL", true, 0);
	  CLExp_asym[n_asym] = dhts->accessValue("CL", false, 0);
	  CLExp_asym_p2[n_asym] = dhts->accessValue("CL", false, 2);
	  CLExp_asym_p1[n_asym] = dhts->accessValue("CL", false, 1);
	  CLExp_asym_n1[n_asym] = dhts->accessValue("CL", false, -1);
	  CLExp_asym_n2[n_asym] = dhts->accessValue("CL", false, -2);
	  
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
	qMuObs_toy[n_toy] = dhts->getQMuFromNLL(nllObsMu1, nllObsMuFree,
						obsPoI, 1);
	
	// Calculate the expected qMu values (see DHTestStat m_dataForExpQMu):
	double expPoI = 0.0;
	double nllExpMu1
	  = dhts->getFitNLL("asimovDataMu0", 1, true, expPoI, false);
	double nllExpMuFree
	  = dhts->getFitNLL("asimovDataMu0", 1, false, expPoI, false);
	qMuExp_toy[n_toy] = dhts->getQMuFromNLL(nllExpMu1, nllExpMuFree, 
						expPoI, 1);
	
	std::cout << "nllObsMu1=" << nllObsMu1 
		  << ", nllObsMuFree=" << nllObsMuFree
		  << ", qMuObs=" << qMuObs_toy[n_toy] << std::endl;
	
	std::cout << "nllExpMu1=" << nllExpMu1
		  << ", nllExpMuFree=" << nllExpMuFree
		  << ", qMuExp=" << qMuExp_toy[n_toy] << std::endl;
	
	
	varValues_toy[n_toy] = crossSection;
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
	DHToyAnalysis *dhta = new DHToyAnalysis(configFile, toyScanOption);
	if (!(dhta->areInputFilesOK())) {
	  std::cout << "DHCLScan: ERROR with toy scan option " << toyScanOption
		    << std::endl;
	  continue;
	}
	
	std::cout << "\nObs: " << std::endl;
	CLObs_toy[i_t] = dhta->calculateCLFromToy(qMuObs_toy[i_t], 0);
	std::cout << "\nExp: " << std::endl;
	CLExp_toy[i_t] = dhta->calculateCLFromToy(qMuExp_toy[i_t], 0);
	std::cout << "\nExp_p2: " << std::endl;
	CLExp_toy_p2[i_t] = dhta->calculateCLFromToy(qMuExp_toy[i_t], 2);
	std::cout << "\nExp_p1: " << std::endl;
	CLExp_toy_p1[i_t] = dhta->calculateCLFromToy(qMuExp_toy[i_t], 1);
	std::cout << "\nExp_n1: " << std::endl;
	CLExp_toy_n1[i_t] = dhta->calculateCLFromToy(qMuExp_toy[i_t], -1);
	std::cout << "\nExp_n2: " << std::endl;
	CLExp_toy_n2[i_t] = dhta->calculateCLFromToy(qMuExp_toy[i_t], -2);
	
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
  
  n_asym--;
  n_toy--;
  
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
  TGraph *gCLExp_asym = new TGraph(n_asym, varValues_asym, CLExp_asym);
  TGraph *gCLObs_asym = new TGraph(n_asym, varValues_asym, CLObs_asym);
  TGraph *gCLExp_toy = new TGraph(n_toy, varValues_toy, CLExp_toy);
  TGraph *gCLObs_toy = new TGraph(n_toy, varValues_toy, CLObs_toy);
  
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
  gCLExp_toy->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
  gCLObs_toy->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
  gCLExp_toy->GetYaxis()->SetTitle("CL");
  gCLObs_toy->GetYaxis()->SetTitle("CL");
  gCLExp_toy->SetLineColor(kBlack);
  gCLObs_toy->SetLineColor(kBlack);
  gCLExp_toy->SetLineStyle(2);
  gCLObs_toy->SetLineStyle(1);
  gCLExp_toy->SetLineWidth(2);
  gCLObs_toy->SetLineWidth(2);
  //gCLObs_toy->GetXaxis()->SetRangeUser(xsMin, xsMax);
  gCLExp_toy->GetYaxis()->SetRangeUser(0.0, 1.0);
  gCLExp_toy_2s->SetFillColor(kYellow);
  gCLExp_toy_1s->SetFillColor(kGreen);

  // Asymptotic graph formatting:
  gCLExp_asym->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
  gCLObs_asym->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
  gCLExp_asym->GetYaxis()->SetTitle("CL");
  gCLObs_asym->GetYaxis()->SetTitle("CL");
  gCLExp_asym->SetLineColor(kBlack);
  gCLObs_asym->SetLineColor(kBlack);
  gCLExp_asym->SetLineStyle(2);
  gCLObs_asym->SetLineStyle(1);
  gCLExp_asym->SetLineWidth(2);
  gCLObs_asym->SetLineWidth(2);
  //gCLObs_asym->GetXaxis()->SetRangeUser(xsMin, xsMax);
  gCLExp_asym->GetYaxis()->SetRangeUser(0.0, 1.0);
  gCLExp_asym_2s->SetFillColor(kYellow);
  gCLExp_asym_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.53,0.19,0.88,0.37);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!config->getBool("DoBlind")) leg.AddEntry(gCLObs_toy,"Obs. Limit","l");
  leg.AddEntry(gCLExp_toy,"Exp. Limit","l");
  leg.AddEntry(gCLExp_toy_1s,"Exp. Limit #pm1#sigma_{exp}","F");
  leg.AddEntry(gCLExp_toy_2s,"Exp. Limit #pm2#sigma_{exp}","F");
  
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
      leg.AddEntry(gCLObs_asym,"Obs. Limit (asymptotic)","l");
    }
    leg.AddEntry(gCLExp_asym,"Exp. Limit (asymptotic)","l");
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
  
  if (options.Contains("asymptotic")) {
    can->Print(Form("%s/scan95CL_asymptotic.eps",outputDir.Data()));
  }
  else if (options.Contains("toy")) {
    can->Print(Form("%s/scan95CL_toy.eps",outputDir.Data()));
  }
  else {
    can->Print(Form("%s/scan95CL_comparison.eps",outputDir.Data()));
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
