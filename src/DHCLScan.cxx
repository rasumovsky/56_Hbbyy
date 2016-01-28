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
  
  TFile wsFile(originFile, "read");
  RooWorkspace *workspace = (RooWorkspace*)wsFile.Get("combinedWS");
  
  // Instantiate the test statistic class for calculations and plots:
  DHTestStat *ts = new DHTestStat(configFile, "new", workspace);
    
  // Arrays to store band information:
  double varValues_asym[100];
  double CLObs_asym[100];
  double CLExp_asym_p2[100];
  double CLExp_asym_p1[100];
  double CLExp_asym[100];
  double CLExp_asym_n1[100];
  double CLExp_asym_n2[100];
  double varValues_toy[100];  
  double CLObs_toy[100];
  double CLExp_toy_p2[100];  
  double CLExp_toy_p1[100];
  double CLExp_toy[100];
  double CLExp_toy_n1[100];
  double CLExp_toy_n2[100];
  
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
      
      // Set output directory for plots:
      ts->setPlotDirectory(outputDir.Data());
      // Set the value of the variable to scan:
      ts->setParam(config->getStr("CLScanVar"), crossSection, true);
      
      // Asymptotic calculation of CL:
      if (options.Contains("asymptotic") || options.Contains("both")) {
	
	// Calculate the 95% CL value:
	ts->calculateNewCL();
	
	// Check that the fit converges:
	if (ts->fitsAllConverged()) {
	  varValues_asym[n_asym] = crossSection;
	  CLObs_asym[n_asym] = ts->accessValue("CL", true, 0);
	  CLExp_asym[n_asym] = ts->accessValue("CL", false, 0);
	  CLExp_asym_p2[n_asym] = ts->accessValue("CL", false, 2);
	  CLExp_asym_p1[n_asym] = ts->accessValue("CL", false, 1);
	  CLExp_asym_n1[n_asym] = ts->accessValue("CL", false, -1);
	  CLExp_asym_n2[n_asym] = ts->accessValue("CL", false, -2);
	  
	  outFile_asym << varValues_asym[n_asym] << " " 
		       << CLObs_asym[n_asym] << " " 
		       << CLExp_asym[n_asym] << " " 
		       << CLExp_asym_p2[n_asym] << " " 
		       << CLExp_asym_p1[n_asym] << " " 
		       << CLExp_asym_n1[n_asym] << " " 
		       << CLExp_asym_n2[n_asym] << std::endl;
	}
	n_asym++;
	ts->clearData();
      }
      
      // Pseudo-experiment calculation of CL:
      if (options.Contains("toy") || options.Contains("both")) {
	
	// Calculate the observed qMu values:
	double obsPoI = 0.0;
	double nllObsMu1 = ts->getFitNLL("obsData", 1, true, obsPoI, false);
	double nllObsMuFree = ts->getFitNLL("obsData", 1, false, obsPoI, false);
	double qMuObs = ts->getQMuFromNLL(nllObsMu1, nllObsMuFree, obsPoI, 1);
	
	// Calculate the expected qMu values (see DHTestStat m_dataForExpQMu):
	double expPoI = 0.0;
	double nllExpMu1=ts->getFitNLL("asimovDataMu0",1,true,expPoI,false);
	double nllExpMuFree=ts->getFitNLL("asimovDataMu0",1,false,expPoI,false);
	double qMuExp = ts->getQMuFromNLL(nllExpMu1, nllExpMuFree, expPoI, 1);
	
	// Then specify file in DHToyAnalysis to load:
	TString toyScanOption = Form("CLScan%d",n_toy);
	DHToyAnalysis *dhta = new DHToyAnalysis(configFile, toyScanOption);
	
	varValues_toy[n_toy] = crossSection;
	CLObs_toy[n_toy] = dhta->calculateCLFromToy(qMuObs, 0);
	CLExp_toy[n_toy] = dhta->calculateCLFromToy(qMuExp, 0);
	CLExp_toy_p2[n_toy] = dhta->calculateCLFromToy(qMuExp, 2);
	CLExp_toy_p1[n_toy] = dhta->calculateCLFromToy(qMuExp, 1);
	CLExp_toy_n1[n_toy] = dhta->calculateCLFromToy(qMuExp, -1);
	CLExp_toy_n2[n_toy] = dhta->calculateCLFromToy(qMuExp, -2);
	
	outFile_toy << varValues_toy[n_toy] << " " 
		    << CLObs_toy[n_toy] << " " 
		    << CLExp_toy[n_toy] << " "
		    << CLExp_toy_p2[n_toy] << " " 
		    << CLExp_toy_p1[n_toy] << " " 
		    << CLExp_toy_n1[n_toy] << " " 
		    << CLExp_toy_n2[n_toy] << std::endl;
	n_toy++;
	delete dhta;
      }
    }
    
    // Close the files that save CL data:
    outFile_asym.close();
    outFile_toy.close();
  }
  
  //----------------------------------------//
  // Plot the results:
  
  // Median expected and observed results:
  TGraph *gCLExp_asym = new TGraph(n_asym, varValues_asym, CLExp_asym);
  TGraph *gCLObs_asym = new TGraph(n_asym, varValues_asym, CLObs_asym);
  TGraph *gCLExp_toy  = new TGraph(n_toy, varValues_asym, CLExp_toy);
  TGraph *gCLObs_toy  = new TGraph(n_toy, varValues_asym, CLObs_toy);
  
  // Also plot the bands:
  TGraphAsymmErrors *gCLExp_asym_2s 
    = new TGraphAsymmErrors(n_asym, varValues_asym, CLExp_asym, 0, 0, 
			    CLExp_asym_n2, CLExp_asym_p2);
  TGraphAsymmErrors *gCLExp_asym_1s
    = new TGraphAsymmErrors(n_asym, varValues_asym, CLExp_asym, 0, 0, 
			    CLExp_asym_n1, CLExp_asym_p1);
  TGraphAsymmErrors *gCLExp_toy_2s 
    = new TGraphAsymmErrors(n_toy, varValues_toy, CLExp_toy, 0, 0, 
			    CLExp_toy_n2, CLExp_toy_p2);
  TGraphAsymmErrors *gCLExp_toy_1s
    = new TGraphAsymmErrors(n_toy, varValues_toy, CLExp_toy, 0, 0, 
			    CLExp_toy_n1, CLExp_toy_p1);
  
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
  gCLObs_toy->GetXaxis()->SetRangeUser(xsMin, xsMax);
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
  gCLObs_asym->GetXaxis()->SetRangeUser(xsMin, xsMax);
  gCLExp_asym_2s->SetFillColor(kYellow);
  gCLExp_asym_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.56,0.18,0.88,0.34);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.05);
  leg.AddEntry(gCLObs_toy,"Obs. Limit","l");
  leg.AddEntry(gCLExp_toy,"Exp. Limit","l");
  leg.AddEntry(gCLExp_toy_1s,"Exp. Limit #pm1#sigma_{exp}","F");
  leg.AddEntry(gCLExp_toy_2s,"Exp. Limit #pm2#sigma_{exp}","F");
  
  // Plotting options:
  if (options.Contains("toy")) {
    gCLExp_toy->Draw("AL");
    gCLExp_toy_2s->Draw("3same");
    gCLExp_toy_1s->Draw("3same");
    gCLExp_toy->Draw("LSAME");
    gCLObs_toy->Draw("LSAME");
  }
  else if (options.Contains("asymptotic")) {
    gCLExp_asym->Draw("AL");
    gCLExp_asym_2s->Draw("3same");
    gCLExp_asym_1s->Draw("3same");
    gCLExp_asym->Draw("LSAME");
    gCLObs_asym->Draw("LSAME");
  }
  else if (options.Contains("both")) {
    gCLExp_toy->Draw("AL");
    gCLExp_toy_2s->Draw("3same");
    gCLExp_toy_1s->Draw("3same");
    gCLExp_toy->Draw("LSAME");
    gCLObs_toy->Draw("LSAME");
    
    gCLExp_asym->SetLineColor(kBlue+2);
    gCLObs_asym->SetLineColor(kBlue+2);
    gCLExp_asym->SetLineStyle(2);
    gCLObs_asym->SetLineStyle(2);
    
    gCLExp_asym->Draw("LSAME");
    gCLObs_asym->Draw("LSAME");
    leg.AddEntry(gCLObs_asym,"Obs. Limit (asymptotic)","l");
    leg.AddEntry(gCLExp_asym,"Exp. Limit (asymptotic)","l");
  }
  leg.Draw("SAME");
  
  // 95% CL Line
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kRed);
  line->DrawLine(xsMin, 0.95, xsMax, 0.95);
  can->Print(Form("%s/scan95CL.eps",outputDir.Data()));
  
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
  delete ts;
  delete workspace;
  delete config;
  wsFile.Close();
  return 0;
}
