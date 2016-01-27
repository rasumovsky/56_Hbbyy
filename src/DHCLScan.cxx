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
  TString options = argv[3];
  
  // Load the analysis configuration file:
  Config *config = new Config(configFile);
  TString jobName = config->getStr("JobName");
  TString anaType = config->getStr("AnalysisType");
  TString outputDir = Form("%s/%s/DHCLScan", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Define the input file, then make a local copy (for remote jobs):
  TString originFile = Form("%s/%s/workspaces/rootfiles/workspaceDH_%s.root",
			    (config->getStr("MasterOutput")).Data(), 
			    jobName.Data(), anaType.Data());
  
  CommonFunc::SetAtlasStyle();
  
  TFile inputFile(originFile, "read");
  RooWorkspace *workspace = (RooWorkspace*)inputFile.Get("combinedWS");
  
  // Instantiate the test statistic class for calculations and plots:
  DHTestStat *ts = new DHTestStat(configFile, "new", workspace);
  
  // Expected and observed results:
  TGraph *gCLExp_asym = new TGraph();
  TGraph *gCLObs_asym = new TGraph();
  TGraph *gCLExp_toy = new TGraph();
  TGraph *gCLObs_toy = new TGraph();
  
  // Loop over number of accepted events:
  int i_e = 0;
  double xsMin = config->getNum("CLScanMin");
  double xsMax = config->getNum("CLScanMax");
  double step = config->getNum("CLScanStep");
  for (double crossSection = xsMin; crossSection < xsMax; crossSection += step){

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
	double asymCLObs = ts->accessValue("CL", true, 0);
	double asymCLExp = ts->accessValue("CL", false, 0);
	double asymCLExpN2 = ts->accessValue("CL", false, -2);
	double asymCLExpN1 = ts->accessValue("CL", false, -1);
	double asymCLExpP1 = ts->accessValue("CL", false, 1);
	double asymCLExpP2 = ts->accessValue("CL", false, 2);
	
	// Fill graphs:
	gCLExp_asym->SetPoint(i_e, crossSection, asymCLExp);
	gCLObs_asym->SetPoint(i_e, crossSection, asymCLObs);
      }
      ts->clearData();
    }
    
    // Pseudo-experiment calculation of CL:
    else if (options.Contains("toy") || options.Contains("both")) {
      
      // Calculate the observed qMu values:
      double obsPoI = 0.0;
      double nllObsMu1 = ts->getFitNLL("obsData", 1, true, obsPoI, false);
      double nllObsMuFree = ts->getFitNLL("obsData", 1, false, obsPoI, false);
      double qMuObs = m_dhts->getQMuFromNLL(nllObsMu1, nllObsMuFree, obsPoI, 1);
      
      // Calculate the expected qMu values (see DHTestStat m_dataForExpQMu):
      double expPoI = 0.0;
      double nllExpMu1 = ts->getFitNLL("asimovDataMu0", 1, true, expPoI, false);
      double nllExpMuFree = ts->getFitNLL("asimovDataMu0",1,false,expPoI,false);
      double qMuExp = m_dhts->getQMuFromNLL(nllExpMu1, nllExpMuFree, expPoI, 1);
      
      /////////////////

      // NEED TO SPECIFY FILE SOMEHOW
      // use i_e I think...
      /////////////////

      // Then specify file in DHToyAnalysis to load:
      DHToyAnalysis *dhta = new DHToyAnalysis(configFileName);
      
      
      double toyCLObs = dhta->calculateCLFromToy(qMuObs, 0);
      double toyCLExp = dhta->calculateCLFromToy(qMuExp, 0);
      double toyCLExpN2 = dhta->calculateCLFromToy(qMuExp, -2);
      double toyCLExpN1 = dhta->calculateCLFromToy(qMuExp, -1);
      double toyCLExpP1 = dhta->calculateCLFromToy(qMuExp, 1);
      double toyCLExpP2 = dhta->calculateCLFromToy(qMuExp, 2);
      
      // Fill graphs:
      gCLExp_toy->SetPoint(i_e, crossSection, toyCLExp);
      gCLObs_toy->SetPoint(i_e, crossSection, toyCLObs);
    }
    
    i_e++;
  }

  TCanvas *can = new TCanvas("can","can");
  can->cd();
  gCLExp->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
  gCLObs->GetXaxis()->SetTitle(config->getStr("CLScanPrintName"));
  gCLExp->GetYaxis()->SetTitle("95% CL");
  gCLObs->GetYaxis()->SetTitle("95% CL");
  gCLExp->SetLineColor(kBlack);
  gCLExp->SetLineWidth(2);
  gCLExp->SetMarkerColor(kBlack);
  gCLExp->SetMarkerColor(kBlack);
  gCLExp->SetMarkerStyle(2);
  gCLObs->SetLineColor(kBlue);
  gCLObs->SetLineWidth(2);
  gCLObs->SetMarkerColor(kBlue);
  gCLObs->SetMarkerStyle(2);
  gCLObs->GetXaxis()->SetRangeUser(xsMin, xsMax);
  gCLObs->Draw("ALP");
  gCLExp->Draw("LPSAME");
  TLegend leg(0.7,0.2,0.88,0.3);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.05);
  leg.AddEntry(gCLExp, "Expected", "P");
  leg.AddEntry(gCLObs, "Observed", "P");
  leg.Draw("SAME");
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
  delete gCLObs;
  delete gCLExp;
  delete ts;
  delete workspace;
  delete config;
  inputFile.Close();
  return 0;
}
