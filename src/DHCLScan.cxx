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
  TGraph *gCLExp = new TGraph();
  TGraph *gCLObs = new TGraph();
  
  // Loop over number of accepted events:
  int i_e = 0;
  double xsMin = config->getNum("CLScanMin");
  double xsMax = config->getNum("CLScanMax");
  double step = config->getNum("CLScanStep");
  for (double crossSection = xsMin; crossSection < xsMax; crossSection += step){
    
    // Set output directory for plots:
    ts->setPlotDirectory(outputDir.Data());
    ts->setParam(config->getStr("CLScanVar"), crossSection, true);
    
    // Calculate the 95% CL value:
    ts->calculateNewCL();
    
    // Check that the fit converges:
    if (ts->fitsAllConverged()) {
      double currExpCL = ts->accessValue("CL", false, 0);
      double currObsCL = ts->accessValue("CL", true, 0);
      // fill graphs (converting fb to pb):
      gCLExp->SetPoint(i_e, crossSection, currExpCL);
      gCLObs->SetPoint(i_e, crossSection, currObsCL);
      i_e++;
    }
    ts->clearData();
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
