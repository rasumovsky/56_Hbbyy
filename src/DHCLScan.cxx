////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHCLScan.cxx                                                        //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 07/08/2015                                                          //
//                                                                            //
//  Includes a main method for using the DHTestStat.cxx class. This is useful //
//  for grid jobs.                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "DHTestStat.h"

double convertToXS(double nEvents) {
  return (nEvents * 0.252155);
}

int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 4) {
    std::cout << "\nUsage: " << argv[0]
	      << " <configFile> <DHSignal> <options>" << std::endl;
    exit(0);
  }
  
  TString configFile = argv[1];
  TString DHSignal = argv[2];
  TString options = argv[3];
  
  // Load the analysis configuration file:
  Config *config = new Config(configFile);
  TString jobName = config->getStr("jobName");
  TString anaType = DHAnalysis::getAnalysisType(config, DHSignal);
  TString outputDir = Form("%s/%s/DHCLScan", 
			   (config->getStr("masterOutput")).Data(),
			   (config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Define the input file, then make a local copy (for remote jobs):
  TString originFile = Form("%s/%s/workspaces/rootfiles/workspaceDH_%s.root",
			    (config->getStr("masterOutput")).Data(), 
			    jobName.Data(), anaType.Data());
  
  CommonFunc::SetAtlasStyle();
  
  TFile inputFile(originFile, "read");
  RooWorkspace *workspace = (RooWorkspace*)inputFile.Get("combinedWS");
  
  // Instantiate the test statistic class for calculations and plots:
  DHTestStat *ts = new DHTestStat(configFile, DHSignal, "new", workspace);
  
  TGraph *gCLExp = new TGraph();
  TGraph *gCLObs = new TGraph();
  
  // Loop over number of accepted events:
  for (int i_e = 0; i_e < 40; i_e++) {
    // Set output directory for plots:
    ts->setPlotDirectory(outputDir.Data());
    ts->setParams("nDH_NonResSR", ((double)i_e+1.0)/2.0, true);
    
    // Calculate the 95% CL value:
    ts->calculateNewCL();
    
    // Check that the fit converges:
    if (ts->fitsAllConverged()) {
      double currExpCL = ts->accessValue("CL", false, 0);
      double currObsCL = ts->accessValue("CL", true, 0);
      // fill graphs
      gCLExp->SetPoint(i_e, convertToXS(i_e), currExpCL);
      gCLObs->SetPoint(i_e, convertToXS(i_e), currObsCL);
    }
    
    ts->clearData();
  }
  
  TCanvas *can = new TCanvas("can","can");
  can->cd();
  gCLExp->GetXaxis()->SetTitle("#sigma_{BSM} [pb]");
  gCLObs->GetXaxis()->SetTitle("#sigma_{BSM} [pb]");
  gCLExp->GetYaxis()->SetTitle("95% CL");
  gCLObs->GetYaxis()->SetTitle("95% CL");
  gCLExp->SetLineColor(kBlack);
  gCLExp->SetMarkerColor(kBlack);
  gCLExp->SetMarkerColor(kBlack);
  gCLExp->SetMarkerStyle(2);
  gCLObs->SetMarkerColor(kBlue);
  gCLObs->SetMarkerStyle(2);
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
  line->DrawLine(gCLObs->GetXaxis()->GetXmin(), 0.95,
		 gCLObs->GetXaxis()->GetXmax(), 0.95);
  can->Print(Form("%s/scan95CL.eps",outputDir.Data()));
  
  std::cout << "DHCLScan: Finished!" << std::endl;
 
  // Instantiate the test statistic class for calculations and plots:
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
