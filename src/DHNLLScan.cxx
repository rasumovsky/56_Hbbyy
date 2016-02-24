////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHNLLScan.cxx                                                        /
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 24/02/2016                                                          //
//                                                                            //
//  Performs an NLL scan.                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "DHTestStat.h"

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
  TString anaType = config->getStr("AnalysisType");
  TString outputDir = Form("%s/%s/DHNLLScan", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Define the input file, then make a local copy (for remote jobs):
  TString originFile = Form("%s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root",
			    (config->getStr("MasterOutput")).Data(), 
			    (config->getStr("JobName")).Data(),
			    anaType.Data());
  
  // Create a nametag for results:
  TString tag = (anaType.EqualTo("Resonant")) ? 
    Form("ResonantMX%d", resonanceMass) : "NonResonant";
  if (options.Contains("NEvents")) tag.Append("_NEvent");
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Instantiate the statistics class:
  DHTestStat *dhts = new DHTestStat(configFile, 
				    config->getStr("TestStatOptions"),
				    NULL);
  dhts->setPlotDirectory(outputDir);
  
  // Non-resonant analysis: specify cross-section:
  dhts->setParam(config->getStr("CrossSectionVar"),
		 config->getNum("CrossSectionValue"), true);
  
  // Resonant analysis: specify cross-section and resonance mass:
  if (anaType.EqualTo("Resonant")) {
    dhts->setParam(config->getStr("ResonanceMassVar"), resonanceMass, true);
  }
  
  // Scan 1: Statistical scan, freezes all systematics! (stat only):
  std::map<int,double> scanStat
    = dhts->scanNLL("statistical", config->getStr("DataToScan"), 
		    config->getStr("ParamToScan"),
		    config->getStrV("SysSources"));
  TGraph *gStat = dhts->nllScanGraph();

  
  // Scan 2: Experimental scan, freezes theory systematics! (stat + sys):
  std::map<int,double> scanExp 
    = dhts->scanNLL("experimental", config->getStr("DataToScan"), 
		    config->getStr("ParamToScan"),
		    config->getStrV("UncertaintiesTheory"));
  TGraph *gSys = dhts->nllScanGraph();
  
  // Scan 3: Nominal scan, don't freeze anything! (stat + sys + th)
  std::vector<TString> null; null.clear();
  std::map<int,double> scanNom
    = dhts->scanNLL("nominal", config->getStr("DataToScan"),
		    config->getStr("ParamToScan"), null);
  TGraph *gAll = dhts->nllScanGraph();
  
  // Check that fits succeeded:
  if (!dhts->fitsAllConverged()) {
    std::cout << "\nDMMaster: ERROR! Statistics fits failed.\n" << std::endl;
  }
  
  delete dhts;
  
  // Check that TGraphs have been recovered:
  if (!gAll || !gStat || !gSys) {
    std::cout << "DHNLLScan: ERROR! Histograms are null." << std::endl;
    exit(0);
  }

  // Draw the graph:
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->cd();
  gStat->SetLineColor(kRed+1);
  gSys->SetLineColor(kBlue+1);
  gAll->SetLineColor(1);
  gAll->SetLineWidth(2);
  gStat->SetLineWidth(2);
  gSys->SetLineWidth(2);
  gAll->SetLineStyle(1);
  gStat->SetLineStyle(2);
  gSys->SetLineStyle(5);
  gAll->GetXaxis()->SetTitle(config->getStr("ParamToScan"));
  gAll->GetYaxis()->SetTitle("-2#DeltaNLL");
  gAll->GetYaxis()->SetRangeUser(0, gAll->GetMaximum());
  gAll->Draw("AL");
  gStat->Draw("LSAME");
  gSys->Draw("LSAME");
    
  // 1 sigma and 2 sigma lines:
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->DrawLine(gAll->GetXaxis()->GetXmin(), 1.0,
		 gAll->GetXaxis()->GetXmax(), 1.0);
  line->DrawLine(gAll->GetXaxis()->GetXmin(), 4.0,
		 gAll->GetXaxis()->GetXmax(), 4.0);
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.20, 0.86, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.32, 0.86, config->getStr("ATLASLabel"));
  t.DrawLatex(0.20, 0.80, Form("#sqrt{s} = 13 TeV: %2.1f fb^{-1}",
			       (config->getNum("AnalysisLuminosity")/
				1000.0)));
  
  TLegend leg(0.55, 0.76, 0.85, 0.9);
  leg.SetFillColor(0);
  leg.SetTextSize(0.05);
  leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.AddEntry(gStat, "Stat.", "l");
  leg.AddEntry(gSys, "Stat. + Syst.", "l");
  leg.AddEntry(gAll, "Stat. + Syst. + Th.", "l");
  leg.Draw("SAME");
  
  // Print the canvas:
  canvas->Print(Form("%s/scanNLL_%s.eps", outputDir.Data(), tag.Data()));
  canvas->Print(Form("%s/scanNLL_%s.C", outputDir.Data(), tag.Data()));
    
  // Print the NLL scan results!
  std::cout << "\n----------------------------------------\n" 
	    << "Printing NLL scan results:" << std::endl;
  std::cout << "Statistics-only: " 
	    << scanStat[0] << " +" << scanStat[1] << " -" << scanStat[-1]
	    << std::endl;
  std::cout << "Statistics + systematics:" 
	    << scanExp[0] << " +" << scanExp[1] << " -" << scanExp[-1]
	    << std::endl;
  std::cout << "Statistics + Systematics + Theory: " 
	    << scanNom[0] << " +" << scanNom[1] << " -" << scanNom[-1] 
	    << "\n" << std::endl;
  
  return 0;
}
