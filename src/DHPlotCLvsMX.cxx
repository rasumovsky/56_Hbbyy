////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DHPlotCLvsMX.cxx                                                    //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/02/2016                                                          //
//                                                                            //
//  Plots the resonant analysis limits as a function of MX.                   //
//                                                                            //
//  Macro options:                                                            //
//  - "toy"        Get limits from toy MC jobs.                               //
//  - "asymptotic" Get asymptotic limits.                                     //
//  - "both"       Compares toy and asymptotic limits.                        //
//  - "NEvents"    Calculate limits in terms of number of events.             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"

/**
   -----------------------------------------------------------------------------
   The main method scans the 95% CL for various signal cross-sections.
   @param configFile - The analysis configuration file.
   @param options - Job options. Can be "toy" or "asymptotic" or "both"
   @param resMass - The resonance mass.
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
  TString outputDir = Form("%s/%s/DHPlotCLvsMX", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  TString inputDir = Form("%s/%s/DHCLScan", 
			  (config->getStr("MasterOutput")).Data(),
			  (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Check that we are doing resonant analysis:
  if (!(config->getStr("AnalysisType")).EqualTo("Resonant")) {
    std::cout << "DHPlotCLvsMX: Cannot plot CL vs MX for nonresonant analysis!"
	      << std::endl;
  }
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
        
  // Arrays to store band information:
  double varValues_asym[100] = {0};
  double CLObs_asym[100]     = {0};
  double CLExp_asym_p2[100]  = {0};
  double CLExp_asym_p1[100]  = {0};
  double CLExp_asym[100]     = {0};
  double CLExp_asym_n1[100]  = {0};
  double CLExp_asym_n2[100]  = {0};

  double varValues_toy[100] = {0};  
  double CLObs_toy[100]     = {0};
  double CLExp_toy_p2[100]  = {0};  
  double CLExp_toy_p1[100]  = {0};
  double CLExp_toy[100]     = {0};
  double CLExp_toy_n1[100]  = {0};
  double CLExp_toy_n2[100]  = {0};
    
  int n_asym = 0;
  int n_toy = 0;
  
  // Scan information:
  std::vector<double> scanMXValues = config->getNumV("MXScanValues");
    
  //----------------------------------------//
  // Open CL values from files (different file for each MX):
  for (int i_m = 0; i_m < (int)scanMXValues.size(); i_m++) {
    int mX = (int)(scanMXValues[i_m]);
    
    // Open the saved CL values from asymptotics:
    if (options.Contains("asymptotic") || options.Contains("both")) {
      TString fileAsymName = (options.Contains("NEvents")) ? 
	Form("%s/limits_asymptotic_ResonantMX%d_NEvent.txt",inputDir.Data(),mX):
	Form("%s/limits_asymptotic_ResonantMX%d.txt", inputDir.Data(), mX);
      ifstream inputFile_asym(fileAsymName);
      if (inputFile_asym.is_open()) {
	TString inName = ""; double inValue = 0.0;
	while (inputFile_asym >> inName >> inValue) {
	  if (!options.Contains("NEvents")) inValue = inValue / 1000.0;
	  // WARNING! The + and - are deliberately swapped here!
	  if (inName.EqualTo("observedCL")) CLObs_asym[n_asym] = inValue;
	  if (inName.EqualTo("expectedCL")) CLExp_asym[n_asym] = inValue;
	  if (inName.EqualTo("expectedCL_n2")) CLExp_asym_p2[n_asym] = inValue;
	  if (inName.EqualTo("expectedCL_n1")) CLExp_asym_p1[n_asym] = inValue;
	  if (inName.EqualTo("expectedCL_p1")) CLExp_asym_n1[n_asym] = inValue;
	  if (inName.EqualTo("expectedCL_p2")) CLExp_asym_n2[n_asym] = inValue;
	}
      }
      inputFile_asym.close();
      varValues_asym[n_asym] = mX;
      n_asym++;
    }
    
    // Open the saved CL values from toys:
    if (options.Contains("toy") || options.Contains("both")) {
      TString fileToyName = (options.Contains("NEvents")) ? 
	Form("%s/limits_toy_ResonantMX%d_NEvent.txt", inputDir.Data(), mX) :
	Form("%s/limits_toy_ResonantMX%d.txt", inputDir.Data(), mX);
      ifstream inputFile_toy(fileToyName);
      if (inputFile_toy.is_open()) {
	TString inName = ""; double inValue = 0.0;
	while (inputFile_toy >> inName >> inValue) {
	  if (!options.Contains("NEvents")) inValue = inValue / 1000.0;
	  // WARNING! The + and - are deliberately swapped here!
	  if (inName.EqualTo("observedCL")) CLObs_toy[n_toy] = inValue;
	  if (inName.EqualTo("expectedCL")) CLExp_toy[n_toy] = inValue;
	  if (inName.EqualTo("expectedCL_n2")) CLExp_toy_p2[n_toy] = inValue;
	  if (inName.EqualTo("expectedCL_n1")) CLExp_toy_p1[n_toy] = inValue;
	  if (inName.EqualTo("expectedCL_p1")) CLExp_toy_n1[n_toy] = inValue;
	  if (inName.EqualTo("expectedCL_p2")) CLExp_toy_n2[n_toy] = inValue;
	}
      }
      inputFile_toy.close();
      varValues_toy[n_toy] = mX;
      n_toy++;
    }
  }
  
  //----------------------------------------//
  // Plot the results:
  
  // Calculate errors for plot:
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
  TGraph *gCLExp_toy  = new TGraph(n_toy,  varValues_toy,  CLExp_toy);
  TGraph *gCLObs_toy  = new TGraph(n_toy,  varValues_toy,  CLObs_toy);
  
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
  gCLExp_toy->GetXaxis()->SetTitle("m_{X} [GeV]");
  gCLObs_toy->GetXaxis()->SetTitle("m_{X} [GeV]");
  gCLExp_toy_2s->GetXaxis()->SetTitle("m_{X} [GeV]");
  if (options.Contains("NEvents")) {
    gCLExp_toy->GetYaxis()
      ->SetTitle("95% CL limit on X#rightarrowhh event yield");
    gCLObs_toy->GetYaxis()
      ->SetTitle("95% CL limit on X#rightarrowhh event yield");
    gCLExp_toy_2s->GetYaxis()
      ->SetTitle("95% CL limit on X#rightarrowhh event yield");
  }
  else {
    gCLExp_toy->GetYaxis()
      ->SetTitle("95% CL limit on #sigma_{X}#timesBR_{X#rightarrowhh} [pb]");
    gCLObs_toy->GetYaxis()
      ->SetTitle("95% CL limit on #sigma_{X}#timesBR_{X#rightarrowhh} [pb]");
    gCLExp_toy_2s->GetYaxis()
      ->SetTitle("95% CL limit on #sigma_{X}#timesBR_{X#rightarrowhh} [pb]");
  }
  
  gCLExp_toy->SetLineColor(kBlack);
  gCLObs_toy->SetLineColor(kBlack);
  gCLExp_toy->SetLineStyle(2);
  gCLObs_toy->SetLineStyle(1);
  gCLExp_toy->SetLineWidth(2);
  gCLObs_toy->SetLineWidth(2);
  gCLExp_toy_2s->SetFillColor(kYellow);
  gCLExp_toy_1s->SetFillColor(kGreen);

  // Asymptotic graph formatting:
  gCLExp_asym->GetXaxis()->SetTitle("m_{X} [GeV]");
  gCLObs_asym->GetXaxis()->SetTitle("m_{X} [GeV]");
  gCLExp_asym_2s->GetXaxis()->SetTitle("m_{X} [GeV]");
  if (options.Contains("NEvents")) {
    gCLExp_asym->GetYaxis()
      ->SetTitle("95% CL limit on #sigma_{X#rightarrowhh} [pb]");
    gCLObs_asym->GetYaxis()
      ->SetTitle("95% CL limit on #sigma_{X#rightarrowhh} [pb]");
    gCLExp_asym_2s->GetYaxis()
      ->SetTitle("95% CL limit on #sigma_{X#rightarrowhh} [pb]");
  }
  else {
    gCLExp_asym->GetYaxis()
      ->SetTitle("95% CL limit on X#rightarrowhh event yield");
    gCLObs_asym->GetYaxis()
      ->SetTitle("95% CL limit on X#rightarrowhh event yield");
    gCLExp_asym_2s->GetYaxis()->SetTitle("95% CL limit on event yield");
  }
  gCLExp_asym->SetLineColor(kBlack);
  gCLObs_asym->SetLineColor(kBlack);
  gCLExp_asym->SetLineStyle(2);
  gCLObs_asym->SetLineStyle(1);
  gCLExp_asym->SetLineWidth(2);
  gCLObs_asym->SetLineWidth(2);
  gCLExp_asym_2s->SetFillColor(kYellow);
  gCLExp_asym_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.61,0.68,0.89,0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!config->getBool("DoBlind")) leg.AddEntry(gCLObs_toy,"Obs. limit","l");
  leg.AddEntry(gCLExp_toy,"Exp. limit","l");
  leg.AddEntry(gCLExp_toy_1s,"Exp. limit #pm1#sigma_{exp}","F");
  leg.AddEntry(gCLExp_toy_2s,"Exp. limit #pm2#sigma_{exp}","F");
  
  // Plotting options:
  if (options.Contains("toy")) {
    //gCLExp_toy_2s->GetYaxis()->SetRangeUser(1.8, 24);
    //gCLExp_toy_2s->GetYaxis()->SetRangeUser(1.8, 12);
    gCLExp_toy_2s->Draw("A3");
    gCLExp_toy->Draw("Lsame");
    gCLExp_toy_1s->Draw("3same");
    gCLExp_toy->Draw("LSAME");
    if (!config->getBool("DoBlind")) gCLObs_toy->Draw("LSAME");
  }
  else if (options.Contains("asymptotic")) {
    gCLExp_asym_2s->Draw("A3");
    gCLExp_asym->Draw("Lsame");
    gCLExp_asym_1s->Draw("3same");
    gCLExp_asym->Draw("LSAME");
    if (!config->getBool("DoBlind")) gCLObs_asym->Draw("LSAME");
  }
  else if (options.Contains("both")) {
    gCLExp_toy_2s->Draw("A3");
    gCLExp_toy->Draw("Lsame");
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
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.2, 0.87, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.32, 0.87, config->getStr("ATLASLabel"));
  t.DrawLatex(0.2, 0.81, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			      (config->getNum("AnalysisLuminosity")/1000.0)));
  
  // Print the canvas:
  if (options.Contains("asymptotic")) {
    if (options.Contains("NEvents")) {
      can->Print(Form("%s/limits_asymptotic_NEvent.eps", outputDir.Data()));
      can->Print(Form("%s/limits_asymptotic_NEvent.C", outputDir.Data()));
    }
    else {
      can->Print(Form("%s/limits_asymptotic.eps", outputDir.Data()));
      can->Print(Form("%s/limits_asymptotic.C", outputDir.Data()));
    }
  }
  else if (options.Contains("toy")) {
    if (options.Contains("NEvents")) {
      can->Print(Form("%s/limits_toy_NEvent.eps", outputDir.Data()));
      can->Print(Form("%s/limits_toy_NEvent.C", outputDir.Data()));
    }
    else {
      can->Print(Form("%s/limits_toy.eps", outputDir.Data()));
      can->Print(Form("%s/limits_toy.C", outputDir.Data()));
    }
  }
  else {
    if (options.Contains("NEvents")) {
      can->Print(Form("%s/limits_comparison_NEvent.eps", outputDir.Data()));
      can->Print(Form("%s/limits_comparison_NEvent.C", outputDir.Data()));
    }
    else {
      can->Print(Form("%s/limits_comparison.eps", outputDir.Data()));
      can->Print(Form("%s/limits_comparison.C", outputDir.Data()));
    }
  }
  
  // Delete pointers, close files, return:
  std::cout << "DHPlotCLvsMX: Finished!" << std::endl;
  delete can;
  delete gCLObs_asym;
  delete gCLExp_asym;
  delete gCLObs_toy;
  delete gCLExp_toy;
  delete gCLExp_asym_2s;
  delete gCLExp_asym_1s;
  delete gCLExp_toy_2s;
  delete gCLExp_toy_1s;
  delete config;
  return 0;
}
