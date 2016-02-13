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
//  Options: Can be "toy" or "asymptotic" or "both"                           //
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
      
      ifstream inputFile_asym(Form("%s/limits_asymptotic_ResonantMX%d.txt",
				   inputDir.Data(), mX));
      if (inputFile_asym.is_open()) {
	TString inName = ""; double inValue = 0.0;
	while (inputFile_asym >> inName >> inValue) {
	  inValue = inValue / 1000.0;
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
      
      ifstream inputFile_toy(Form("%s/limits_toy_ResonantMX%d.txt",
				  inputDir.Data(), mX));
      if (inputFile_toy.is_open()) {
	TString inName = ""; double inValue = 0.0;
	while (inputFile_toy >> inName >> inValue) {
	  inValue = inValue / 1000.0;
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
  gCLExp_toy->GetXaxis()->SetTitle("M_{X} [GeV]");
  gCLObs_toy->GetXaxis()->SetTitle("M_{X} [GeV]");
  gCLExp_toy_2s->GetXaxis()->SetTitle("M_{X} [GeV]");
  gCLExp_toy->GetYaxis()->SetTitle("95% CL #sigma_{X#rightarrowhh} [pb]");
  gCLObs_toy->GetYaxis()->SetTitle("95% CL #sigma_{X#rightarrowhh} [pb]");
  gCLExp_toy_2s->GetYaxis()->SetTitle("95% CL #sigma_{X#rightarrowhh} [pb]");
  gCLExp_toy->SetLineColor(kBlack);
  gCLObs_toy->SetLineColor(kBlack);
  gCLExp_toy->SetLineStyle(2);
  gCLObs_toy->SetLineStyle(1);
  gCLExp_toy->SetLineWidth(2);
  gCLObs_toy->SetLineWidth(2);
  gCLExp_toy_2s->SetFillColor(kYellow);
  gCLExp_toy_1s->SetFillColor(kGreen);

  // Asymptotic graph formatting:
  gCLExp_asym->GetXaxis()->SetTitle("M_{X} [GeV]");
  gCLObs_asym->GetXaxis()->SetTitle("M_{X} [GeV]");
  gCLExp_asym_2s->GetXaxis()->SetTitle("M_{X} [GeV]");
  gCLExp_asym->GetYaxis()->SetTitle("95% CL #sigma_{X#rightarrowhh} [pb]");
  gCLObs_asym->GetYaxis()->SetTitle("95% CL #sigma_{X#rightarrowhh} [pb]");
  gCLExp_asym_2s->GetYaxis()->SetTitle("95% CL #sigma_{X#rightarrowhh} [pb]");
  gCLExp_asym->SetLineColor(kBlack);
  gCLObs_asym->SetLineColor(kBlack);
  gCLExp_asym->SetLineStyle(2);
  gCLObs_asym->SetLineStyle(1);
  gCLExp_asym->SetLineWidth(2);
  gCLObs_asym->SetLineWidth(2);
  gCLExp_asym_2s->SetFillColor(kYellow);
  gCLExp_asym_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.60,0.73,0.88,0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!config->getBool("DoBlind")) leg.AddEntry(gCLObs_toy,"Obs. Limit","l");
  leg.AddEntry(gCLExp_toy,"Exp. Limit","l");
  leg.AddEntry(gCLExp_toy_1s,"Exp. Limit #pm1#sigma_{exp}","F");
  leg.AddEntry(gCLExp_toy_2s,"Exp. Limit #pm2#sigma_{exp}","F");
  
  // Plotting options:
  if (options.Contains("toy")) {
    //gCLExp_toy->Draw("AL");
    //gCLExp_toy_2s->Draw("3same");
    gCLExp_toy_2s->Draw("A3");
    gCLExp_toy->Draw("Lsame");
    gCLExp_toy_1s->Draw("3same");
    gCLExp_toy->Draw("LSAME");
    if (!config->getBool("DoBlind")) gCLObs_toy->Draw("LSAME");
  }
  else if (options.Contains("asymptotic")) {
    //gCLExp_asym->Draw("AL");
    //gCLExp_asym_2s->Draw("3same");
    gCLExp_asym_2s->Draw("A3");
    gCLExp_asym->Draw("Lsame");
    gCLExp_asym_1s->Draw("3same");
    gCLExp_asym->Draw("LSAME");
    if (!config->getBool("DoBlind")) gCLObs_asym->Draw("LSAME");
  }
  else if (options.Contains("both")) {
    //gCLExp_toy->Draw("AL");
    //gCLExp_toy_2s->Draw("3same");
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
      leg.AddEntry(gCLObs_asym,"Obs. Limit (asymptotic)","l");
    }
    leg.AddEntry(gCLExp_asym,"Exp. Limit (asymptotic)","l");
  }
  gPad->RedrawAxis();
  leg.Draw("SAME");
    
  // Print the canvas:
  if (options.Contains("asymptotic")) {
    can->Print(Form("%s/limits_asymptotic.eps", outputDir.Data()));
  }
  else if (options.Contains("toy")) {
    can->Print(Form("%s/limits_toy.eps", outputDir.Data()));
  }
  else can->Print(Form("%s/limits_comparison.eps", outputDir.Data()));
  
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
