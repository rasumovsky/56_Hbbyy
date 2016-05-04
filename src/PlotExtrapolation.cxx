////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PlotExtrapolation.cxx                                               //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 04/05/2016                                                          //
//                                                                            //
//  Plot some extrapolated test statistics for the 2016 dataset.              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// C++ libraries:
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

// ROOT libraries:
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TROOT.h"
#include "TString.h"
#include "Rtypes.h"

// Package libraries:
#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
  
/**
   -----------------------------------------------------------------------------
   Produces a plot of the Myy cut efficiency and fits it for each category.
   @param configFile - The analysis configuration file.
   @param options - The plot options.
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 2) { 
    std::cout << "\nUsage: " << argv[0] << "<configFile>" << std::endl;
    exit(0);
  }
  
  TString configFileName = argv[1];
    
  // Set the ATLAS Style:
  CommonFunc::SetAtlasStyle();
  
  // Create config class:
  Config *m_config = new Config(configFileName);
  TString outputDir = Form("%s/%s/PlotExtrapolation", 
			   (m_config->getStr("MasterOutput")).Data(), 
			   (m_config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  /*
  double zXsec[6] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};
  double z0_med[6] = {1.177, 1.512, 1.873, 2.452, 2.996, 4.923};
  double z0_P1[6] = {0.373, 0.508, 0.673, 0.807, 1.11, 1.69};
  double z0_N1[6] = {0.482, 0.650, 0.846, 0.979, 1.346, 2.10};
  */  
  double zXsec[5] = {0.5, 1.0, 2.0, 4.0, 16.0};
  double z0_med[5] = {1.177, 1.512, 1.873, 2.452, 4.923};
  double z0_P1[5] = {0.373, 0.508, 0.673, 0.807, 1.69};
  double z0_N1[5] = {0.482, 0.650, 0.846, 0.979, 2.10};
  
  TGraph *zMed = new TGraph(5, zXsec, z0_med);
  zMed->SetLineColor(kBlue+1);
  TGraphAsymmErrors *zBand = new TGraphAsymmErrors(5, zXsec, z0_med, 0, 0, 
						   z0_N1, z0_P1);
  zBand->GetXaxis()->SetTitle("luminosity in 2016 [fb^{-1}]");
  zBand->GetYaxis()->SetTitle("Z_{0} [#sigma]");
  zBand->SetFillColor(kBlue-10);
  /* 
  double clXsec[6] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};
  double CL_med[6] = {0.552135, 0.710045, 0.844945, 0.891268, 0.971635, 0.9992};
  double CL_P1[6] = {0.101, 0.1005, 0.0792, 0.0738, 0.0240, 0.000741};
  double CL_N1[6] = {0.158, 0.183, 0.190, 0.258, 0.183, 0.0276};
  */

  double clXsec[5] = {0.5, 1.0, 2.0, 8.0, 16.0};
  double CL_med[5] = {0.552135, 0.710045, 0.844945, 0.971635, 0.9992};
  double CL_P1[5] = {0.101, 0.1005, 0.0792, 0.0240, 0.000741};
  double CL_N1[5] = {0.158, 0.183, 0.190, 0.183, 0.0276};

  TGraph *clMed = new TGraph(5, clXsec, CL_med);
  clMed->SetLineColor(kRed+1);
  TGraphAsymmErrors *clBand = new TGraphAsymmErrors(5, clXsec, CL_med, 0, 0, 
						    CL_N1, CL_P1);
  clBand->GetXaxis()->SetTitle("Luminosity in 2016 [fb^{-1}]");
  clBand->GetYaxis()->SetTitle("CL_{S} exclusion");
  clBand->GetYaxis()->SetRangeUser(0.4, 1.0);
  clBand->SetFillColor(kRed-10);
  

  TCanvas *can = new TCanvas("can", "can");
  can->cd();
  
  TLegend zLeg(0.6, 0.2, 0.9, 0.36);
  zLeg.SetBorderSize(0);
  zLeg.SetFillColor(0);
  zLeg.SetTextSize(0.04);
  zLeg.AddEntry(zMed, "#hat{#sigma}_{X#rightarrowhh}^{#sqrt{s}=8TeV}", "L");
  zLeg.AddEntry(zBand, "#hat{#sigma}_{X#rightarrowhh}^{#sqrt{s}=8TeV} #pm #Delta#hat{#sigma}_{X#rightarrowhh}^{#sqrt{s}=8TeV}", "F");
  
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->SetLineColor(kRed+1);
  
  zBand->Draw("A3");
  zMed->Draw("LSAME");
  line->DrawLine(zBand->GetXaxis()->GetXmin(), 5.0,
		 zBand->GetXaxis()->GetXmax(), 5.0);
  zLeg.Draw("SAME");
  can->Print(Form("%s/Extrapolate_z0.eps", outputDir.Data()));
  
  TLegend clLeg(0.6, 0.2, 0.9, 0.36);
  clLeg.SetBorderSize(0);
  clLeg.SetFillColor(0);
  clLeg.SetTextSize(0.04);
  clLeg.AddEntry(clMed, "#hat{#sigma}_{X#rightarrowhh}^{#sqrt{s}=8TeV}", "L");
  clLeg.AddEntry(clBand, "#hat{#sigma}_{X#rightarrowhh}^{#sqrt{s}=8TeV} #pm #Delta#hat{#sigma}_{X#rightarrowhh}^{#sqrt{s}=8TeV}", "F");
  
  clBand->Draw("A3");
  clMed->Draw("LSAME");
  line->DrawLine(zBand->GetXaxis()->GetXmin(), 0.95,
		 zBand->GetXaxis()->GetXmax(), 0.95);
  clLeg.Draw("SAME");
  can->Print(Form("%s/Extrapolate_CL.eps", outputDir.Data()));
  
  return 0;
}
