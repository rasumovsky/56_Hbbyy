////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PlotMyyCutEff.cxx                                                   //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 18/01/2016                                                          //
//                                                                            //
//  This class fits the efficiency values for the Myy cut in the di-Higgs     //
//  resonant search.                                                          //
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
  TString outputDir = Form("%s/%s/MyyCutEff", 
			   (m_config->getStr("MasterOutput")).Data(), 
			   (m_config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  double massValues[4] = {275, 300, 325, 350};
  double massError[4] = {0.1, 0.1, 0.1, 0.1};
  
  //double value_jj[4] = {0.8327, 0.8446, 0.8818, 0.8846};//MC
  double value_jj[4] = {0.7949, 0.8666, 0.9104, 0.9370};//Fit
  double error_jj[4] = {0.0378, 0.0220, 0.0286, 0.0524};
  
  //double value_bj[4] = {0.8663, 0.8722, 0.8792, 0.8864};//MC
  double value_bj[4] = {0.9041, 0.8791, 0.9036, 0.9222};//Fit
  double error_bj[4] = {0.0378, 0.0069, 0.0244, 0.0358};
  
  //double value_bb[4] = {0.8722, 0.8798, 0.8909, 0.8924};//MC
  double value_bb[4] = {0.8691, 0.8753, 0.9015, 0.8993};//Fit
  double error_bb[4] = {0.0081, 0.0045, 0.0106, 0.0069};
  
  // Create TGraphs:
  TGraphErrors *graph_jj
    = new TGraphErrors(4, massValues, value_jj, massError, error_jj);
  graph_jj->SetMarkerColor(kGreen+1);
  graph_jj->SetLineColor(kGreen+1);
  graph_jj->SetLineWidth(2);
  TGraphErrors *graph_bj
    = new TGraphErrors(4, massValues, value_bj, massError, error_bj);
  graph_bj->SetMarkerColor(kBlue+1);
  graph_bj->SetLineColor(kBlue+1);
  graph_bj->SetLineWidth(2);
  TGraphErrors *graph_bb
    = new TGraphErrors(4, massValues, value_bb, massError, error_bb);
  graph_bb->SetMarkerColor(kRed+1);
  graph_bb->SetLineColor(kRed+1);
  graph_bb->SetLineWidth(2);
  
  // Create fit functions:
  TF1 *fit_jj = new TF1("fit_jj", "pol1", 270, 355);
  fit_jj->SetLineColor(kGreen+1);
  graph_jj->Fit(fit_jj);
  TF1 *fit_bj = new TF1("fit_bj", "pol1", 270, 355);
  fit_bj->SetLineColor(kBlue+1);
  graph_bj->Fit(fit_bj);
  TF1 *fit_bb = new TF1("fit_bb", "pol1", 270, 355);
  fit_bb->SetLineColor(kRed+1);
  graph_bb->Fit(fit_bb);
  
  // Plot the results:
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  can->cd();
  
  graph_jj->GetXaxis()->SetTitle("M_{X} [GeV]");
  graph_jj->GetYaxis()->SetTitle("m_{#gamma#gamma} cut efficiency");
  graph_jj->Draw("AEP");
  graph_bj->Draw("EPSAME");
  graph_bb->Draw("EPSAME");
  
  fit_jj->Draw("LSAME");
  fit_bj->Draw("LSAME");
  fit_bb->Draw("LSAME");
  
  // Create a legend to which we will add items:
  TLegend leg(0.20, 0.75, 0.50, 0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.03);
  leg.AddEntry(graph_jj, "0 b-tag category", "lep");
  leg.AddEntry(graph_bj, "1 b-tag category", "lep");
  leg.AddEntry(graph_bb, "2 b-tag category", "lep");
  leg.Draw("SAME");
  
  TLatex l; l.SetNDC(); l.SetTextColor(kBlack);
  l.SetTextFont(42); l.SetTextSize(0.04);
  l.SetTextColor(kGreen+1);
  l.DrawLatex(0.65, 0.30,
	      Form("#epsilon_{m_{#gamma#gamma}}=%2.2f + %3.4f#upointM_{X}",
		   fit_jj->GetParameter(0), fit_jj->GetParameter(1)));
  l.SetTextColor(kBlue+1);
  l.DrawLatex(0.65, 0.26, 
	      Form("#epsilon_{m_{#gamma#gamma}}=%2.2f + %3.4f#upointM_{X}",
		   fit_bj->GetParameter(0), fit_bj->GetParameter(1)));
  l.SetTextColor(kRed+1);
  l.DrawLatex(0.65, 0.22,
	      Form("#epsilon_{m_{#gamma#gamma}}=%2.2f + %3.4f#upointM_{X}",
		   fit_bb->GetParameter(0), fit_bb->GetParameter(1)));
  
  
  can->Print(Form("%s/MyyCutEffFit.eps", outputDir.Data()));
  return 0;
}
