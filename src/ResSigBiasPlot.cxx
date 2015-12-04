////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  MakeBiasPlot.cxx                                                          //
//                                                                            //
//  Created: Andrew Hard 02/12/2015                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// C++ libraries:
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

// ROOT libraries:
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"

// Package libraries:
#include "CommonFunc.h"

/**
   -----------------------------------------------------------------------------
   Main method creates a plot of biases from signal modeling.
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc != 1) {
    std::cout << "\nUsage: " << argv[0] << std::endl;
    exit(0);
  }
  CommonFunc::SetAtlasStyle();
  
  double masses[4] = {275, 300, 325, 350};
  
  // Dummy for axis ranges (since TGraphs are impossible to rescale...):
  double mass_dummy[4] = {-1.0, 1.0, 2.0, 3.0};
  //double norm_dummy[4] = {-1.0, 1.0, 3.0, 6.0};
  double norm_dummy[4] = {-1.0, 1.0, 2.0, 3.0};
  
  // CBGA:
  double mass_CBGA[4] = {0.2109, 0.0713, 0.1960, 0.3923};
  double norm_CBGA[4] = {-0.0630, -0.0648, -0.0571, -0.0289};
  
  // DoubleCB:
  double mass_DoubleCB[4] = {-0.023774, 0.0587716, 0.451878, 0.0847704};
  double norm_DoubleCB[4] = {-0.0828842, -0.0868228, -0.0903737, -0.0546082};
  
  // CBPlusVoigt:
  double mass_CBPlusVoigt[4] = {2.5914, 2.71807, 1.88794, 1.72865};
  double norm_CBPlusVoigt[4] = {-0.0367975, 5.67886, -0.0348815, 4.60176};
  
  // GAx3:
  double mass_GAx3[4] = {0.131433, 0.0268657, 0.0548338, 0.0359159};
  double norm_GAx3[4] = {0.0878434, -0.0345971, -0.197138, -0.00965631};
  
  // Create graphs:
  TGraph *gMass_dummy = new TGraph(4, masses, mass_dummy);
  TGraph *gMass_CBGA = new TGraph(4, masses, mass_CBGA);
  TGraph *gMass_DoubleCB = new TGraph(4, masses, mass_DoubleCB);
  TGraph *gMass_CBPlusVoigt = new TGraph(4, masses, mass_CBPlusVoigt);
  TGraph *gMass_GAx3 = new TGraph(4, masses, mass_GAx3);

  TGraph *gNorm_dummy = new TGraph(4, masses, norm_dummy);
  TGraph *gNorm_CBGA = new TGraph(4, masses, norm_CBGA);
  TGraph *gNorm_DoubleCB = new TGraph(4, masses, norm_DoubleCB);
  TGraph *gNorm_CBPlusVoigt = new TGraph(4, masses, norm_CBPlusVoigt);
  TGraph *gNorm_GAx3 = new TGraph(4, masses, norm_GAx3);
  
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  can->cd();
  
  // Format histograms:
  gMass_dummy->SetLineColor(0);
  gMass_CBGA->SetLineColor(kRed+1);
  gMass_DoubleCB->SetLineColor(kBlue+1);
  gMass_CBPlusVoigt->SetLineColor(kGreen+1);
  gMass_GAx3->SetLineColor(kOrange+1);
  gNorm_dummy->SetLineColor(0);
  gNorm_CBGA->SetLineColor(kRed+1);
  gNorm_DoubleCB->SetLineColor(kBlue+1);
  gNorm_CBPlusVoigt->SetLineColor(kGreen+1);
  gNorm_GAx3->SetLineColor(kOrange+1);
  
  gMass_CBGA->SetLineWidth(2);
  gMass_DoubleCB->SetLineWidth(2);
  gMass_CBPlusVoigt->SetLineWidth(2);
  gMass_GAx3->SetLineWidth(2);
  gNorm_CBGA->SetLineWidth(2);
  gNorm_DoubleCB->SetLineWidth(2);
  gNorm_CBPlusVoigt->SetLineWidth(2);
  gNorm_GAx3->SetLineWidth(2);
  
  gMass_dummy->GetXaxis()->SetTitle("M_{X} [GeV]");
  gNorm_dummy->GetXaxis()->SetTitle("M_{X} [GeV]");
  
  gMass_dummy->GetYaxis()->SetTitle("Mass Bias [%]");
  gNorm_dummy->GetYaxis()->SetTitle("Normalization Bias [%]");
  
  TLegend leg(0.61, 0.74, 0.9, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  //leg.SetTextSize(0.04);
  leg.AddEntry(gMass_CBGA, "Crystal Ball + Gaussian", "L");
  leg.AddEntry(gMass_DoubleCB, "Double-sided Crystal Ball", "L");
  leg.AddEntry(gMass_CBPlusVoigt, "Crystal Ball + Voigt", "L");
  leg.AddEntry(gMass_GAx3, "Triple Gaussian", "L");
  
  TLine *line = new TLine();
  line->SetLineStyle(3);
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
  
  // Plot mass biases:
  gMass_dummy->Draw("AL");
  line->DrawLine(gMass_dummy->GetXaxis()->GetXmin(), 0.0, 
		 gMass_dummy->GetXaxis()->GetXmax(), 0.0);
  gMass_CBPlusVoigt->Draw("LSAME");
  gMass_CBGA->Draw("LSAME");
  gMass_DoubleCB->Draw("LSAME");
  gMass_GAx3->Draw("LSAME");
  leg.Draw("SAME");
  can->Print("mass_bias.eps");
  can->Clear();
  
  // Plot normalization biases:
  gNorm_dummy->Draw("AL");
  line->DrawLine(gMass_dummy->GetXaxis()->GetXmin(), 0.0, 
		 gMass_dummy->GetXaxis()->GetXmax(), 0.0);
  gNorm_CBPlusVoigt->Draw("LSAME");
  gNorm_CBGA->Draw("LSAME");
  gNorm_DoubleCB->Draw("LSAME");
  gNorm_GAx3->Draw("LSAME");
  leg.Draw("SAME");
  can->Print("norm_bias.eps");
  
  return 0;
}
