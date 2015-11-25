////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  measureSignalBias.cxx                                                     //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 22/11/2015                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This main method loads a pre-existing signal parameterization and performs//
//  toy MC bias studies.                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamTools/HggTwoSidedCBPdf.h"
#include "HGamTools/SigParam.h"
#include "HGamTools/AtlasStyle.h"

#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"

/**
   -----------------------------------------------------------------------------
   Extract information on the average, minimum, and maximum values from a vector
   of doubles. 
   @param values - The vector of doubles to explore. 
   @param avg - The average variable passed by reference.
   @param min - The minimum variable passed by reference.
   @param max - The maximum variable passed by reference.
*/
void getAvgMinMax(std::vector<double> values, double& avg, double& min, 
		  double& max) {
  avg = 0.0;
  for (int i_e = 0; i_e < (int)values.size(); i_e++) {
    avg += values[i_e];
    if (values[i_e] < min || i_e == 0) min = values[i_e];
    if (values[i_e] > max || i_e == 0) max = values[i_e];
  }
  avg = avg / (double)values.size();
}

/**
   -----------------------------------------------------------------------------
   Create a LaTex formatted version of the variable name for printing:
*/
TString varPrintName(TString varName) {
  varName.ReplaceAll("frac", "fraction_{");
  varName.ReplaceAll("Nom","");
  varName.ReplaceAll("nCB","n_{CB");
  varName.ReplaceAll("width","w_{");
  varName.ReplaceAll("sigma","#sigma_{");
  varName.ReplaceAll("alpha","#alpha_{");
  varName.ReplaceAll("mu", "#mu_{");
  varName += "}";
  return varName;
}

/**
   -----------------------------------------------------------------------------
   The main method for this utility. Provide 1 argument - the location of the 
   config (.cfg) file, which should be stored in the data/ directory.
*/
int main(int argc, char *argv[])
{
  // Check that the config file location is provided.
  if (argc < 2) HG::fatal("No arguemnts provided");
  HG::Config *settings = new HG::Config(TString(argv[1]));

  // Print configuration for benefit of user:
  std::cout << "measureSignalBias will run with parameters:"
	    << std::endl;
  settings->printDB();
  
  // Set the function type:
  TString function = settings->getStr("SignalFunctionalForm");
  
  // Check that output directory exists:
  TString outputDir = settings->getStr("IODirectory");
  system(Form("mkdir -vp %s/Bias", outputDir.Data()));
  
  // Set the ATLAS Style for plots:
  SetAtlasStyle();
  
  // Choose the category for bias studies:
  int category = settings->getInt("Category");
  
  //Settings for toy Monte Carlo:
  int seed = settings->getInt("RandomSeed");
  
  // Create canvas for plotting:
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  
  // Get the mass point to investigate:
  double mass = settings->getNum("ResonanceMass");
  
  // Store normalization from fit and data:
  std::map<TString,std::vector<double> > biasVals; biasVals.clear();
  biasVals["MassFit"].clear();
  biasVals["MassData"].clear();
  biasVals["NormFit"].clear();
  biasVals["NormData"].clear();
  
  // Now get bias on a point.
  for (int i_t = 0; i_t < settings->getInt("NumberOfToys"); i_t++) {
    SigParam *spi
      = new SigParam(settings->getStr("SampleName"),outputDir+"/Parameterized");
    // Load the signal parameterization from file each time:
    spi->loadParameterization(outputDir+"/Parameterized",
			      settings->getStr("SampleName"));
    // Then toss and fit toy, storing bias data for each:
    std::vector<double> currVals
      = spi->doBiasTest(mass, category, "mctoy", seed+i_t);
    if ((int)currVals.size() >= 4) {
      if (currVals[2] > 0.0001) {
	biasVals["MassData"].push_back(currVals[0]);
	biasVals["MassFit"].push_back(currVals[1]);
	biasVals["NormData"].push_back(currVals[2]);
	biasVals["NormFit"].push_back(currVals[3]);
      }
    }
    delete spi;
  }
  
  // Create output file to store the bias measurements:
  std::ofstream biasOutput; 
  biasOutput.open(Form("%s/Bias/biasValues_Mass%2.2f.txt",
		       outputDir.Data(), mass));
  
  // Loop over the measurements (mass and normalization):
  std::vector<TString> names; names.clear();
  names.push_back("Mass");
  names.push_back("Norm");
  for (int i_n = 0; i_n < (int)names.size(); i_n++) {
    
    // Then plot the bias on the normalization:
    double fitAvg = 0.0; double fitMin = 0.0; double fitMax = 0.0;
    double dataAvg = 0.0; double dataMin = 0.0; double dataMax = 0.0;
    getAvgMinMax(biasVals[Form("%sFit",names[i_n].Data())],
		 fitAvg, fitMin, fitMax);
    getAvgMinMax(biasVals[Form("%sData",names[i_n].Data())], 
		 dataAvg, dataMin, dataMax);
    
    double bias = 100 * ((fitAvg - dataAvg) / dataAvg);
    
    std::cout << names[i_n] << "\t "
	      << dataAvg << "\t " << dataMin << "\t " << dataMax
	      << "\t " << fitAvg << "\t " << fitMin << "\t " << fitMax
	      << "\t " << bias << std::endl;
    biasOutput << names[i_n] << "\t " << fitAvg << "\t " << dataAvg
	       << "\t " << bias << std::endl;
    
    TH1F *hFit = new TH1F("hFit", "hFit", 50, fitMin, fitMax);
    TH1F *hData = new TH1F("hData", "hData", 50, fitMin, fitMax);
    
    for (int i_e = 0; 
	 i_e < (int)(biasVals[Form("%sFit",(names[i_n]).Data())]).size();
	 i_e++) {
      hFit->Fill(biasVals[Form("%sFit",(names[i_n]).Data())][i_e]);
      hData->Fill(biasVals[Form("%sData",(names[i_n]).Data())][i_e]);
    }
    hFit->GetXaxis()->SetTitle(names[i_n]);
    hFit->SetFillColor(kBlue-7);
    hFit->SetFillStyle(3345);
    hFit->SetLineColor(kBlue+1);
    hFit->GetYaxis()->SetRangeUser(0, (1.3 * hFit->GetMaximum()));
    hData->SetFillColor(kRed-7);
    hData->SetFillStyle(3354);
    hData->SetLineColor(kRed+1);
    can->cd(); 
    can->Clear();
    hFit->Draw("hist");
    if (!names[i_n].Contains("Mass")) hData->Draw("histSAME");
    
    // Draw average value:
    TLine *lFitAvg = new TLine();
    lFitAvg->SetLineStyle(2);
    lFitAvg->SetLineWidth(2);
    lFitAvg->SetLineColor(kBlue+1);
    lFitAvg->DrawLine(fitAvg, 0, fitAvg, hFit->GetMaximum());
    // Draw origin value:
    TLine *lDataAvg = new TLine();
    lDataAvg->SetLineStyle(2);
    lDataAvg->SetLineWidth(2);
    lDataAvg->SetLineColor(kRed+1);
    lDataAvg->DrawLine(dataAvg, 0, dataAvg, hFit->GetMaximum());
    
    TLegend leg(0.55, 0.78, 0.89, 0.92);
    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(hFit, "Fit Distribution", "F");
    leg.AddEntry(lFitAvg, Form("Fit Mean = %2.3f", fitAvg), "L");
    if (!names[i_n].Contains("Mass")) {
      leg.AddEntry(hData, "Truth Distribution", "F");
    }
    leg.AddEntry(lDataAvg, Form("Truth Value = %2.3f", dataAvg), "L");
    leg.Draw("SAME");
    
    can->Print(Form("%s/Bias/biasDist_%s_Mass%2.2f.eps",
		    outputDir.Data(), names[i_n].Data(), mass));
    delete hFit;
    delete hData;
    delete lFitAvg;
    delete lDataAvg;
  }

  biasOutput.close();
  delete can;
  return 0;
}
