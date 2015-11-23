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
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Set the ATLAS Style for plots:
  SetAtlasStyle();
  
  // Choose the category for bias studies:
  int category = settings->getInt("Category");
  
  //Settings for toy Monte Carlo:
  int seed = settings->getInt("RandomSeed");
  
  // Create canvas for plotting:
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  
  // Get the mass point to investigate:
  double masses = settings->getNum("ResonanceMass");
  
  // Store parameter names:
  std::vector<TString> pars; pars.clear();
  // Store parameter original values:
  std::vector<double> originVals; originVals.clear();
  // Store data for histograms:
  std::map<TString,std::vector<double> > storeVals; storeVals.clear();
  
  // Load the parameterized fit:
  SigParam *spi
    = new SigParam(settings->getStr("SampleName"),outputDir+"/Parameterized");
  // Load the signal parameterization from file each time:
  spi->loadParameterization(outputDir+"/Parameterized",
			    settings->getStr("SampleName"));
  
  // Now get bias on a point.
  for (int i_t = 0; i_t < settings->getInt("NumberOfToys"); i_t++) {
    SigParam *spt
      = new SigParam(settings->getStr("SampleName"), outputDir+"/Individual");
    // Load the signal parameterization from file each time:
    spt->loadParameterization(outputDir+"/Individual",
			      settings->getStr("SampleName"));
    // The first time, load parameter information and also initial values:
    if (i_t == 0) {
      pars =spt->variablesForFunction(settings->getStr("SignalFunctionalForm"));
      for (int i_p = 0; i_p < (int)pars.size(); i_p++) {
	originVals
	  .push_back(spt->getParameterValue(pars[i_p], mass, category));
	//originVals
	//.push_back(spi->getParameterizedValue(pars[i_p], mass, category));
	storeVals[pars[i_p]].clear();
      }
    }
    
    // Then toss and fit toy:
    if (spt->generateAndFitData(300, 2, "mctoy", seed+i_t)) {
      // Retrieve parameter valuess for the fit to this toy:
      for (int i_p = 0; i_p < (int)pars.size(); i_p++) {
	storeVals[pars[i_p]]
	  .push_back(spt->getParameterValue(pars[i_p], mass, category));
      }
    }
    delete spt;
  }
  
  // Store the average, minimum, and maximum value of each parameter:
  double avg[100] = {0.0};
  double min[100] = {0.0};
  double max[100] = {0.0};
  
  // Loop over variables and extract average, minimum, maximum:
  std::ofstream biasOutput; 
  biasOutput.open(Form("%s/biasValues_Mass%2.2f.txt",
		       outputDir.Data(), mass));
  
  for (int i_p = 0; i_p < (int)pars.size(); i_p++) {
    getAvgMinMax(storeVals[pars[i_p]], avg[i_p], min[i_p], max[i_p]);
    double bias = 100 * ((originVals[i_p] - avg[i_p]) / avg[i_p]);
    std::cout << "Var= " << pars[i_p] << ", avg=" << avg[i_p] << ", origin=" 
	      << originVals[i_p] << ", bias=" << bias << "%" << std::endl;
    biasOutput << "Var+ " << pars[i_p] << ", avg=" << avg[i_p] << ", origin=" 
	       << originVals[i_p] << ", bias=" << bias << "%" << std::endl;
  }
  biasOutput.close();
  
  // Create and plot histograms for each variable distribution:
  for (int i_p = 0; i_p < (int)pars.size(); i_p++) {
    TH1F *currHist = new TH1F(pars[i_p], pars[i_p], 50, min[i_p], max[i_p]);
    for (int i_e = 0; i_e < (int)storeVals[pars[i_p]].size(); i_e++) {
      currHist->Fill(storeVals[pars[i_p]][i_e]);
    }
    currHist->GetXaxis()
      ->SetTitle(Form("%s value",(varPrintName(pars[i_p])).Data()));
    currHist->SetFillColor(kBlue-10);
    currHist->SetLineColor(kBlue+4);
    can->cd(); 
    can->Clear();
    currHist->Draw("hist");
    
    // Draw average value:
    TLine *lineAvg = new TLine();
    lineAvg->SetLineStyle(2);
    lineAvg->SetLineWidth(2);
    lineAvg->SetLineColor(4);
    lineAvg->DrawLine(avg[i_p], 0, avg[i_p], currHist->GetMaximum());
    // Draw origin value:
    TLine *lineFit = new TLine();
    lineFit->SetLineStyle(2);
    lineFit->SetLineWidth(2);
    lineFit->SetLineColor(2);
    lineFit->DrawLine(originVals[i_p], 0, originVals[i_p],
		      currHist->GetMaximum());
    
    TLegend leg(0.55, 0.82, 0.89, 0.92);
    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(lineAvg,Form("Toy Mean = %2.3f", avg[i_p]), "L");
    leg.AddEntry(lineFit,Form("Fit val = %2.3f", originVals[i_p]), "L");
    leg.Draw("SAME");
    
    // Print the canvas:
    can->Print(Form("%s/biasDist_%s.eps",outputDir.Data(),pars[i_p].Data()));
  }
  delete spi;
  delete can;
  return 0;
}

