////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DMToyAnalysis.cxx                                                         //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 24/06/2015                                                          //
//                                                                            //
//  This program compares test statistic values from pseudo-experiment        //
//  ensembles and asymptotic formulae.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// statistic == "QMu", "QMuTilde" 

#include "DMToyAnalysis.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the DMToyAnalysis class.
*/
DMToyAnalysis::DMToyAnalysis(TString newJobName, TString newDMSignal,
			     TString newCateScheme, TString newOptions) {
  jobName = newJobName;
  DMSignal = newDMSignal;
  cateScheme = newCateScheme;
  options = newOptions;
  
  // set input and output directories:
  outputDir = Form("%s/%s/DMToyAnalysis", DMAnalysis::masterOutput.Data(),
		   jobName.Data());
  TString toyDir = Form("%s/%s/DMPseudoExp", DMAnalysis::masterOutput.Data(),
			jobName.Data());
  TString wsFileName = Form("%s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root",
			    DMAnalysis::masterOutput.Data(), jobName.Data(), 
			    DMSignal.Data());
  
  // Set the internal (private) variable initial conditions:
  namesGlobs.clear();
  namesNuisParams.clear();
  nGlobs = 0;
  nNuisParams = 0;
  nBins = 500;
  binMin = 0;
  binMax = 20;
  
  // Create output directory:
  system(Form("mkdir -vp %s",outputDir.Data()));
  
  // Set ATLAS style template:
  CommonFunc::SetAtlasStyle();
  
  // Add all of the individual pseudoexperiment files together:
  TString toyFileMu0 = Form("%s/toy_mu0.root",toyDir.Data());
  TString toyFileMu1 = Form("%s/toy_mu1.root",toyDir.Data());
  system(Form("hadd -f %s %s/single_files/toy_mu0*",
	      toyFileMu0.Data(), toyDir.Data()));
  system(Form("hadd -f %s %s/single_files/toy_mu1*", 
	      toyFileMu1.Data(), toyDir.Data()));
  
  // Create temporary file lists:
  TString listMu0 = "temp_list_mu0.txt";
  TString listMu1 = "temp_list_mu1.txt";
  system(Form("echo %s | tee %s", toyFileMu0.Data(), listMu0.Data()));
  system(Form("echo %s | tee %s", toyFileMu1.Data(), listMu1.Data()));
  
  // Create TChain of ntuples:
  TChain* chainMu0 = CommonFunc::MakeChain("toy", listMu0, "badfile");
  TChain* chainMu1 = CommonFunc::MakeChain("toy", listMu1, "badfile");
  DMToyTree *treeMu0 = new DMToyTree(chainMu0);
  DMToyTree *treeMu1 = new DMToyTree(chainMu1);
  
  // Store the toy data:
  fillToyHistograms(0, treeMu0);
  fillToyHistograms(1, treeMu1);
  
  // Get the Asimov form of the test statistic:
  TFile workspaceFile(wsFileName, "read");
  workspace = (RooWorkspace*)workspaceFile.Get("combinedWS");
  dmts = new DMTestStat(jobName, DMSignal, cateScheme, "new", workspace);
  
  // Get the asymptotic test statistic distribution:
  getAsymptoticForm("QMu");// THIS SHOULD BE GENERALIZED!!!
  
  // Plot the results:
  plotProfiledMu();
  plotTestStat("QMu");
  plotTestStat("Q0");
  plotTestStatComparison("QMu");
  plotTestStatComparison("Q0");
  
  // Remove the temporary file lists:
  system(Form("rm %s",listMu0.Data()));
  system(Form("rm %s",listMu1.Data()));
  
  std::cout << "DMToyAnalysis: Finished!" << std::endl;
  std::cout << "\t" << treeMu0->fChain->GetEntries()
	    << " mu=0 pseudo experiments were analyzed" << std::endl;
  std::cout << "\t" << treeMu1->fChain->GetEntries()
	    << " mu=1 pseudo experiments were analyzed" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Fill the histograms containing toy Data.
   @param muValue - the mu hypothesis under which the toys were generated.
   @param toyTree - the TTree containing the pseudo data.
*/
void DMToyAnalysis::fillToyHistograms(int muValue, DMToyTree *toyTree) {
  int nEvents = toyTree->fChain->GetEntries();
  std::cout << "DMToyAnalysis: Loop over " << nEvents
	    << " pseudoexperiments with mu = " << muValue << std::endl;
  
  // Instantiate the histograms:
  hMuProfiled[muValue] = new TH1F(Form("hMuProfiled%d",muValue),
				  Form("hMuProfiled%d",muValue), 
				  nBins, -2.0, 4.0);
  hQMu[muValue] = new TH1F(Form("hQMu%d",muValue),Form("hQMu%d",muValue),
			   nBins, binMin, binMax);
  hQ0[muValue] = new TH1F(Form("hQ0%d",muValue),Form("hQ0%d",muValue),
			  nBins, binMin, binMax);
  
  for (int i_p = 0; i_p < 20; i_p++) {
    hNuisMu0[i_p][muValue] = new TH1F(Form("hNuisMu0_%d",muValue),
				      Form("hNuisMu0_%d",muValue), 100, -5, 5);
    hNuisMu1[i_p][muValue] = new TH1F(Form("hNuisMu1_%d",muValue),
				      Form("hNuisMu1_%d",muValue), 100, -5, 5);
    hNuisMuFree[i_p][muValue] = new TH1F(Form("hNuisMuFree_%d",muValue),
					 Form("hNuisMuFree_%d",muValue),
					 100, -5, 5);
    
    hGlobsMu0[i_p][muValue] = new TH1F(Form("hGlobsMu0_%d",muValue),
				       Form("hGlobsMu0_%d",muValue),100,-5,5);
    hGlobsMu1[i_p][muValue] = new TH1F(Form("hGlobsMu1_%d",muValue),
				       Form("hGlobsMu1_%d",muValue),100,-5,5);
    hGlobsMuFree[i_p][muValue] = new TH1F(Form("hGlobsMuFree_%d",muValue),
					  Form("hGlobsMuFree_%d",muValue),
					  100, -5, 5);
  }
  
  // Loop over events in the TTree:
  bool isFirstLoop = true;
  for (int i_e = 0; i_e < nEvents; i_e++) {
    toyTree->fChain->GetEntry(i_e);
    
    // Only plot successful fits:
    if (!(toyTree->convergedMu0 && toyTree->convergedMu1 &&
	  toyTree->convergedMuFree)) continue;
    
    // Get the test statistic values:
    double valueQMu = dmts->getQMuFromNLL(toyTree->nllMu1, toyTree->nllMuFree, 
					  toyTree->muDMVal, 1);
    double valueQ0 = dmts->getQ0FromNLL(toyTree->nllMu0, toyTree->nllMuFree,
					toyTree->muDMVal);
    
    // Fill histograms for the test statistics and POI:
    hQMu[muValue]->Fill(valueQMu);
    hQ0[muValue]->Fill(valueQ0);
    hMuProfiled[muValue]->Fill(toyTree->muDMVal);
    
    // Fill the nuisance parameter histograms:
    if (isFirstLoop) nNuisParams = (int)((*toyTree->namesNP).size());
    // loop over the nuisance parameters in the tree:
    for (int i_p = 0; i_p < nNuisParams; i_p++) {
      hNuisMu0[i_p][muValue]->Fill((*toyTree->valuesNPMu0)[i_p] );
      hNuisMu1[i_p][muValue]->Fill((*toyTree->valuesNPMu1)[i_p] );
      hNuisMuFree[i_p][muValue]->Fill((*toyTree->valuesNPMuFree)[i_p]);
      if (isFirstLoop) namesNuisParams.push_back((*toyTree->namesNP)[i_p]);
    }
    
    // Fill the global observable histograms:
    if (isFirstLoop) nGlobs = (int)((*toyTree->namesGlobs).size());
    // loop over the nuisance parameters in the tree:
    for (int i_p = 0; i_p < nGlobs; i_p++) {
      hGlobsMu0[i_p][muValue]->Fill((*toyTree->valuesGlobsMu0)[i_p] );
      hGlobsMu1[i_p][muValue]->Fill((*toyTree->valuesGlobsMu1)[i_p] );
      hGlobsMuFree[i_p][muValue]->Fill((*toyTree->valuesGlobsMuFree)[i_p]);
      if (isFirstLoop) namesNuisParams.push_back((*toyTree->namesNP)[i_p]);
    }
    
    isFirstLoop = false;
  }
  
  // Then scale the histograms:
  hQMu[muValue]->Scale(1.0 / hQMu[muValue]->Integral(1, nBins));
  hQ0[muValue]->Scale(1.0 / hQ0[muValue]->Integral(1, nBins));
}

/**
   -----------------------------------------------------------------------------
   Get the asymptotic form of the test statistic. Stored in hAsimov histogram.
   @param statistic - the test statistic for asymptotic formula.
*/
void DMToyAnalysis::getAsymptoticForm(TString statistic) {
  std::cout << "DMToyAnalysis: Constructing Asymptotic form of "
	    << statistic << "." << std::endl;
  
  // The histogram to contain the asymptotic form:
  hAsymptotic = new TH1F("hAsymptotic", "hAsymptotic", nBins, binMin, binMax);
  
  // First get the value from fitting Asimov data (asimovDataMu0):
  double muHat = 0.0;
  double nllMu1 = dmts->getFitNLL("asimovDataMu0", 1.0, true, muHat);
  double nllMu0 = dmts->getFitNLL("asimovDataMu0", 0.0, true, muHat);
  double nllMuHat = dmts->getFitNLL("asimovDataMu0", 0.0, false, muHat);
  double qMu = dmts->getQMuFromNLL(nllMu1, nllMuHat, muHat, 1);
  double qMuTilde = dmts->getQMuTildeFromNLL(nllMu1, nllMu0, nllMuHat, muHat,1);
  
  double asimovTestStat = 0.0;
  if (statistic.EqualTo("QMu")) asimovTestStat = qMu;
  else if (statistic.EqualTo("QMuTilde")) asimovTestStat = qMuTilde;
  
  // Construct the test statistic function:
  for (int i_b = 1; i_b <= 1000; i_b++) {
    double q = hAsymptotic->GetBinCenter(i_b);
    if (statistic.EqualTo("QMu")) {
      hAsymptotic->SetBinContent(i_b, dmts->functionQMu(q));
    }
    else if (statistic.EqualTo("QMuTilde")) {
      hAsymptotic->SetBinContent(i_b, dmts->functionQMuTilde(q,asimovTestStat));
    }
  }
  // Then scale:
  hAsymptotic->SetBinContent(0, 0);
  hAsymptotic->SetBinContent(nBins+1, 0);
  hAsymptotic->Scale(1.0 / hAsymptotic->Integral(1, nBins));
  hAsymptotic->SetBinContent(1,1+hAsymptotic->GetBinContent(1));// WHY?
  hAsymptotic->Scale(1.0 / hAsymptotic->Integral(1,nBins));
}

/**
   -----------------------------------------------------------------------------
   Retrieve the functional form of the asymptotic test statistic approximation.
*/
TH1F* DMToyAnalysis::getAsymptoticHist() {
  return hAsymptotic;
}

/**
   -----------------------------------------------------------------------------
   Get the histogram of a particular global observable.
   @param paramName - the name of the global observable.
   @param fitType - the type of fit generating the distribution.
   @param toyMu - the mu value used to generate the toy data that was fitted.
*/
TH1F* DMToyAnalysis::getGlobsHist(TString paramName, TString fitType,
				  int toyMu) {
  int paramIndex = 0;
  for (paramIndex = 0; paramIndex < (int)namesGlobs.size(); paramIndex++) {
    if (TString(namesGlobs[paramIndex]).Contains(paramName)) break;
  }
  if (fitType.EqualTo("Mu0")) return hGlobsMu0[paramIndex][toyMu];
  else if (fitType.EqualTo("Mu1")) return hGlobsMu1[paramIndex][toyMu];
  else if (fitType.EqualTo("MuFree")) return hGlobsMuFree[paramIndex][toyMu];
  else return NULL;
}

/**
   -----------------------------------------------------------------------------
   Get the histogram of a particular nuisance parameter.
   @param paramName - the name of the nuisance parameter.
   @param fitType - the type of fit generating the distribution.
   @param toyMu - the mu value used to generate the toy data that was fitted.
*/
TH1F* DMToyAnalysis::getNuisHist(TString paramName, TString fitType,
				 int toyMu) {
  int paramIndex = 0;
  for (paramIndex = 0; paramIndex < (int)namesNuisParams.size(); paramIndex++) {
    if (TString(namesNuisParams[paramIndex]).Contains(paramName)) break;
  }
  if (fitType.EqualTo("Mu0")) return hNuisMu0[paramIndex][toyMu];
  else if (fitType.EqualTo("Mu1")) return hNuisMu1[paramIndex][toyMu];
  else if (fitType.EqualTo("MuFree")) return hNuisMuFree[paramIndex][toyMu];
  else return NULL;
}

/**
   -----------------------------------------------------------------------------
   Get the signal strength histogram.
   @param toyMu - the mu value used to generate the toy data that was fitted.
*/
TH1F* DMToyAnalysis::getMuHist(int toyMu) {
  return hMuProfiled[toyMu];
}

/**
   -----------------------------------------------------------------------------
   Get the test statistic histogram.
   @param statistic - the name of the test statistic: Q0, QMu, QMuTilde.
   @param toyMu - the mu value used to generate the toy data that was fitted.
*/
TH1F* DMToyAnalysis::getStatHist(TString statistic, int toyMu) {
  if (statistic.EqualTo("Q0")) return hQ0[toyMu];
  else if (statistic.EqualTo("QMu")) return hQMu[toyMu];
  else if (statistic.EqualTo("QMuTilde")) return hQMuTilde[toyMu];
  else return NULL;
}

/**
   -----------------------------------------------------------------------------
   Plot the distributions of nuisance parameters and global observables
   @param paramName - the name of the global observable or nuisance parameter.
   @param paramType - global observable "Globs" or nuisance parameter "Nuis"
   @param toyMu - the mu value used to generate the toy data that was fitted.
*/
void DMToyAnalysis::plotParameter(TString paramName, TString paramType,
				  int toyMu) {
  
  TH1F *histMu0 = NULL; TH1F *histMu1 = NULL; TH1F *histMuFree = NULL;
  if (paramType.EqualTo("Globs")) {
    histMu0 = getGlobsHist(paramName, "Mu0", toyMu);
    histMu1 = getGlobsHist(paramName, "Mu1", toyMu);
    histMuFree = getGlobsHist(paramName, "MuFree", toyMu);
  }
  else if (paramType.EqualTo("Nuis")) {
    histMu0 = getNuisHist(paramName, "Mu0", toyMu);
    histMu1 = getNuisHist(paramName, "Mu1", toyMu);
    histMuFree = getNuisHist(paramName, "MuFree", toyMu);
  }
  
  TCanvas *can = new TCanvas("can", "can",800, 800);
  can->cd();
  
  // Format histograms:
  double min = 0.0001;
  double max = 10;
  histMu0->GetYaxis()->SetRangeUser(min,max);
  histMu0->Scale(1.0 / histMu0->Integral());
  histMu1->Scale(1.0 / histMu1->Integral());
  histMuFree->Scale(1.0 / histMuFree->Integral());
  histMu0->SetLineColor(kRed);
  histMu1->SetLineColor(kBlue);
  histMuFree->SetLineColor(kGreen+2);
  histMu0->SetLineWidth(3);
  histMu1->SetLineWidth(3);
  histMuFree->SetLineWidth(3);
  histMu0->SetLineStyle(1);
  histMu1->SetLineStyle(2);
  histMuFree->SetLineStyle(4);
  
  // Format axis titles:
  histMu0->GetYaxis()->SetTitle("Fraction of toys");
  if (paramName.Contains("nBkg") && !paramType.EqualTo("Globs")) {
    histMu0->GetXaxis()->SetTitle("Background normalization");
  }
  else {
    histMu0->GetXaxis()->SetTitle("parameter pulls [#sigma]");
  }
  
  // Draw histograms:
  gPad->SetLogy();
  histMu0->Draw("hist");
  histMu1->Draw("histSAME");
  histMuFree->Draw("histSAME");
  
  TLegend leg(0.2,0.71,0.4,0.86);
  leg.SetBorderSize(0);
  leg.SetTextSize(0.04);
  leg.SetFillColor(0);
  leg.AddEntry(histMu0,"#mu=0 fixed","l");
  leg.AddEntry(histMu1,"#mu=1 fixed","l");
  leg.AddEntry(histMuFree,"#mu floating","l");
  leg.Draw("SAME");
  
  TLatex typeText;
  typeText.SetNDC();
  typeText.SetTextColor(kBlack);
  typeText.SetTextFont(42); 
  typeText.SetTextSize(0.05);
  typeText.DrawLatex(0.2, 0.87, paramName);
  
  if (!paramName.Contains("nBkg")) {
    TLine *line[5];
    for (int i_l = 0; i_l < 5; i_l++) {
      line[i_l] = new TLine();
      line[i_l]->SetLineStyle(2); 
      line[i_l]->SetLineWidth(1);
      line[i_l]->SetLineColor(kBlack);
      line[i_l]->DrawLine(i_l-2, min, i_l-2, max);
    }
  }
  
  can->Print(Form("%s/plot_%s_toy%i.eps", outputDir.Data(), paramName.Data(), 
		  toyMu));
  can->Clear();
  gPad->SetLogy(0);
}

/**
   -----------------------------------------------------------------------------
   Plot the profiled values of the signal strength.  
*/
void DMToyAnalysis::plotProfiledMu() {
 
  TCanvas *can = new TCanvas("can", "can",800, 800);
  can->cd();
    
  TH1F *histMu0 = getMuHist(0);
  TH1F *histMu1 = getMuHist(1);
  histMu0->SetLineColor(kBlue);
  histMu1->SetLineColor(kRed);
  histMu0->SetLineWidth(2);
  histMu1->SetLineWidth(2);

  gPad->SetLogy();
  histMu0->Draw("");
  histMu1->Draw("SAME");
  
  TLegend leg(0.56,0.74,0.88,0.86);
  leg.SetBorderSize(0);
  leg.SetTextSize(0.04);
  leg.SetFillColor(0);
  leg.AddEntry(histMu0,"#mu=0 toy MC","l");
  leg.AddEntry(histMu1,"#mu=1 toy MC","l");
  leg.Draw("SAME");
  
  // Also get the mean signal strengths from the toy profiling:
  double meanMu0 = histMu0->GetMean();
  TLine *lineMu0 = new TLine();
  lineMu0->SetLineStyle(2);
  lineMu0->SetLineWidth(3);
  lineMu0->SetLineColor(kBlue); 
  lineMu0->DrawLine(meanMu0, histMu0->GetYaxis()->GetXmin(),
		    meanMu0, histMu0->GetYaxis()->GetXmax());
  double meanMu1 = histMu1->GetMean();
  TLine *lineMu1 = new TLine();
  lineMu1->SetLineStyle(2);
  lineMu1->SetLineWidth(3);
  lineMu1->SetLineColor(kRed);
  lineMu1->DrawLine(meanMu1, histMu0->GetYaxis()->GetXmin(),
		    meanMu1, histMu0->GetYaxis()->GetXmax());
  TLatex textMu0;
  textMu0.SetNDC();
  textMu0.SetTextColor(kBlue);
  textMu0.SetTextFont(42);
  textMu0.SetTextSize(0.04); 
  textMu0.DrawLatex(0.56, 0.65,Form("#mu=0 toy mean = %2.2f",meanMu0));
  TLatex textMu1;
  textMu1.SetNDC();
  textMu1.SetTextColor(kRed); 
  textMu1.SetTextFont(42);
  textMu1.SetTextSize(0.04);
  textMu1.DrawLatex(0.56, 0.7,Form("#mu=1 toy mean = %2.2f",meanMu1));
  
  can->Print(Form("%s/plot_profiledMu.eps", outputDir.Data()));
  can->Clear();
  gPad->SetLogy(0);  
}

/**
   -----------------------------------------------------------------------------
   Plot the distributions of nuisance parameters and global observables
   @param statistic -the name of the test statistic to plot.
*/
void DMToyAnalysis::plotTestStat(TString statistic) {
  TH1F *hStatMu0 = getStatHist(statistic, 0);
  TH1F *hStatMu1 = getStatHist(statistic, 1);

  TCanvas *can = new TCanvas("can", "can",800, 800);
  can->cd();
  
  hStatMu0->GetXaxis()->SetTitle(printStatName(statistic));
  hStatMu1->GetXaxis()->SetTitle(printStatName(statistic));
  
  hStatMu0->GetYaxis()->SetTitle("Normalized entries");  
  hStatMu1->GetYaxis()->SetTitle("Normalized entries");
  
  hStatMu0->SetLineColor(kBlue);
  hStatMu1->SetLineColor(kRed);
  
  // Draw test statistic histograms:
  gPad->SetLogy();
  hStatMu0->GetYaxis()->SetRangeUser(0.0001,1.0);
  hStatMu0->Draw("");
  hStatMu1->Draw("SAME");
  
  hAsymptotic->SetLineColor(kBlack);
  hAsymptotic->Draw("SAME");
  
  TLegend leg(0.49, 0.76, 0.84, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(hStatMu1, "#mu=1 toy MC","l");
  leg.AddEntry(hStatMu0, "#mu=0 toy MC","l");
  leg.AddEntry(hAsymptotic, "Asyptotic distribution","l");
  leg.Draw("SAME");
  
  can->Print(Form("%s/plot_%s.eps", outputDir.Data(), statistic.Data()));
  can->Clear();
  gPad->SetLogy(0);
}

/**
   -----------------------------------------------------------------------------
   Draw the asymptotic function along with toy test statistic, and comparisons.
   @param statistic - the name of the test statistic.
*/
void DMToyAnalysis::plotTestStatComparison(TString statistic) {
  
  // Construct canvas and sub-pads:
  TCanvas *can = new TCanvas("can", "can", 1000, 1200);
  can->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.6, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.4, 1.0, 0.6);
  TPad *pad3 = new TPad("pad2", "pad2", 0.0, 0.2, 1.0, 0.4);
  TPad *pad4 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.2);
  pad1->SetBottomMargin(0.00001);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.00001);
  pad3->SetTopMargin(0.00001);
  pad3->SetBottomMargin(0.00001);
  pad4->SetTopMargin(0.00001);
  pad4->SetBottomMargin(0.4);
  pad1->SetBorderMode(0);
  pad2->SetBorderMode(0);
  pad3->SetBorderMode(0);
  pad4->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();
  
  // Get the toy and asymptotic distributions:
  TH1F *hStatMu1 = getStatHist(statistic, 1);
  //hAsymptotic;
  
  // Calculate the histograms to plot:
  TH1F *hIntegralToy = new TH1F("hIntegralToy", "hIntegralToy", 
				nBins, binMin, binMax);
  TH1F *hIntegralAsym = new TH1F("hIntegralAsym", "hIntegralAsym",
				 nBins, binMin, binMax);
  TH1F *hRatio = new TH1F("hRatio", "hRatio", nBins, binMin, binMax);
  TH1F *hSignificanceToy = new TH1F("hSignificanceToy", "hSignificanceToy",
				    nBins, binMin, binMax);
  TH1F *hSignificanceAsym = new TH1F("hSignificanceAsym", "hSignificanceAsym",
				     nBins, binMin, binMax);
  for (int i_b = 1; i_b <= nBins; i_b++) {
    double valueAsym = hAsymptotic->Integral(i_b, nBins);
    double valueToy = hStatMu1->Integral(i_b, nBins);
    hIntegralToy->SetBinContent(i_b, valueToy);
    hIntegralAsym->SetBinContent(i_b, valueAsym);
    hSignificanceToy->SetBinContent(i_b, -1.0*TMath::NormQuantile(valueToy));
    hSignificanceAsym->SetBinContent(i_b, -1.0*TMath::NormQuantile(valueAsym));
    
    double ratio = (TMath::NormQuantile(valueToy) == 0) ? 1.0 :
      TMath::NormQuantile(valueAsym) / TMath::NormQuantile(valueToy);
    
    hRatio->SetBinContent(i_b, ratio);
  }
   
  // Pad 1: Draw test statistic under mu=1 hypothesis and asymptotic function.
  pad1->cd();
  hAsymptotic->SetLineColor(kBlue);
  hStatMu1->Draw("");
  hAsymptotic->Draw("SAME");
  gPad->SetLogy();
  gPad->SetLogx();
  
  TLegend leg(0.65, 0.65, 0.9, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.06);
  TString printName = printStatName(statistic);
  leg.AddEntry(hStatMu1,Form("Toy MC %s",printName.Data()),"l");
  leg.AddEntry(hAsymptotic,Form("Asymptotic %s",printName.Data()),"l");
  leg.Draw("SAME");
  TLatex textToy;
  textToy.SetNDC();
  textToy.SetTextColor(kRed); 
  textToy.SetTextFont(42);
  textToy.SetTextSize(0.05);
  textToy.DrawLatex(0.2, 0.2, Form("Toy Mean = %2.2f", hStatMu1->GetMean()));
  TLatex textAsym;
  textAsym.SetNDC();
  textAsym.SetTextColor(kBlue); 
  textAsym.SetTextFont(42);
  textAsym.SetTextSize(0.05);
  textAsym.DrawLatex(0.2, 0.15, Form("Asymptotics Mean = %2.2f", 
				     hAsymptotic->GetMean()));
  
  // Pad 2: Draw the CDFs for asymptotics and pseudo-experiments
  pad2->cd();
  hIntegralToy->SetLineColor(kRed);
  hIntegralAsym->SetLineColor(kBlue);
  hIntegralToy->GetXaxis()->SetTitle(printName);
  hIntegralAsym->GetXaxis()->SetTitle(printName);
  hIntegralToy->GetYaxis()->SetTitle("CDF");
  hIntegralAsym->GetYaxis()->SetTitle("CDF");
  hIntegralToy->GetYaxis()->SetTitleSize(0.1);
  hIntegralToy->GetYaxis()->SetLabelSize(0.1);
  hIntegralToy->GetYaxis()->SetTitleOffset(0.65);
  hIntegralToy->Draw("");
  hIntegralAsym->Draw("SAME");
  gPad->SetLogy();
  gPad->SetLogx();
  
  // Pad 3: Plot the calculated significance corresponding to the CDF on Pad 2.
  pad3->cd();
  hSignificanceToy->SetLineWidth(2);
  hSignificanceAsym->SetLineWidth(2);
  hSignificanceToy->SetLineColor(kRed);
  hSignificanceAsym->SetLineColor(kBlue);
  hSignificanceToy->GetXaxis()->SetTitle(printName);
  hSignificanceAsym->GetXaxis()->SetTitle(printName);
  hSignificanceToy->GetYaxis()->SetTitle("Z [#sigma]");
  hSignificanceAsym->GetYaxis()->SetTitle("Z [#sigma]");
  hSignificanceToy->SetBinContent(1, 0);
  hSignificanceAsym->SetBinContent(1, 0);
  hSignificanceToy->GetYaxis()->SetTitleSize(0.1);
  hSignificanceToy->GetYaxis()->SetLabelSize(0.1);
  hSignificanceToy->GetYaxis()->SetTitleOffset(0.65);
  hSignificanceToy->GetYaxis()->SetNdivisions(6);
  hSignificanceToy->Draw("");
  hSignificanceAsym->Draw("SAME");
  gPad->SetLogx();
  
  // Pad 4: Calculate the ratio of Z values from Pad 3.
  pad4->cd();
  hRatio->SetLineColor(kBlack);
  hRatio->SetLineWidth(2);
  hRatio->GetYaxis()->SetTitle("Z_{Asym.} / Z_{Toy}");
  hRatio->GetXaxis()->SetTitle(printName);
  hRatio->GetYaxis()->SetTitleSize(0.1);
  hRatio->GetYaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetTitleOffset(0.65);
  hRatio->GetYaxis()->SetNdivisions(6);
  hRatio->GetXaxis()->SetTitleSize(0.1);
  hRatio->GetXaxis()->SetLabelSize(0.1);
  hRatio->GetXaxis()->SetTitleOffset(0.9);
  hRatio->Draw();
  gPad->SetLogx();
  TLine *line2 = new TLine();
  line2->SetLineStyle(2);
  line2->SetLineWidth(1);
  line2->SetLineColor(kBlack);
  line2->DrawLine(binMin, 1.0, binMax, 1.0);
  
  can->Print(Form("%s/plot_comp_%s.eps", outputDir.Data(), statistic.Data()));
  can->Print(Form("%s/plot_comp_%s.png", outputDir.Data(), statistic.Data()));
  can->Clear();
}

/**
   -----------------------------------------------------------------------------
   Convert the statistic name into a LaTex formatted name.
   @param statistic - the name of the test statistic.
*/
TString DMToyAnalysis::printStatName(TString statistic) {
  if (statistic.EqualTo("Q0")) return TString("q_{0}");
  else if (statistic.EqualTo("QMu")) return TString("q_{#mu}");
  else if (statistic.EqualTo("QMuTilde")) return TString("#tilde{q}_{#mu}");
  else return TString("q");
}

