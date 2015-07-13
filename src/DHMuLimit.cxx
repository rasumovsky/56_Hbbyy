/*
Author: Aaron Armbruster
Date:   2012-05-25
Email:  armbrusa@umich.edu
Description: Script to run asymptotic CLs.

--------
00-01-00
-First version with updated bands

--------
00-01-01
-Fixed problem in asimov data creation that affected +1,2sigma bands

--------
00-01-02
-(Re)added support for non-sim pdfs (still need to be extended)
-Fixed default doFit arg of makeAsimovData
-Added better output for unresolved fit failures
-Improved retry loop for fit failures


/////////////////////
//////PREAMBLE///////
/////////////////////

The script uses an iterative process to find the crossing of qmu with the qmu95(mu/sigma) curve, 
where qmu95(mu/sigma) is found assuming asymptotic formula for the distribution of the 
test statistic f(qmu|mu') (arxiv 1007.1727) and of the test statistic qmu (or tilde)

The sequence is

mu_i+1 = mu_i - gamma_i*(mu_i - mu'_i)

where gamma_i is a dynamic damping factor used for convergence (nominal gamma_i = 1), and mu'_i is 
determined by extrapolating the test statistic to the qmu95 curve assuming qmu is parabolic:

qmu'_i = (mu'_i - muhat)^2 / sigma_i^2 = qmu95(mu'_i / sigma_i)

where sigma_i is determined by computing qmu_i (not '):

sigma_i = (mu_i - muhat) / sqrt(qmu_i)

At the crossing qmu_N = qmu95 the assumption that qmu is a parabola goes away, 
so we're not ultimately dependent on this assumption beyond its use in the asymptotic formula.

The sequence ends when the relative correction factor gamma*(mu_i - mu'_i) / mu_i is less than some
specified precision (0.005 by default)




///////////////////////////
//////AFTER RUNNING////////
///////////////////////////


The results will be printed as well as stored in a root file in the folder 'root/<folder>', where <folder>
is specified by you (default 'test')

The root file has a 7-bin TH1D named 'limit', where each bin is filled with the upper limit values in this order:

1: Observed
2: Median
3: +2 sigma
4: +1 sigma
5: -1 sigma
6: -2 sigma
7: mu=0 fit status (only meaningful if asimov data is generated within the macro)

It will also store the result of the old bands procedure in a TH1D named 'limit_old'. 




//////////////////////////


This version is functionally fully consistent with the previous tag.

NOTE: The script runs significantly faster when compiled
*/




#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "statistics.h"
#include "DHAnalysis.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

TString inputDHSignal;

//band configuration
bool betterBands           = 1; // (recommendation = 1) improve bands by using a more appropriate asimov dataset for those points
bool betterNegativeBands   = 0; // (recommendation = 0) improve also the negative bands
bool profileNegativeAtZero = 0; // (recommendation = 0) profile asimov for negative bands at zero


//other configuration
string defaultMinimizer    = "Minuit";     // or "Minuit"
int defaultPrintLevel      = -1;            // Minuit print level
int defaultStrategy        = 1;             // Minimization strategy. 0-2. 0 = fastest, least robust. 2 = slowest, most robust
bool killBelowFatal        = 1;             // In case you want to suppress RooFit warnings further, set to 1
bool doBlind               = 0;             // in case your analysis is blinded
bool conditionalExpected   = 1 && !DHAnalysis::doBlind; // Profiling mode for Asimov data: 0 = conditional MLEs, 1 = nominal MLEs
bool doTilde               = 0;             // bound mu at zero if true and do the \tilde{q}_{mu} asymptotics
bool doExp                 = 1;             // compute expected limit
bool doObs                 = 1 && !DHAnalysis::doBlind; // compute observed limit
double precision           = 0.005;         // % precision in mu that defines iterative cutoff
bool verbose               = 1;             // 1 = very spammy
bool usePredictiveFit      = 0;             // experimental, extrapolate best fit nuisance parameters based on previous fit results
bool extrapolateSigma      = 0;             // experimantal, extrapolate sigma based on previous fits
int maxRetries             = 3;             // number of minimize(fcn) retries before giving up





//don't touch!
map<RooNLLVar*, double> map_nll_muhat;
map<RooNLLVar*, double> map_muhat;
map<RooDataSet*, RooNLLVar*> map_data_nll;
map<RooNLLVar*, string> map_snapshots;
map<RooNLLVar*, map<double, double> > map_nll_mu_sigma;
RooWorkspace* w = NULL;
ModelConfig* mc = NULL;
RooDataSet* data = NULL;
RooRealVar* firstPOI = NULL;
RooNLLVar* asimov_0_nll = NULL;
RooNLLVar* obs_nll = NULL;
int nrMinimize=0;
int direction=1;
int global_status=0;
double target_CLs=0.05;

//main
void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName,
		       const char* modelConfigName,
		       const char* dataName,
		       const char* asimovDataName,
		       string folder,
		       double CL,
		       TString option);

 //for backwards compatibility
void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName = "combined",
		       const char* modelConfigName = "ModelConfig",
		       const char* dataName = "obsData",
		       const char* asimovDataName = "AsimovDataMu0",
		       const char* conditionalSnapshot = "conditionalGlobs_0",
		       const char* nominalSnapshot = "nominalGlobs",
		       string folder = "test",
		       double CL = 0.95,
		       TString option="");

double getLimit(RooNLLVar* nll, double initial_guess = 0);
double getSigma(RooNLLVar* nll, double mu, double muhat, double& qmu);
double getQmu(RooNLLVar* nll, double mu);
void saveSnapshot(RooNLLVar* nll, double mu);
void loadSnapshot(RooNLLVar* nll, double mu);
void doPredictiveFit(RooNLLVar* nll, double mu1, double m2, double mu);
RooNLLVar* createNLL(RooDataSet* _data);
double getNLL(RooNLLVar* nll);
double findCrossing(double sigma_obs, double sigma, double muhat);
void setMu(double mu);
double getQmu95_brute(double sigma, double mu);
double getQmu95(double sigma, double mu);
double calcCLs(double qmu_tilde, double sigma, double mu);
double calcPmu(double qmu_tilde, double sigma, double mu);
double calcPb(double qmu_tilde, double sigma, double mu);
double calcDerCLs(double qmu, double sigma, double mu);
int minimize(RooNLLVar* nll);
int minimize(RooAbsReal* nll);

void unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter);
RooDataSet* makeAsimovData(bool doConditional, RooNLLVar* conditioning_nll, double mu_val, string* mu_str = NULL, string* mu_prof_str = NULL, double mu_val_profile = -999, bool doFit = true);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// runAsymptoticsCLs:

void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName,
		       const char* modelConfigName,
		       const char* dataName,
		       const char* asimovDataName,
		       const char* conditionalSnapshot,
		       const char* nominalSnapshot,
		       string folder,
		       double CL,
		       TString option)
{
  conditionalSnapshot = ""; // warningless compile
  nominalSnapshot = "";     // warningless compile

  runAsymptoticsCLs(infile, workspaceName, modelConfigName, dataName, asimovDataName, folder, CL, option);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// runAsymptoticsCLs:

void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName,
		       const char* modelConfigName,
		       const char* dataName,
		       const char* asimovDataName,
		       string folder,
		       double CL,
		       TString option)
{
  cout << "runAsymptoticsCLs( " << infile << ", " << workspaceName << ", " << modelConfigName << ", " << dataName << ", " << asimovDataName << ", " << folder << ", " << CL << ", " << option << " )" << endl;
  TStopwatch timer;
  timer.Start();
  
  if (killBelowFatal) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(defaultMinimizer.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(defaultStrategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(defaultPrintLevel);
  
  //check inputs
  TFile f(infile);
  w = (RooWorkspace*)f.Get(workspaceName);
  if (!w)
  {
    cout << "ERROR::Workspace: " << workspaceName << " doesn't exist!" << endl;
    return;
  }

  mc = (ModelConfig*)w->obj(modelConfigName);
  if (!mc)
  {
    cout << "ERROR::ModelConfig: " << modelConfigName << " doesn't exist!" << endl;
    return;
  }
  firstPOI = (RooRealVar*)mc->GetParametersOfInterest()->first();

  data = (RooDataSet*)w->data(dataName);
  if (!data)
  {
    cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
    return;
  }
  
  if(option.Contains("nosys")){
    if(w->set("nuisAll")) statistics::constSet((RooArgSet*)w->set("nuisAll"),true);
  }
  obs_nll = createNLL(data);
  map_snapshots[obs_nll] = "nominalGlobs";
  map_data_nll[data] = obs_nll;
  w->saveSnapshot("nominalGlobs",*mc->GetGlobalObservables());
  w->saveSnapshot("nominalNuis",*mc->GetNuisanceParameters());
  
  global_status=0;
  RooDataSet* asimovData_0 = (RooDataSet*)w->data(asimovDataName);
  if( !asimovData_0 )
  {
    asimovData_0 = makeAsimovData(conditionalExpected, obs_nll, 0);
  }
  int asimov0_status=global_status;
  asimov_0_nll = createNLL(asimovData_0);
  map_snapshots[asimov_0_nll] = "conditionalGlobs_0";
  map_data_nll[asimovData_0] = asimov_0_nll;
  setMu(0);
  map_muhat[asimov_0_nll] = 0;
  saveSnapshot(asimov_0_nll, 0);
  w->loadSnapshot("conditionalNuis_0");
  w->loadSnapshot("conditionalGlobs_0");
  map_nll_muhat[asimov_0_nll] = asimov_0_nll->getVal();

  
  target_CLs=1-CL;

  double med_limit = doExp ? getLimit(asimov_0_nll, 1.0) : 1.0;
  int med_status=global_status;

  double sigma = med_limit/sqrt(3.84); // pretty close
  double mu_up_p4_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 4), 1) + 4);
  double mu_up_p3_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 3), 1) + 3);
  double mu_up_p2_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 2), 1) + 2);
  double mu_up_p1_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 1), 1) + 1);
  double mu_up_n1_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-1), 1) - 1);
  double mu_up_n2_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-2), 1) - 2);
  double mu_up_n3_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-3), 1) - 3);
  double mu_up_n4_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-4), 1) - 4);

  double mu_up_p4 = mu_up_p4_approx;
  double mu_up_p3 = mu_up_p3_approx;
  double mu_up_p2 = mu_up_p2_approx;
  double mu_up_p1 = mu_up_p1_approx;
  double mu_up_n1 = mu_up_n1_approx;
  double mu_up_n2 = mu_up_n2_approx;
  double mu_up_n3 = mu_up_n3_approx;
  double mu_up_n4 = mu_up_n4_approx;

  firstPOI->setRange(-5*sigma, 5*sigma);
  map<int, int> N_status;
  if (betterBands && doExp) // no better time than now to do this
  {
    //find quantiles, starting with +2, since we should be at +1.96 right now

    double init_targetCLs = target_CLs;
    firstPOI->setRange(-5*sigma, 5*sigma);
    for (int N=2;N>=-2;N--)
    {
      if (N < 0 && !betterNegativeBands) continue;
      if (N == 0) continue;
      target_CLs=2*(1-ROOT::Math::gaussian_cdf(fabs(N))); // change this so findCrossing looks for sqrt(qmu95)=2
      if (N < 0) direction = -1;

      //get the acual value
      double NtimesSigma = getLimit(asimov_0_nll, N*med_limit/sqrt(3.84)); // use N * sigma(0) as an initial guess
      N_status[N] += global_status;
      sigma = NtimesSigma/N;
      cout << endl;
      cout << "Found N * sigma = " << N << " * " << sigma << endl;

      string muStr,muStrPr;
      w->loadSnapshot("conditionalGlobs_0");
      double pr_val = NtimesSigma;
      if (N < 0 && profileNegativeAtZero) pr_val = 0;
      RooDataSet* asimovData_N = makeAsimovData(1, asimov_0_nll, NtimesSigma, &muStr, &muStrPr, pr_val, 0);
      

      RooNLLVar* asimov_N_nll = createNLL(asimovData_N);
      map_data_nll[asimovData_N] = asimov_N_nll;
      map_snapshots[asimov_N_nll] = "conditionalGlobs"+muStrPr;
      w->loadSnapshot(map_snapshots[asimov_N_nll].c_str());
      w->loadSnapshot(("conditionalNuis"+muStrPr).c_str());
      setMu(NtimesSigma);

      double nll_val = asimov_N_nll->getVal();
      saveSnapshot(asimov_N_nll, NtimesSigma);
      map_muhat[asimov_N_nll] = NtimesSigma;
      if (N < 0 && doTilde)
      {
	setMu(0);
	firstPOI->setConstant(1);
	nll_val = getNLL(asimov_N_nll);
      }
      map_nll_muhat[asimov_N_nll] = nll_val;

      target_CLs = init_targetCLs;
      direction=1;
      double initial_guess = findCrossing(NtimesSigma/N, NtimesSigma/N, NtimesSigma);
      double limit = getLimit(asimov_N_nll, initial_guess);
      N_status[N] += global_status;

      if (N == 2) mu_up_p2 = limit;
      else if (N == 1) mu_up_p1 = limit;
      else if (N ==-1) mu_up_n1 = limit;
      else if (N ==-2) mu_up_n2 = limit;
      //return;
    }
    direction = 1;
    target_CLs = init_targetCLs;

  }

  w->loadSnapshot("conditionalNuis_0");
  double obs_limit = doObs ? getLimit(obs_nll, med_limit) : 0;
  int obs_status=global_status;

  bool hasFailures = false;
  if (obs_status != 0 || med_status != 0 || asimov0_status != 0) hasFailures = true;
  for (map<int, int>::iterator itr=N_status.begin();itr!=N_status.end();itr++)
  {
    if (itr->second != 0) hasFailures = true;
  }
  if (hasFailures)
  {
    cout << "--------------------------------" << endl;
    cout << "Unresolved fit failures detected" << endl;
    cout << "Asimov0:  " << asimov0_status << endl;
    for (map<int, int>::iterator itr=N_status.begin();itr!=N_status.end();itr++)
    {
      cout << "+" << itr->first << "sigma:  " << itr->first << endl;
    }
    cout << "Median:   " << med_status << endl;
    cout << "Observed: " << obs_status << endl;
    cout << "--------------------------------" << endl;
  }

  if (betterBands) cout << "Guess for bands" << endl;
  cout << "+2sigma:  " << mu_up_p2_approx << endl;
  cout << "+1sigma:  " << mu_up_p1_approx << endl;
  cout << "-1sigma:  " << mu_up_n1_approx << endl;
  cout << "-2sigma:  " << mu_up_n2_approx << endl;
  if (betterBands)
  {
    cout << endl;
    cout << "Correct bands" << endl;
    cout << "+2sigma:  " << mu_up_p2 << endl;
    cout << "+1sigma:  " << mu_up_p1 << endl;
    cout << "-1sigma:  " << mu_up_n1 << endl;
    cout << "-2sigma:  " << mu_up_n2 << endl;

  }

  cout << "Median:   " << med_limit << endl;
  cout << "Observed: " << obs_limit << endl;
  cout << endl;
  
  TFile fout(Form("%s/file_CLs_%s.root",folder.c_str(),inputDHSignal.Data()),"recreate");
  
  ofstream ftxt(Form("%s/text_CLs_%s.txt",folder.c_str(),inputDHSignal.Data()));
  ftxt << "CLs"
       << "\t" << inputDHSignal
       << "\t" << obs_limit
       << "\t" << med_limit
       << "\t" << mu_up_p2
       << "\t" << mu_up_p1
       << "\t" << mu_up_n1
       << "\t" << mu_up_n2
       << endl;
  ftxt.close();

  TH1D* h_lim = new TH1D("limit","limit",7,0,7);
  h_lim->SetBinContent(1, obs_limit);
  h_lim->SetBinContent(2, med_limit);
  h_lim->SetBinContent(3, mu_up_p2);
  h_lim->SetBinContent(4, mu_up_p1);
  h_lim->SetBinContent(5, mu_up_n1);
  h_lim->SetBinContent(6, mu_up_n2);
  h_lim->SetBinContent(7, global_status);

  h_lim->GetXaxis()->SetBinLabel(1, "Observed");
  h_lim->GetXaxis()->SetBinLabel(2, "Expected");
  h_lim->GetXaxis()->SetBinLabel(3, "+2sigma");
  h_lim->GetXaxis()->SetBinLabel(4, "+1sigma");
  h_lim->GetXaxis()->SetBinLabel(5, "-1sigma");
  h_lim->GetXaxis()->SetBinLabel(6, "-2sigma");
  h_lim->GetXaxis()->SetBinLabel(7, "Global status"); // do something with this later

  TH1D* h_lim_old = new TH1D("limit_old","limit_old",7,0,7); // include also old approximation of bands
  h_lim_old->SetBinContent(1, obs_limit);
  h_lim_old->SetBinContent(2, med_limit);
  h_lim_old->SetBinContent(3, mu_up_p2_approx);
  h_lim_old->SetBinContent(4, mu_up_p1_approx);
  h_lim_old->SetBinContent(5, mu_up_n1_approx);
  h_lim_old->SetBinContent(6, mu_up_n2_approx);
  h_lim_old->SetBinContent(7, global_status);

  h_lim_old->GetXaxis()->SetBinLabel(1, "Observed");
  h_lim_old->GetXaxis()->SetBinLabel(2, "Expected");
  h_lim_old->GetXaxis()->SetBinLabel(3, "+2sigma");
  h_lim_old->GetXaxis()->SetBinLabel(4, "+1sigma");
  h_lim_old->GetXaxis()->SetBinLabel(5, "-1sigma");
  h_lim_old->GetXaxis()->SetBinLabel(6, "-2sigma");
  h_lim_old->GetXaxis()->SetBinLabel(7, "Global status"); 

  fout.Write();
  fout.Close();

  cout << "Finished with " << nrMinimize << " calls to minimize(nll)" << endl;
  timer.Print();
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// getLimit:

double getLimit( RooNLLVar* nll, double initial_guess )
{
  cout << "------------------------" << endl;
  cout << "getLimit( " << nll->GetName() << ", " << initial_guess << " )" << endl;
  
  //get initial guess based on muhat and sigma(muhat)
  firstPOI->setConstant(0);
  global_status=0;

  if (nll == asimov_0_nll) 
  {
    setMu(0);
    firstPOI->setConstant(1);
  }

  double muhat;
  if (map_nll_muhat.find(nll) == map_nll_muhat.end())
  {
    double nll_val = getNLL(nll);
    muhat = firstPOI->getVal();
    saveSnapshot(nll, muhat);
    map_muhat[nll] = muhat;
    if (muhat < 0 && doTilde)
    {
      setMu(0);
      firstPOI->setConstant(1);
      nll_val = getNLL(nll);
    }

    map_nll_muhat[nll] = nll_val;
  }
  else
  {
    muhat = map_muhat[nll];
  }

  if (muhat < 0.1 || initial_guess != 0) setMu(initial_guess);
  double qmu,qmuA;
  double sigma_guess = getSigma( asimov_0_nll, firstPOI->getVal(), 0, qmu );
  double sigma_b = sigma_guess;
  double mu_guess = findCrossing( sigma_guess, sigma_b, muhat );
  double pmu = calcPmu(qmu, sigma_b, mu_guess);
  double pb = calcPb(qmu, sigma_b, mu_guess);
  double CLs = calcCLs(qmu, sigma_b, mu_guess);
  double qmu95 = getQmu95(sigma_b, mu_guess);
  setMu(mu_guess);

  cout << "Initial guess:  " << mu_guess << endl;
  cout << "Sigma(obs):     " << sigma_guess << endl;
  cout << "Sigma(mu,0):    " << sigma_b << endl;
  cout << "muhat:          " << muhat << endl;
  cout << "qmu95:          " << qmu95 << endl;
  cout << "qmu:            " << qmu << endl;
  cout << "pmu:            " << pmu << endl;
  cout << "1-pb:           " << pb << endl;
  cout << "CLs:            " << CLs << endl;
  cout << endl;

  int nrDamping = 1;
  map<double, double> guess_to_corr;
  double damping_factor = 1.0;
  //double damping_factor_pre = damping_factor;
  int nrItr = 0;
  double mu_pre = muhat;
  double mu_pre2 = muhat;
  while (fabs(mu_pre-mu_guess) > precision*mu_guess*direction)
  {
    cout << "----------------------" << endl;
    cout << "Starting iteration " << nrItr << " of " << nll->GetName() << endl;
    // do this to avoid comparing multiple minima in the conditional and unconditional fits
    if (nrItr == 0) loadSnapshot(nll, muhat);
    else if (usePredictiveFit) doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
    else loadSnapshot(asimov_0_nll, mu_pre);

    sigma_guess=getSigma(nll, mu_guess, muhat, qmu);
    saveSnapshot(nll, mu_guess);


    if (nll != asimov_0_nll)
    {
      if (nrItr == 0) loadSnapshot(asimov_0_nll, map_nll_muhat[asimov_0_nll]);
      else if (usePredictiveFit) 
      {
	if (nrItr == 1) doPredictiveFit(nll, map_nll_muhat[asimov_0_nll], mu_pre, mu_guess);
	else doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
      }
      else loadSnapshot(asimov_0_nll, mu_pre);

      sigma_b=getSigma(asimov_0_nll, mu_guess, 0, qmuA);
      saveSnapshot(asimov_0_nll, mu_guess);
    }
    else
    {
      sigma_b=sigma_guess;
      qmuA=qmu;
    }

    double corr = damping_factor*(mu_guess - findCrossing(sigma_guess, sigma_b, muhat));
    for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++)
    {
      if (fabs(itr->first - (mu_guess-corr)) < direction*mu_guess*0.02 && fabs(corr) > direction*mu_guess*precision) 
      {
	damping_factor *= 0.8;
	cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
	if (nrDamping++ > 10)
	{
	  nrDamping = 1;
	  damping_factor = 1.0;
	}
	corr *= damping_factor;
	break;
      }
    }

    //subtract off the difference in the new and damped correction
    guess_to_corr[mu_guess] = corr;
    mu_pre2 = mu_pre;
    mu_pre = mu_guess;
    mu_guess -= corr;


    pmu = calcPmu(qmu, sigma_b, mu_pre);
    pb = calcPb(qmu, sigma_b, mu_pre);
    CLs = calcCLs(qmu, sigma_b, mu_pre);
    qmu95 = getQmu95(sigma_b, mu_pre);


    cout << "NLL:            " << nll->GetName() << endl;
    cout << "Previous guess: " << mu_pre << endl;
    cout << "Sigma(obs):     " << sigma_guess << endl;
    cout << "Sigma(mu,0):    " << sigma_b << endl;
    cout << "muhat:          " << muhat << endl;
    cout << "pmu:            " << pmu << endl;
    cout << "1-pb:           " << pb << endl;
    cout << "CLs:            " << CLs << endl;
    cout << "qmu95:          " << qmu95 << endl;
    cout << "qmu:            " << qmu << endl;
    cout << "qmuA0:          " << qmuA << endl;
    cout << "Precision:      " << direction*mu_guess*precision << endl;
    cout << "Correction:    "  << (-corr<0?" ":"") << -corr << endl;
    cout << "New guess:      " << mu_guess << endl;
    cout << endl;

      

    nrItr++;
    if (nrItr > 25)
    {
      cout << "Infinite loop detected in getLimit(). Please intervene." << endl;
      break;
    }
  }
  
  cout << "Found limit for nll " << nll->GetName() << ": " << mu_guess << endl;
  cout << "Finished in " << nrItr << " iterations." << endl;
  cout << endl;
  return mu_guess;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// getSigma:

double getSigma( RooNLLVar* nll, double mu, double muhat, double& qmu )
{
  cout << "getSigma( " << nll->GetName() << ", " << mu << ", " << muhat << ", " << qmu << " )" << endl;
  qmu = getQmu( nll, mu );
  if( verbose ) cout << "qmu = " << qmu << endl;
  if( mu*direction < muhat ) return fabs(mu-muhat)/sqrt(qmu);
  else if (muhat < 0 && doTilde) return sqrt(mu*mu-2*mu*muhat*direction)/sqrt(qmu);
  else return (mu-muhat)*direction/sqrt(qmu);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// getQmu:

double getQmu(RooNLLVar* nll, double mu)
{
  cout << "getQmu( " << nll->GetName() << ", " << mu << " )" << endl;
  double nll_muhat = map_nll_muhat[nll];
  bool isConst = firstPOI->isConstant();
  firstPOI->setConstant(1);
  setMu(mu);
  double nll_val = getNLL(nll);
  firstPOI->setConstant(isConst);
  //cout << "  qmu = 2 * (" << nll_val << " - " << nll_muhat << ")" << endl;
  return 2*(nll_val-nll_muhat);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void saveSnapshot(RooNLLVar* nll, double mu)
{
  stringstream snapshotName;
  snapshotName << nll->GetName() << "_" << mu;
  w->saveSnapshot(snapshotName.str().c_str(), *mc->GetNuisanceParameters());
}

void loadSnapshot(RooNLLVar* nll, double mu)
{
  stringstream snapshotName;
  snapshotName << nll->GetName() << "_" << mu;
  w->loadSnapshot(snapshotName.str().c_str());
}

void doPredictiveFit(RooNLLVar* nll, double mu1, double mu2, double mu)
{
  if (fabs(mu2-mu) < direction*mu*precision*4)
  {
    loadSnapshot(nll, mu2);
    return;
  }

//extrapolate to mu using mu1 and mu2 assuming nuis scale linear in mu
  const RooArgSet* nuis = mc->GetNuisanceParameters();
  int nrNuis = nuis->getSize();
  double* theta_mu1 = new double[nrNuis];
  double* theta_mu2 = new double[nrNuis];

  TIterator* itr = nuis->createIterator();
  RooRealVar* var=NULL;
  int counter = 0;
  loadSnapshot(nll, mu1);
  while ((var == (RooRealVar*)itr->Next()))
  {
    theta_mu1[counter++] = var->getVal();
  }

  itr->Reset();
  counter = 0;
  loadSnapshot(nll, mu2);
  while ((var == (RooRealVar*)itr->Next()))
  {
    theta_mu2[counter++] = var->getVal();
  }

  itr->Reset();
  counter = 0;
  while ((var == (RooRealVar*)itr->Next()))
  {
    double m = (theta_mu2[counter] - theta_mu1[counter])/(mu2-mu1);
    double b = theta_mu2[counter] - m*mu2;
    double theta_extrap = m*mu+b;
    
    var->setVal(theta_extrap);
    counter++;
  }

  delete itr;
  delete theta_mu1;
  delete theta_mu2;
}

RooNLLVar* createNLL(RooDataSet* _data)
{
  RooArgSet nuis = *mc->GetNuisanceParameters();
  RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*_data, Constrain(nuis));
  return nll;
}

double getNLL(RooNLLVar* nll)
{
  string snapshotName = map_snapshots[nll];
  if (snapshotName != "") w->loadSnapshot(snapshotName.c_str());
  minimize(nll);
  double val = nll->getVal();
  w->loadSnapshot("nominalGlobs");
  return val;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// findCrossing:

double findCrossing(double sigma_obs, double sigma, double muhat)
{
  cout << "findCrossing( " << sigma_obs << ", " << sigma << ", " << muhat << " )" << endl;
  double mu_guess = muhat + ROOT::Math::gaussian_quantile(1-target_CLs,1)*sigma_obs*direction;
  int nrItr = 0;
  int nrDamping = 1;

  map<double, double> guess_to_corr;
  double damping_factor = 1.0;
  double mu_pre = mu_guess - 10*mu_guess*precision;
  while (fabs(mu_guess-mu_pre) > direction*mu_guess*precision)
  {
    double sigma_obs_extrap = sigma_obs;
    double eps = 0;
    if (extrapolateSigma)
    {
      //map<double, double> map_mu_sigma = map_nll_mu_sigma[nll];
    }

    mu_pre = mu_guess;

    double qmu95 = getQmu95(sigma, mu_guess);
    double qmu;
    qmu = 1./sigma_obs_extrap/sigma_obs_extrap*(mu_guess-muhat)*(mu_guess-muhat);
    if (muhat < 0 && doTilde) qmu = 1./sigma_obs_extrap/sigma_obs_extrap*(mu_guess*mu_guess-2*mu_guess*muhat);

    double dqmu_dmu = 2*(mu_guess-muhat)/sigma_obs_extrap/sigma_obs_extrap - 2*qmu*eps;

    double corr = damping_factor*(qmu-qmu95)/dqmu_dmu;
    for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++)
    {
      if (fabs(itr->first - mu_guess) < direction*mu_guess*precision) 
      {
	damping_factor *= 0.8;
	if (verbose) cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
	if (nrDamping++ > 10)
	{
	  nrDamping = 1;
	  damping_factor = 1.0;
	}
	corr *= damping_factor;
	break;
      }
    }
    guess_to_corr[mu_guess] = corr;

    mu_guess = mu_guess - corr;
    nrItr++;
    if (nrItr > 100)
    {
      cout << "Infinite loop detected in findCrossing. Please intervene." << endl;
      exit(1);
    }
    if (verbose) cout << "mu_guess = " << mu_guess << ", mu_pre = " << mu_pre << ", qmu = " << qmu << ", qmu95 = " << qmu95 << ", sigma_obs_extrap = " << sigma_obs_extrap << ", sigma = " << sigma << ", direction*mu*prec = " << direction*mu_guess*precision << endl;
  }

  return mu_guess;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// setMu:

void setMu(double mu)
{
  cout << "setMu( " << mu << " )" << endl;
  if (mu != mu)
  {
    cout << "  ERROR::POI gave nan. Please intervene." << endl;
    exit(1);
  }
  if( mu > 0 && firstPOI->getMax() < mu ) firstPOI->setMax(2*mu);
  if( mu < 0 && firstPOI->getMin() > mu ) firstPOI->setMin(2*mu);
  firstPOI->setVal(mu);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// getQmu95_brute:

double getQmu95_brute(double sigma, double mu)
{
  double step_size = 0.001;
  double start = step_size;
  if (mu/sigma > 0.2) start = 0;
  for (double qmu=start;qmu<20;qmu+=step_size)
  {
    double CLs = calcCLs(qmu, sigma, mu);

    if (CLs < target_CLs) return qmu;
  }

  return 20;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// getQmu95:

double getQmu95(double sigma, double mu)
{
  double qmu95 = 0;
  //no sane man would venture this far down into |mu/sigma|
  double target_N = ROOT::Math::gaussian_cdf(1-target_CLs,1);
  if (fabs(mu/sigma) < 0.25*target_N)
  {
    qmu95 = 5.83/target_N;
  }
  else
  {
    map<double, double> guess_to_corr;
    double qmu95_guess = pow(ROOT::Math::gaussian_quantile(1-target_CLs,1),2);
    int nrItr = 0;
    int nrDamping = 1;
    double damping_factor = 1.0;
    double qmu95_pre = qmu95_guess - 10*2*qmu95_guess*precision;
    while (fabs(qmu95_guess-qmu95_pre) > 2*qmu95_guess*precision)
    {
      qmu95_pre = qmu95_guess;
      if (verbose)
      {
	cout << "qmu95_guess = " << qmu95_guess << endl;
	cout << "CLs = " << calcCLs(qmu95_guess, sigma, mu) << endl;
	cout << "Derivative = " << calcDerCLs(qmu95_guess, sigma, mu) << endl;
      }

      double corr = damping_factor*(calcCLs(qmu95_guess, sigma, mu)-target_CLs)/calcDerCLs(qmu95_guess, sigma, mu);
      for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++)
      {
	if (fabs(itr->first - qmu95_guess) < 2*qmu95_guess*precision) 
	{
	  damping_factor *= 0.8;
	  if (verbose) cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
	  if (nrDamping++ > 10)
	  {
	    nrDamping = 1;
	    damping_factor = 1.0;
	  }
	  corr *= damping_factor;
	}
      }

      guess_to_corr[qmu95_guess] = corr;
      qmu95_guess = qmu95_guess - corr;

      if (verbose)
      {
	cout << "next guess = " << qmu95_guess << endl; 
	cout << "precision = " << 2*qmu95_guess*precision << endl;
	cout << endl;
      }
      nrItr++;
      if (nrItr > 200)
      {
	cout << "Infinite loop detected in getQmu95. Please intervene." << endl;
	exit(1);
      }
    }
    qmu95 = qmu95_guess;
  }

  if (qmu95 != qmu95) 
  {
    qmu95 = getQmu95_brute(sigma, mu);
  }
  if (verbose) cout << "Returning qmu95 = " << qmu95 << endl;

  return qmu95;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// calcCLs:

double calcCLs(double qmu_tilde, double sigma, double mu)
{
  double pmu = calcPmu(qmu_tilde, sigma, mu);
  double pb = calcPb(qmu_tilde, sigma, mu);
  if (verbose)
  {
    cout << "pmu = " << pmu << endl;
    cout << "pb = " << pb << endl;
  }
  if (pb == 1) return 0.5;
  return pmu/(1-pb);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// calcPmu:

double calcPmu(double qmu, double sigma, double mu)
{
  double pmu;
  if (qmu < mu*mu/(sigma*sigma) || !doTilde)
  {
    pmu = 1-ROOT::Math::gaussian_cdf(sqrt(qmu));
  }
  else
  {
    pmu = 1-ROOT::Math::gaussian_cdf((qmu+mu*mu/(sigma*sigma))/(2*fabs(mu/sigma)));
  }
  if (verbose) cout << "for pmu, qmu = " << qmu << ", sigma = " << sigma<< ", mu = " << mu << ", pmu = " << pmu << endl;
  return pmu;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// calcPb:

double calcPb(double qmu, double sigma, double mu)
{
  if (qmu < mu*mu/(sigma*sigma) || !doTilde)
  {
    return 1-ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qmu));
  }
  else
  {
    return 1-ROOT::Math::gaussian_cdf((mu*mu/(sigma*sigma) - qmu)/(2*fabs(mu/sigma)));
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// calcDerCLs:

double calcDerCLs(double qmu, double sigma, double mu)
{
  double dpmu_dq = 0;
  double d1mpb_dq = 0;

  if (qmu < mu*mu/(sigma*sigma))
  {
    double zmu = sqrt(qmu);
    dpmu_dq = -1./(2*sqrt(qmu*2*TMath::Pi()))*exp(-zmu*zmu/2);
  }
  else 
  {
    double zmu = (qmu+mu*mu/(sigma*sigma))/(2*fabs(mu/sigma));
    dpmu_dq = -1./(2*fabs(mu/sigma))*1./(sqrt(2*TMath::Pi()))*exp(-zmu*zmu/2);
  }

  if (qmu < mu*mu/(sigma*sigma))
  {
    double zb = fabs(mu/sigma)-sqrt(qmu);
    d1mpb_dq = -1./sqrt(qmu*2*TMath::Pi())*exp(-zb*zb/2);
  }
  else
  {
    double zb = (mu*mu/(sigma*sigma) - qmu)/(2*fabs(mu/sigma));
    d1mpb_dq = -1./(2*fabs(mu/sigma))*1./(sqrt(2*TMath::Pi()))*exp(-zb*zb/2);
  }

  double pb = calcPb(qmu, sigma, mu);
  return dpmu_dq/(1-pb)-calcCLs(qmu, sigma, mu)/(1-pb)*d1mpb_dq;
}

int minimize(RooNLLVar* nll)
{
  nrMinimize++;
  RooAbsReal* fcn = (RooAbsReal*)nll;
  return minimize(fcn);
}

int minimize(RooAbsReal* fcn)
{
  static int nrItr = 0;
  
  
  int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
  int save_strat = strat;
  RooMinimizer minim(*fcn);
  minim.setStrategy(strat);
  minim.setPrintLevel(printLevel);
  
  
  int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  
  
  //up the strategy
  if (status != 0 && status != 1 && strat < 2)
  {
    strat++;
    cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
    minim.setStrategy(strat);
    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  }

  if (status != 0 && status != 1 && strat < 2)
  {
    strat++;
    cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
    minim.setStrategy(strat);
    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  }
  

  // //switch minuit version and try again
  if (status != 0 && status != 1)
  {
    string minType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
    string newMinType;
    if (minType == "Minuit2") newMinType = "Minuit";
    else newMinType = "Minuit2";
  
    cout << "Switching minuit type from " << minType << " to " << newMinType << endl;
  
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(newMinType.c_str());
    strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    minim.setStrategy(strat);

    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


    if (status != 0 && status != 1 && strat < 2)
    {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    if (status != 0 && status != 1 && strat < 2)
    {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minType.c_str());
  }

  if (status != 0 && status != 1)
  {
    nrItr++;
    if (nrItr > maxRetries)
    {
      nrItr = 0;
      global_status++;
      cout << "WARNING::Fit failure unresolved with status " << status << endl;
      return status;
    }
    else
    {
      if (nrItr == 0) // retry with mu=0 snapshot
      {
	w->loadSnapshot("conditionalNuis_0");
	return minimize(fcn);
      }
      else if (nrItr == 1) // retry with nominal snapshot
      {
	w->loadSnapshot("nominalNuis");
	return minimize(fcn);
      }
    }
  }

  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);


  if (nrItr != 0) cout << "Successful fit" << endl;
  nrItr=0;
  return status;
}

void unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter)
{
  if (counter > 50)
  {
    cout << "ERROR::Couldn't unfold constraints!" << endl;
    cout << "Initial: " << endl;
    initial.Print("v");
    cout << endl;
    cout << "Final: " << endl;
    final.Print("v");
    exit(1);
  }
  TIterator* itr = initial.createIterator();
  RooAbsPdf* pdf;
  while ((pdf = (RooAbsPdf*)itr->Next()))
  {
    RooArgSet nuis_tmp = nuis;
    RooArgSet constraint_set(*pdf->getAllConstraints(obs, nuis_tmp, false));
    //if (constraint_set.getSize() > 1)
    //{
    string className(pdf->ClassName());
    if (className != "RooGaussian" && className != "RooLognormal" && className != "RooGamma" && className != "RooPoisson" && className != "RooBifurGauss")
    {
      counter++;
      unfoldConstraints(constraint_set, final, obs, nuis, counter);
    }
    else
    {
      final.add(*pdf);
    }
  }
  delete itr;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// makeAsimovData:

RooDataSet* makeAsimovData(bool doConditional, RooNLLVar* conditioning_nll, double mu_val, string* mu_str, string* mu_prof_str, double mu_val_profile, bool doFit)
{
  if (mu_val_profile == -999) mu_val_profile = mu_val;


  cout << "Creating asimov data at mu = " << mu_val << ", profiling at mu = " << mu_val_profile << endl;
  
  ////////////////////
  //make asimov data//
  ////////////////////
  RooAbsPdf* combPdf = mc->GetPdf();

  int _printLevel = 0;

  stringstream muStr;
  muStr << setprecision(5);
  muStr << "_" << mu_val;
  if (mu_str) *mu_str = muStr.str();

  stringstream muStrProf;
  muStrProf << setprecision(5);
  muStrProf << "_" << mu_val_profile;
  if (mu_prof_str) *mu_prof_str = muStrProf.str();

  RooRealVar* mu = (RooRealVar*)mc->GetParametersOfInterest()->first();//w->var("mu");
  mu->setVal(mu_val);

  RooArgSet mc_obs = *mc->GetObservables();
  RooArgSet mc_globs = *mc->GetGlobalObservables();
  RooArgSet mc_nuis = *mc->GetNuisanceParameters();

  //pair the nuisance parameter to the global observable
  RooArgSet mc_nuis_tmp = mc_nuis;
  RooArgList nui_list("ordered_nuis");
  RooArgList glob_list("ordered_globs");
  RooArgSet constraint_set_tmp(*combPdf->getAllConstraints(mc_obs, mc_nuis_tmp, false));
  RooArgSet constraint_set;
  int counter_tmp = 0;
  unfoldConstraints(constraint_set_tmp, constraint_set, mc_obs, mc_nuis_tmp, counter_tmp);

  TIterator* cIter = constraint_set.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)cIter->Next()))
  {
    RooAbsPdf* pdf = (RooAbsPdf*)arg;
    if (!pdf) continue;
    TIterator* nIter = mc_nuis.createIterator();
    RooRealVar* thisNui = NULL;
    RooAbsArg* nui_arg;
    while ((nui_arg = (RooAbsArg*)nIter->Next()))
    {
      if (pdf->dependsOn(*nui_arg))
      {
	thisNui = (RooRealVar*)nui_arg;
	break;
      }
    }
    delete nIter;
    
    
    //need this incase the observable isn't fundamental. 
    //in this case, see which variable is dependent on the nuisance parameter and use that.
    RooArgSet* components = pdf->getComponents();
    components->remove(*pdf);
    if (components->getSize())
    {
      TIterator* itr1 = components->createIterator();
      RooAbsArg* arg1;
      while ((arg1 = (RooAbsArg*)itr1->Next()))
      {
	TIterator* itr2 = components->createIterator();
	RooAbsArg* arg2;
	while ((arg2 = (RooAbsArg*)itr2->Next()))
	{
	  if (arg1 == arg2) continue;
	  if (arg2->dependsOn(*arg1))
	  {
	    components->remove(*arg1);
	  }
	}
	delete itr2;
      }
      delete itr1;
    }
    if (components->getSize() > 1)
    {
      cout << "ERROR::Couldn't isolate proper nuisance parameter" << endl;
      return NULL;
    }
    else if (components->getSize() == 1)
    {
      thisNui = (RooRealVar*)components->first();
    }



    TIterator* gIter = mc_globs.createIterator();
    RooRealVar* thisGlob = NULL;
    RooAbsArg* glob_arg;
    while ((glob_arg = (RooAbsArg*)gIter->Next()))
    {
      if (pdf->dependsOn(*glob_arg))
      {
	thisGlob = (RooRealVar*)glob_arg;
	break;
      }
    }
    delete gIter;

    if (!thisNui || !thisGlob)
    {
      cout << "WARNING::Couldn't find nui or glob for constraint: " << pdf->GetName() << endl;
      continue;
    }

    if (_printLevel >= 1) cout << "Pairing nui: " << thisNui->GetName() << ", with glob: " << thisGlob->GetName() << ", from constraint: " << pdf->GetName() << endl;

    nui_list.add(*thisNui);
    glob_list.add(*thisGlob);
  }
  delete cIter;
  
  
  
  //save the snapshots of nominal parameters, but only if they're not already saved
  w->saveSnapshot("tmpGlobs",*mc->GetGlobalObservables());
  w->saveSnapshot("tmpNuis",*mc->GetNuisanceParameters());
  if (!w->loadSnapshot("nominalGlobs"))
  {
    cout << "nominalGlobs doesn't exist. Saving snapshot." << endl;
    w->saveSnapshot("nominalGlobs",*mc->GetGlobalObservables());
  }
  else w->loadSnapshot("tmpGlobs");
  if (!w->loadSnapshot("nominalNuis"))
  {
    cout << "nominalNuis doesn't exist. Saving snapshot." << endl;
    w->saveSnapshot("nominalNuis",*mc->GetNuisanceParameters());
  }
  else w->loadSnapshot("tmpNuis");

  RooArgSet nuiSet_tmp(nui_list);

  mu->setVal(mu_val_profile);
  mu->setConstant(1);
  if (doConditional && doFit)
  {
    minimize(conditioning_nll);
  }
  mu->setConstant(0);
  mu->setVal(mu_val);



  //loop over the nui/glob list, grab the corresponding variable from the tmp ws, and set the glob to the value of the nui
  int nrNuis = nui_list.getSize();
  if (nrNuis != glob_list.getSize())
  {
    cout << "ERROR::nui_list.getSize() != glob_list.getSize()!" << endl;
    return NULL;
  }

  for (int i=0;i<nrNuis;i++)
  {
    RooRealVar* nui = (RooRealVar*)nui_list.at(i);
    RooRealVar* glob = (RooRealVar*)glob_list.at(i);
    glob->setVal(nui->getVal());
  }
  
  //save the snapshots of conditional parameters
  w->saveSnapshot(("conditionalGlobs"+muStrProf.str()).c_str(),*mc->GetGlobalObservables());
  w->saveSnapshot(("conditionalNuis" +muStrProf.str()).c_str(),*mc->GetNuisanceParameters());

  if (!doConditional)
  {
    w->loadSnapshot("nominalGlobs");
    w->loadSnapshot("nominalNuis");
  }

  if (_printLevel >= 1) cout << "Making asimov" << endl;
  //make the asimov data (snipped from Kyle)
  mu->setVal(mu_val);

  int iFrame=0;

  const char* weightName="weightVar";
  RooArgSet obsAndWeight;
  
  obsAndWeight.add(*mc->GetObservables());
  
  RooRealVar* weightVar = NULL;
  if (!(weightVar = w->var(weightName)))
  {
    w->import(*(new RooRealVar(weightName, weightName, 1,0,10000000)));
    weightVar = w->var(weightName);
  }
  obsAndWeight.add(*w->var(weightName));
  w->defineSet("obsAndWeight",obsAndWeight);
  
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  // MAKE ASIMOV DATA FOR OBSERVABLES

  RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());

  RooDataSet* asimovData;
  if (!simPdf)
  {
    // Get pdf associated with state from simpdf
    RooAbsPdf* pdftmp = mc->GetPdf();
	
    // Generate observables defined by the pdf associated with this state
    RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

    if (_printLevel >= 1)
    {
      obstmp->Print();
    }

    asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight),WeightVar(*weightVar));

    RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
    double expectedEvents = pdftmp->expectedEvents(*obstmp);
    double thisNorm = 0;
    for(int jj=0; jj<thisObs->numBins(); ++jj){
      thisObs->setBin(jj);

      thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
      if (thisNorm*expectedEvents <= 0)
      {
	cout << "WARNING::Detected bin with zero expected events (" << thisNorm*expectedEvents << ") ! Please check your inputs. Obs = " << thisObs->GetName() << ", bin = " << jj << endl;
      }
      if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18)) asimovData->add(*mc->GetObservables(), thisNorm*expectedEvents);
    }
    
    if (_printLevel >= 1)
    {
      asimovData->Print();
      cout <<"sum entries "<<asimovData->sumEntries()<<endl;
    }
    if(asimovData->sumEntries()!=asimovData->sumEntries()){
      cout << "sum entries is nan"<<endl;
      exit(1);
    }
         
    w->import(*asimovData);

    if (_printLevel >= 1)
    {
      asimovData->Print();
      cout << endl;
    }
  }
  else
  {
    map<string, RooDataSet*> asimovDataMap;

    
    //try fix for sim pdf
    RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;
    int nrIndices = 0;
    while((tt=(RooCatType*) iter->Next())) {
      nrIndices++;
    }
    for (int i=0;i<nrIndices;i++){
      channelCat->setIndex(i);
      iFrame++;
      // Get pdf associated with state from simpdf
      RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;
	
      // Generate observables defined by the pdf associated with this state
      RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

      if (_printLevel >= 1)
      {
	obstmp->Print();
	cout << "on type " << channelCat->getLabel() << " " << iFrame << endl;
      }

      RooDataSet* obsDataUnbinned = new RooDataSet(Form("combAsimovData%d",iFrame),Form("combAsimovData%d",iFrame),RooArgSet(obsAndWeight,*channelCat),WeightVar(*weightVar));
      RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
      double expectedEvents = pdftmp->expectedEvents(*obstmp);
      double thisNorm = 0;
      for(int jj=0; jj<thisObs->numBins(); ++jj){
	thisObs->setBin(jj);

	thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
	if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18)) obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
      }
    
      if (_printLevel >= 1)
      {
	obsDataUnbinned->Print();
	cout <<"sum entries "<<obsDataUnbinned->sumEntries()<<endl;
      }
      if(obsDataUnbinned->sumEntries()!=obsDataUnbinned->sumEntries()){
	cout << "sum entries is nan"<<endl;
	exit(1);
      }
      
      asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;//tempData;

      if (_printLevel >= 1)
      {
	cout << "channel: " << channelCat->getLabel() << ", data: ";
	obsDataUnbinned->Print();
	cout << endl;
      }
    }

    asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));
    w->import(*asimovData);
  }
  w->loadSnapshot("nominalGlobs");
  return asimovData;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// main:

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cout << "Usage: " << argv[0] << " <jobName> <DHSignal> <option>" << endl;
    return 0;
  }
  
  TString jobName = argv[1];
  inputDHSignal = argv[2];
  TString option = argv[3];// can be "highCL", "nosys"
  
  TString inputDir
    = Form("%s/%s", DHAnalysis::masterOutput.Data(), jobName.Data());
  TString inputFileName
    = Form("%s/DHWorkspace/rootfiles/workspaceDH.root", inputDir.Data());
  TString outputDir = Form("%s/DHMuLimit/single_files", inputDir.Data());
  
  // Make the output directory if it doesn't already exist:
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Copy the input file locally:
  TString localInputFileName = Form("workspaceDH.root", inputDHSignal.Data());
  system(Form("cp %s %s", inputFileName.Data(), localInputFileName.Data()));
  
  TString wname = "combinedWS";
  TString aname = "";//keep blank so that program creates its own Asimov data.
  
  // Model name and data name are dependent on the signal model tested:
  TString currAna = DHAnalysis::getAnalysisType(inputDHSignal);
  TString mname = Form("modelConfig_%s",currAna.Data());
  TString dname = Form("obsData_%s", currAna.Data());
  
  double CL = 0.95;
  if (option.Contains("highCL")) CL = 0.99;
  std::cout << "REGTEST: calculating " << CL*100.0 << "% CL limit" << std::endl;
  runAsymptoticsCLs(localInputFileName.Data(), wname.Data(), mname.Data(),
		    dname.Data(), aname.Data(), outputDir.Data(), CL, option);
  
  // Remove the local input file copy when job completes.
  system(Form("rm %s", localInputFileName.Data()));
}
