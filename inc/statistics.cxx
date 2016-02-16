////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Author: Haichen Wang                                                      //
//  Copyright (C) 2008                                                        //
//  All Rights Reserved                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "statistics.h"

using namespace RooFit;
using namespace RooStats;

ClassImp(statistics)

statistics::statistics(){
  //other configuration
  defaultMinimizer    = "Minuit2";     // or "Minuit"
  defaultPrintLevel      = -1;             // Minuit print level
  defaultStrategy        = 0;             // Minimization strategy. 0-2. 0 = fastest, least robust. 2 = slowest, most robust
  killBelowFatal        = 1;             // In case you want to suppress RooFit warnings further, set to 1
  doBlind               = 0;             // in case your analysis is blinded
  conditionalExpected   = 1 && !doBlind; // Profiling mode for Asimov data: 0 = conditional MLEs, 1 = nominal MLEs
  doTilde               = 1;             // bound mu at zero if true and do the \tilde{q}_{mu} asymptotics
  doExp                 = 1;             // compute expected limit
  doObs                 = 1 && !doBlind; // compute observed limit
  precision           = 0.005;         // % precision in mu that defines iterative cutoff
  verbose               = 0;             // 1 = very spammy
  usePredictiveFit      = 0;             // experimental, extrapolate best fit nuisance parameters based on previous fit results
  extrapolateSigma      = 0;             // experimantal, extrapolate sigma based on previous fits
  maxRetries             = 3;             // number of minimize(fcn) retries before giving up
  fRandom=new TRandom3();
  if (killBelowFatal) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(defaultMinimizer.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(defaultStrategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(defaultPrintLevel);
}

statistics::~statistics(){
  SafeDelete(fRandom);
}

RooNLLVar* statistics::createNLL(RooAbsData* _data, ModelConfig* _mc)
{
  RooArgSet nuis = *_mc->GetNuisanceParameters();
  RooNLLVar* nll = (RooNLLVar*)_mc->GetPdf()->createNLL(*_data, Constrain(nuis), Extended(_mc->GetPdf()->canBeExtended()));
  return nll;
}

RooFitResult* statistics::minimize(RooAbsReal* fcn, TString option, RooArgSet *minosVars, bool m_save)
{
  option.ToLower();

  bool doHesse=option.Contains("hesse");
  bool doMinos=option.Contains("minos");
  bool m_2sigma=option.Contains("2sigma");

  int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
  int save_strat = strat;
  RooMinimizer minim(*fcn);
  minim.setStrategy(strat);
  minim.setPrintLevel(printLevel);
  if(m_2sigma) minim.setErrorLevel(2);
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

  cout << "status is " << status << endl;

//switch minuit version and try again
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
    cout << "WARNING::Fit failure unresolved with status " << status << endl;
    return NULL;
  }
  
  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);

  if(doHesse) minim.hesse();
  if(doMinos){
    if(minosVars==NULL) minim.minos();
    else minim.minos(*minosVars);
  }
  if(m_save) return minim.save();
  else return NULL;
}

RooFitResult* statistics::minimize(RooNLLVar* nll, TString option, RooArgSet *minosVars)
{
  RooAbsReal* fcn = (RooAbsReal*)nll;
  return minimize(fcn,option,minosVars);
}

void statistics::recoverSet(RooArgSet* set, RooArgSet* snapshot){
  TIterator *iter = set -> createIterator();
  RooRealVar* parg = NULL;
  TIterator *iter_ss = snapshot -> createIterator();
  RooRealVar* parg_ss = NULL;  
  while((parg_ss=(RooRealVar*)iter_ss->Next())&&(parg=(RooRealVar*)iter->Next())){
    parg->setVal(parg_ss->getVal());
  }
}

void statistics::constSet(RooArgSet* set, bool flag, RooArgSet* snapshot){
  TIterator *iter = set -> createIterator();
  RooRealVar* parg = NULL;
//   if(snapshot!=NULL) recoverSet(set,snapshot);
  if(snapshot!=NULL) *set=*snapshot;
  while((parg=(RooRealVar*)iter->Next()) ){
    parg->setConstant(flag); 
  }
  SafeDelete(iter);
}

void statistics::randomizeSet(RooArgSet* set, int seed, bool protection){
  TIterator *iter = set -> createIterator();
  RooRealVar* parg = NULL;

  fRandom -> SetSeed(seed) ;   
  while((parg=(RooRealVar*)iter->Next()) ){
    double r =  fRandom->Gaus(0,1);
    if(protection&&r<-5) r=-5;
    if(protection&&r>5) r=5;
    parg->setVal(r); 
  }
  SafeDelete(iter);
}

void statistics::retrieveSet(RooWorkspace*w, RooArgSet* set, RooArgSet* snapshot){
  TIterator *iter = snapshot -> createIterator();
  RooRealVar* parg = NULL;
  while((parg=(RooRealVar*)iter->Next())){
    if((bool)w->obj(parg->GetName()))
      set -> add( *(RooRealVar*)w->obj(parg->GetName()) );
  }
}

RooDataSet* statistics::histToDataSet(TH1* h, RooRealVar* x, RooRealVar* w, RooCategory* c){
  double xmin=x->getMin();
  double xmax=x->getMax();
  RooDataSet *histData =new RooDataSet(h->GetName(),h->GetTitle(),RooArgSet(*x,*w),WeightVar(*w));
  int nbin=h->GetNbinsX();
  for( int ibin = 1 ; ibin <= nbin ; ibin ++ ) {
    double lowedge=h->GetBinLowEdge(ibin);
    double highedge=h->GetBinLowEdge(ibin)+h->GetBinWidth(ibin);
    if(highedge<xmin) continue;
    if(lowedge>xmax) continue;
    x->setVal(h->GetBinCenter(ibin));
    double weight = h->GetBinContent(ibin);
    double weighterr = h->GetBinError(ibin);
    weighterr=0;
    w->setVal(weight);
    histData -> add( RooArgSet(*x,*w) , weight, weighterr);
  }
  histData->Print(); 
  cout<<c<<endl;
  if(c){
    cout<<"creating channel list"<<endl;
    map<string,RooDataSet*> datamap;
    datamap[c->getLabel()]=histData;
    histData=new RooDataSet(h->GetName(),h->GetTitle(),RooArgSet(*x,*w),Index(*c),Import(datamap),WeightVar(*w));
  }
  return histData;
}

void statistics::randomizeSet(RooAbsPdf* pdf, RooArgSet* globs, int seed){
  if(seed>=0) RooRandom::randomGenerator() -> SetSeed(seed) ; // This step is necessary
  RooDataSet *one=pdf->generate(*globs, 1);
  const RooArgSet *values=one->get(0);
  RooArgSet *allVars=pdf->getVariables();
  *allVars=*values;
  delete one;
  delete allVars;
}

double statistics::pvalueError(double pvalue, int ntoy){
  return sqrt(pvalue*(1-pvalue)/double(ntoy));
}

double statistics::pvalueFromToy(vector<double> teststat, double threshold){
  std::sort(teststat.begin(),teststat.end());
  int ntoy=teststat.size();
  int marker=-1;

  for(int itoy=0;itoy<ntoy;itoy++){
//     cout<<toy[itoy]<<endl;
//     if(!isfinite(toy[itoy])) cout<<"NaN!"<<endl; 
    if(teststat[itoy]>threshold){ marker=itoy; break;}
  }
  double pvalue=double(marker)/double(ntoy);
  if(pvalue>0.5) pvalue=1-pvalue;
  return pvalue;
}

map<string,double> statistics::expFromToy(vector<double> teststat){
  double median=0;
  double mean=0;
  double n2s=0;
  double n1s=0;
  double p1s=0;
  double p2s=0;

  int num = (int)teststat.size();

  sort(teststat.begin(),teststat.end());

  for(int index=0; index<num; index++){
    double frac = double(index)/double(num);
    if (frac < LB2S) n2s = teststat[index];
    if (frac < LB1S) n1s = teststat[index];
    if (frac < B0S) median = teststat[index];
    if (frac < UB1S) p1s = teststat[index];
    if (frac < UB2S) p2s = teststat[index];
  }
  map<string,double> results;
  results["n2s"]=n2s;
  results["n1s"]=n1s;
  results["median"]=median;
  results["p1s"]=p1s;
  results["p2s"]=p2s;

  return results;
}
