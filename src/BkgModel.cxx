////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: BkgModel.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/06/2015                                                          //
//                                                                            //
//  This class implements a broad range of fit functions to use for the       //
//  background model in the H->yy analysis.                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "BkgModel.h"

/**
   -----------------------------------------------------------------------------
   Initialize the BkgModel class and make a new RooCategory.
   @param newOptions - The job options ("New", "FromFile")
   @param newObservable - The RooRealVar to be used in fits (m_obs).
   @returns void.
*/
BkgModel::BkgModel(RooRealVar *newObs) {
  std::cout << "\nBkgModel::Initializing..." << std::endl;
  
  // Assign the observable based on inputs:
  if (!newObs) newObs = new RooRealVar("m_yy", "m_yy", 105.0, 160.0);
  setObservable(newObs);
  
  std::cout << "\nBkgModel::Initialized!" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Add the chosen background model to the workspace provided. Also add the 
   associated nuisance parameters to the nuisParams set.
   @param workspace - The workspace to which the PDFs will be added.
   @param nuisParams - The set of nuisance parameters to which the background
                       parameters will be added.
   @param function - the functional form (PDF) for the model.
*/
void BkgModel::addBkgToCateWS(RooWorkspace *&workspace, RooArgSet *&nuisParams,
			      TString function) {
  std::cout << "BkgModel: Adding " << function
	    << " background model to the workspace" << std::endl;
  
  // First, call the function to get the background PDF:
  RooAbsPdf* currBkgModel = getBkgPDF(function);
  
  // Then add it to the workspace:
  workspace->import(*currBkgModel);
  
  // Then add the parameters to the workspace:
  RooArgSet *currArgs = currBkgModel->getVariables();
  TIterator *iterArgs = currArgs->createIterator();
  RooRealVar* currIter = NULL;
  while ((currIter = (RooRealVar*)iterArgs->Next())) nuisParams->add(*currIter);
  
  // Finally, include a normalization parameter for the background:
  workspace->factory("nBkg[100,0,1000000]");
  nuisParams->add(*workspace->var("nBkg"));
  
  std::cout << "BkgModel: Finished adding background model" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Get the background PDF using name of function.
   @param function - the functional form (PDF) for the model.
   @returns The corresponding background PDF for the analysis.
*/
RooAbsPdf* BkgModel::getBkgPDF(TString function) {
  std::cout << "BkgModel: Constructing " << function << std::endl;
  
  int order = getOrderFromFunc(function);
  
  // Set the range of the m_obs variable from DMHeader.h:
  RooConstVar min("min", "min", m_obsMin);
  RooConstVar max("max", "max", m_obsMax);
  
  // Pointers to the background PDFs:
  RooAbsPdf *background = NULL;
  RooBernsteinM *bern = NULL;
  RooGenericPdf *exppol = NULL;
  RooLandau *landau = NULL;
  
  // Background fit variables:
  RooRealVar *pVar[10];
  RooRealVar *cVar[10];
  RooRealVar *lVar[2];
  RooArgList *bkgArgs = new RooArgList();
  if (function.Contains("Exppol")) bkgArgs->add(*m_obs);

  TString expFitFormat = "exp(";
  
  // Loop over the order of the function to define parameters:
  for (int i_p = 0; i_p <= order; i_p++) {
    // Parameters for Bernstein polynomial:
    if (function.Contains("Bern")) {
      if (i_p == 0) {
	pVar[i_p] = new RooRealVar(Form("pVar%d",i_p), 
				   Form("pVar%d",i_p), 1);
	pVar[i_p]->setConstant(true);
      }
      else {
	pVar[i_p] = new RooRealVar(Form("pVar%d",i_p), Form("pVar%d",i_p),
				   0.1, 0.0, 10.0);
	pVar[i_p]->setConstant(false);
      }
      bkgArgs->add(*pVar[i_p]);
    }
    // Parameters for exponential polynomial:
    else if (function.Contains("Exppol") && i_p < order) {
      cVar[i_p] = new RooRealVar(Form("cVar%d",i_p), Form("cVar%d",i_p),
				 0.0, -1.0, 1.0 );
      cVar[i_p]->setConstant(false);
      bkgArgs->add(*cVar[i_p]);
      expFitFormat += Form("@%d",i_p+1);
      for (int i_t = 0; i_t <= i_p; i_t++) {
	expFitFormat += "*(@0-100)";
      }
    }
  }
  expFitFormat += ")";
  
  std::cout << "Printing bkgArgs:" << std::endl;
  bkgArgs->Print("v");

  // Construct the desired PDF:
  if (function.Contains("Bern")) {
    bern = new RooBernsteinM("bkgPdf", "bkgPdf", *m_obs, *bkgArgs, &min, &max);
    background = bern;
  }
  else if (function.Contains("Exppol")) {
    exppol = new RooGenericPdf("bkgPdf", expFitFormat, *bkgArgs);
    background = exppol;
  }
  else if (function.EqualTo("Landau")) {
    lVar[0] = new RooRealVar("lVar0", "lVar0", 300.0, 200.0, 1000.0);
    lVar[1] = new RooRealVar("lVar1", "lVar1", 50.0, 10.0, 200.0);
    lVar[0]->setConstant(false);
    lVar[1]->setConstant(false);
    landau = new RooLandau("bkgPdf","bkgPdf", *m_obs, *lVar[0], *lVar[1]);
    background = landau;
  }
  
  // Returns the constructed PDF:
  std::cout << "BkgModel: Finished constructing " << function << std::endl;
  return background;
}

/**
   -----------------------------------------------------------------------------
   Returns a pointer to the mass observable used in the dataset.
   @returns pointer to the observable (m_obs).
*/
RooRealVar* BkgModel::getMassObservable() {
  return m_obs;
}

/**
   -----------------------------------------------------------------------------
   Set the pointer to the observable. 
   @param newObs - The new RooRealVar observable to use for datasets. 
   @returns void.
*/
void BkgModel::setObservable(RooRealVar *newObs) {
  m_obs = newObs;
  m_obsName = m_obs->GetName();
  m_obsMin = m_obs->getMin();
  m_obsMax = m_obs->getMax();
  std::cout << "BkgModel: observable = " << m_obsName << "[" << m_obsMin << ", "
	    << m_obsMax << "]" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Get the order of the function from the name.
   @param function - the name of the function.
   @returns - the order of the function.
*/
int BkgModel::getOrderFromFunc(TString function) {
  if (function.EqualTo("Landau")) return 1;
  for (int i_o = 0; i_o < 10; i_o++) {
    if (function.Contains(Form("O%d",i_o))) return i_o;
  }
  return 0;
}
