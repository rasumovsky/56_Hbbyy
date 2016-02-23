////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DHNP.cxx                                                                  //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 23/02/2016                                                          //
//                                                                            //
//  This program stores data for a single nuisance parameter from many fits.  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DHNP.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the DHNP class.
   @param nuisName- The name of the nuisance parameter.
*/
DHNP::DHNP(TString nuisName) {
  m_nuisName = nuisName;
  m_nuis["Mu0"] = new RooRealVar(Form("%sMu0",m_nuisName.Data()),
				 Form("%sMu0",m_nuisName.Data()), 0.0);
  m_nuis["Mu1"] = new RooRealVar(Form("%sMu1",m_nuisName.Data()),
				 Form("%sMu1",m_nuisName.Data()), 0.0);
  m_nuis["MuFree"] = new RooRealVar(Form("%sMuFree",m_nuisName.Data()),
				    Form("%sMuFree",m_nuisName.Data()), 0.0);
  m_uncertaintyOnPoI = 0.0;
}

/**
   -----------------------------------------------------------------------------
   Get the nuisance parameter name.
   @return - The name of the nuisance parameter.
*/
TString DHNP::getName() {
  return m_nuisName;
}

/**
   -----------------------------------------------------------------------------
   Get the nuisance parameter following the specified fit.
   @param fit - The name of the fit.
   @return - The post-fit nuisance parameter.
*/
RooRealVar *DHNP::getNuis(TString fit) {
  return m_nuis[fit];
}

/**
   -----------------------------------------------------------------------------
   Get the uncertainty on the PoI due to this parameter.
   @return - The uncertainty on the PoI due to this nuisance parameter.
*/
double DHNP::getPoIUncertainty() {
  return m_uncertaintyOnPoI;
}

/**
   -----------------------------------------------------------------------------
   Set the nuisance parameter for a given fit.
   @param fit - The name of the fit.
   @param nuisPar - A pointer to the nuisance parameter.
*/  
void DHNP::setNuis(TString fit, RooRealVar *nuisPar) {
  m_nuis[fit]
    = (RooRealVar*)nuisPar->clone(Form("%s%s", m_nuisName.Data(), fit.Data()));
  
}

/**
   -----------------------------------------------------------------------------
   Set the uncertainty on the PoI due to this parameter.
   @param uncertaintyOnPoI - The uncertainty on the parameter of interest.
*/  
void DHNP::setPoIUncertainty(double uncertaintyOnPoI) {
  m_uncertaintyOnPoI = uncertaintyOnPoI;  
}
