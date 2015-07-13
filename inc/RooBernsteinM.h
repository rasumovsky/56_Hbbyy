/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooBernstein.h 28259 2009-04-16 16:21:16Z wouter $
 * Authors:                                                                  *
 *   Kyle Cranmer
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_BERNSTEINM
#define ROO_BERNSTEINM

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooConstVar;
class RooArgList ;

class RooBernsteinM : public RooAbsPdf {
public:

  RooBernsteinM() ;
  RooBernsteinM(const char *name, const char *title,
		RooAbsReal& _x, const RooArgList& _coefList, RooConstVar* a = 0, RooConstVar* b = 0) ;

  RooBernsteinM(const RooBernsteinM& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooBernsteinM(*this, newname); }
  inline virtual ~RooBernsteinM() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  Double_t getXMin() const { return m_xmin; }
  Double_t getXMax() const { return m_xmax; }

private:

  RooRealProxy _x;
  RooListProxy _coefList ;

  Double_t evaluate() const;

  // JBdV
  Double_t m_xmin, m_xmax;

  ClassDef(RooBernsteinM,1) // Bernstein polynomial PDF

};

#endif
