/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooBernstein.cxx 44982 2012-07-10 08:36:13Z moneta $
 * Authors:                                                                  *
 *   Kyle Cranmer
 *                                                                           *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Bernstein basis polynomials are positive-definite in the range [0,1].
// In this implementation, we extend [0,1] to be the range of the parameter.
// There are n+1 Bernstein basis polynomials of degree n.
// Thus, by providing N coefficients that are positive-definite, there 
// is a natural way to have well bahaved polynomail PDFs.
// For any n, the n+1 basis polynomials 'form a partition of unity', eg.
//  they sum to one for all values of x. See
// http://www.idav.ucdavis.edu/education/CAGDNotes/Bernstein-Polynomials.pdf
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>
#include "TMath.h"
#include "RooBernsteinM.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooArgList.h"

using namespace std;

ClassImp(RooBernsteinM);


//_____________________________________________________________________________
RooBernsteinM::RooBernsteinM()
{
}


//_____________________________________________________________________________
RooBernsteinM::RooBernsteinM(const char* name, const char* title, 
			     RooAbsReal& x, const RooArgList& coefList, RooConstVar *a, RooConstVar *b): 
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients","List of coefficients",this)
{
  bool init = false;
  // Constructor
  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "RooBernsteinM::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
	   << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    _coefList.add(*coef);
  }
  delete coefIter ;

  if (a && b) {
    init = true;
    m_xmin = a->getVal();
    m_xmax = b->getVal();
  }
  if (!init) {
    // JBdV, the full range, used to define the reduce var (x-xmin)/(xmax-xmin)
    RooRealVar *xx = (RooRealVar*)(&x);
    m_xmin = xx->getMin();
    m_xmax = xx->getMax();
  }

  //cout << "RooBernstein defined with range " << m_xmin << " " << m_xmax << endl;
}


//_____________________________________________________________________________
RooBernsteinM::RooBernsteinM(const RooBernsteinM& other, const char* name) :
  RooAbsPdf(other, name), 
  _x("x", this, other._x), 
  _coefList("coefList",this,other._coefList)
{
  m_xmin = other.getXMin();
  m_xmax = other.getXMax();

  //cout << "Creating a Bern from a mother " << m_xmin << " " << m_xmax << endl;
}


//_____________________________________________________________________________
Double_t RooBernsteinM::evaluate() const 
{

  //JBdVDouble_t xmin = _x.min();
  //JBdVDouble_t x = (_x - xmin) / (_x.max() - xmin); // rescale to [0,1]
  Double_t x = (_x-m_xmin) / (m_xmax - m_xmin);  
  Int_t degree = _coefList.getSize() - 1; // n+1 polys of degree n
  RooFIter iter = _coefList.fwdIterator();

  if(degree == 0) {

    return ((RooAbsReal *)iter.next())->getVal();

  } else if(degree == 1) {

    Double_t a0 = ((RooAbsReal *)iter.next())->getVal(); // c0
    Double_t a1 = ((RooAbsReal *)iter.next())->getVal()-a0; // c1 - c0
    return a1 * x + a0;

  } else if(degree == 2) {

    Double_t a0 = ((RooAbsReal *)iter.next())->getVal(); // c0
    Double_t a1 = 2 * (((RooAbsReal *)iter.next())->getVal() - a0); // 2 * (c1 - c0)
    Double_t a2 = ((RooAbsReal *)iter.next())->getVal() - a1 - a0; // c0 - 2 * c1 + c2
    return (a2 * x + a1) * x + a0;

  } else if(degree > 2) {

    Double_t t = x;
    Double_t s = 1 - x;

    Double_t result = ((RooAbsReal *)iter.next())->getVal() * s;    
    for(Int_t i = 1; i < degree; i++) {
      result = (result + t * TMath::Binomial(degree, i) * ((RooAbsReal *)iter.next())->getVal()) * s;
      t *= x;
    }
    result += t * ((RooAbsReal *)iter.next())->getVal(); 

    return result;
  }

  // in case list of arguments passed is empty
  return TMath::SignalingNaN();
}


//_____________________________________________________________________________
Int_t RooBernsteinM::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  // No analytical calculation available (yet) of integrals over subranges
  if (rangeName && strlen(rangeName)) {
    return 0 ;
  }

  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}


//_____________________________________________________________________________
Double_t RooBernsteinM::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;
  Double_t xmin = _x.min(rangeName); Double_t xmax = _x.max(rangeName);
  Int_t degree= _coefList.getSize()-1; // n+1 polys of degree n
  Double_t norm(0) ;

  RooFIter iter = _coefList.fwdIterator() ;
  Double_t temp=0;
  for (int i=0; i<=degree; ++i){
    // for each of the i Bernstein basis polynomials
    // represent it in the 'power basis' (the naive polynomial basis)
    // where the integral is straight forward.
    temp = 0;
    for (int j=i; j<=degree; ++j){ // power basis≈ß
      //JBdVtemp += pow(-1.,j-i) * TMath::Binomial(degree, j) * TMath::Binomial(j,i) / (j+1);
      // this line could be buggy.
      temp += pow(-1.,j-i) * TMath::Binomial(degree, j) * TMath::Binomial(j,i) / (j+1) * ( TMath::Power(xmax-m_xmin,j+1) - TMath::Power(xmin-m_xmin,j+1) ) / TMath::Power(m_xmax-m_xmin,j);
    }
    temp *= ((RooAbsReal*)iter.next())->getVal(); // include coeff
    norm += temp; // add this basis's contribution to total
  }
  //JBdVnorm *= xmax-xmin;

  return norm;
}
