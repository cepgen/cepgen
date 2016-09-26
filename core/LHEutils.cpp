#include "LHEutils.h"

HEPRUP::HEPRUP(const int nprup_) :
  idwtup(2), nprup(nprup_)
{
  std::fill_n(this->idbmup, 2, 2212);
  std::fill_n(this->ebmup, 2, 3.5e3);
  std::fill_n(this->pdfgup, 2, 0);
  std::fill_n(this->pdfsup, 2, 10042);

  this->nprup = nprup_;

  xsecup = new double[nprup];
  xerrup = new double[nprup];
  xmaxup = new double[nprup];
  lprup = new int[nprup];

  for (int i=0; i<this->nprup; i++) {
    xsecup[i] = 0.;
    xerrup[i] = 0.;
    xmaxup[i] = 0.;
    lprup[i] = 0;
  }
}
HEPRUP::~HEPRUP()
{
  delete[] xsecup;
  delete[] xerrup;
  delete[] xmaxup;
  delete[] lprup;
}

HEPEUP::HEPEUP(const int nup_) :
  nup(nup_), idprup(0), xwgtup(1.), scalup(-1.), aqedup(1./137.), aqcdup(0.13)
{
  idup = new int[nup];
  istup = new int[nup];
  for (int i=0; i<2; i++) {
    mothup[i] = new int[nup];
    icolup[i] = new int[nup];
  }
  for (int p=0; p<5; p++) {
    pup[p] = new double[nup];
  }
  vtimup = new double[nup];
  spinup = new double[nup];
}

HEPEUP::~HEPEUP()
{
  delete[] idup;
  delete[] istup;
  for (int i=0; i<2; i++) {
    delete[] mothup[i];
    delete[] icolup[i];
  }
  for (int p=0; p<5; p++) {
    delete[] pup[p];
  }
  delete[] vtimup;
  delete[] spinup;
}
