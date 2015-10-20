#ifndef LHEutils_h
#define LHEutils_h

#include <iostream>
#include <iomanip>
//#include <cstdlib>
//#include <cmath>
//#include <fstream>

/**
 * @brief Generic user-process interface for events generator
 */

/**
 * @brief User-process run information
 */
class HEPRUP {
 public:
  HEPRUP(const int nprup_=1);
  ~HEPRUP();
  /**
   * @brief ID of the beam 1 and 2 particles according to the Particle Data Group convention
   */
  int idbmup[2];
  /**
   * @brief Energy in GeV of the beam 1 and 2 particles
   */
  double ebmup[2];
  /**
   * @brief Author group for beam 1 and 2, according to PDFLIB
   */
  int pdfgup[2];
  /**
   * @brief PDF set ID for beam 1 and 2, according to PDFLIB
   */
  int pdfsup[2];

  int idwtup;
  int nprup;
  double *xsecup;
  double *xerrup;
  double *xmaxup;
  int *lprup;
};

/**
 * @brief User-process event information
 */
class HEPEUP {
 public:
  HEPEUP(const int nup_=500);
  ~HEPEUP();
  /**
   * @brief Maximum number of particle entries
   */
  static const int maxnup = 500;
  /**
   * @brief Number of particle entries in this event
   */
  int nup;
  /**
   * @brief ID of the process in this event
   */
  int idprup;
  /**
   * @brief Event weight
   */
  double xwgtup;
  /**
   * @brief Scale of the event in GeV, as used for the calculation of PDFs
   */
  double scalup;
  /**
   * @brief QED coupling \f$\alpha_{\mathrm{QED}}\f$ used for this event
   */
  double aqedup;
  /**
   * @brief QCD coupling \f$\alpha_{\mathrm{QCD}}\f$ used for this event
   */
  double aqcdup;

  /**
   * @brief Particle ID according to the Particle Data Group convention
   */
  int* idup;
  /**
   * @brief Status code
   */
  int* istup;
  /**
   * @brief Index of first and last mother
   */
  int* mothup[2];
  /**
   * @brief Index for the colour flow line passing through the colour (resp. anti-colour) of the particle
   */
  int* icolup[2];
  /**
   * @brief Lab-frame momentum of the particle, in GeV
   */
  double* pup[5];
  /**
   * @brief Invariant lifetime \f$c\tau\f$ in mm
   */
  double* vtimup;
  /**
   * @brief Cosine of the angle between the spin-vector of the particle and the 3-momentum of the decaying particle, in the lab frame
   */
  double* spinup;
};

#endif
