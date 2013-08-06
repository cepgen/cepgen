#ifndef _UTILS_H
#define _UTILS_H

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "gnuplot.h"
#include "particle.h"
//#include "gamgam.h"

#define MAX_HISTOS 20
/** @brief Electromagnetic coupling constant \f$\alpha_{em}=\frac{e^2}{4\pi\epsilon_0\hbar c}\f$ */
#define alphaF 1./137.04
/** @brief \f$\frac{1}{(\hbar c)^2}~[\mathrm b^{-1}]\f$? */
#define muBarn 1./389.39
#define pi 3.1415926535897932384626434
//#define sconst 2.1868465E10
#define sconst 3.89351824E8
#define sconstb 2.1868465E10
//#define RANMAX 1e8

double GetMassFromPDGId(int);

/**
 * Computes the proton structure function (F.W Brasse et al., DESY 76/11 (1976),
 *   http://dx.doi.org/10.1016/0550-3213(76)90231-5)
 * @cite Brasse1976413
 */
bool PSF(double,double,double*,double*,double*);
/**
 * @brief Defines modified variables of integration to avoid peaks integrations
 *  (see paper from Vermaseren for details : Two-photon processes,
 *  Nucl.Phys.B229 (1983) 347-371)
 * Returns a set of two modified variables of integration to maintain
 *  the stability of the integrant. These two new variables are :
 * - \f$y_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\f$
 *   the new variable
 * - \f$\mathrm dy_{out} = x_{min}\left(\frac{x_{max}}{x_{min}}\right)^{exp}\log\frac{x_{min}}{x_{max}}\f$,
 *   the new variable's differential form
 * @brief Redefines the variables of integration in order to avoid the
 *  strong peaking of the integrant.
 * @param expo_ Exponant
 * @param xmin_ Minimal value of the variable
 * @param xmax_ Maximal value of the variable
 * @param out_ The new variable definition
 * @param dout_ The differential variant of the new variable definition
 * \note This method overrides the set of `mapxx` subroutines in ILPAIR,
 *  with a slight difference according to the sign of the
 *  \f$\mathrm dy_{out}\f$ parameter :
 *  - left unchanged :
 * > `mapw2`, `mapxq`, `mapwx`, `maps2`
 *  - opposite sign :
 * > `mapt1`, `mapt2`
 */
void Map(double,double,double,double*,double*);
void Mapla(double,double,int,double,double,double*,double*);
//void Symmetrise(double, double, double*, double*);

/**
 * @brief List of input parameters used to start and run the simulation
 *  job.
 * @note The default parameters are derived from GMUINI in LPAIR
 */
class InputParameters {
  public:
    InputParameters();
    ~InputParameters();
    /**
     * @brief Dumps the input parameters in the console
     */
    void Dump();
    /**
     * @brief Reads content from config file to load the variables
     * @param inFile_ Name of the configuration file to load
     */
    bool ReadConfigFile(std::string);
    /**
     * @brief Stores the full run configuration to an external config file
     * @param outFile_ Name of the configuration file to create
     */
    bool StoreConfigFile(std::string);
    int ncvg; // ??
    /** @brief Number of Vegas integrations */
    int itvg;
    /** @brief First incoming particle's momentum (in GeV/c) */
    double in1p;
    /** @brief Second incoming particle's momentum (in GeV/c) */
    double in2p;
    /**
     * The first incoming particle type and kind of interaction :
     * - 1 - electron,
     * - 2 - proton elastic,
     * - 3 - proton inelastic without parton treatment,
     * - 4 - proton inelastic in parton model
     * @brief First particle's mode
     * @note Was named PMOD in ILPAIR
     */
    int p1mod;
    /**
     * @brief Second particle's mode
     * @note Was named EMOD in ILPAIR
     */
    int p2mod;
    /**
     * The particle code of produced leptons :
     * - 11 - for \f$e^+e^-\f$ pairs
     * - 13 - for \f$\mu^+\mu^-\f$ pairs
     * - 15 - for \f$\tau^+\tau^-\f$ pairs
     * @brief PDG id of the outgoing leptons
     */
    int pair;
    /**
     * Set of cuts to apply on the outgoing leptons in order to restrain the
     * available kinematic phase space :
     * - 0 - No cuts at all (for the total cross section)
     * - 1 - Vermaserens' hypothetical detector cuts : for both leptons,
     *   + \f$\frac{|p_z|}{|\mathbf p|}\leq\f$ 0.75 and \f$p_T\geq\f$ 1 GeV, or
     *   + 0.75 \f$<\frac{|p_z|}{|\mathbf p|}\leq\f$ 0.95 and \f$p_z>\f$ 1 GeV,
     * - 2 - Cuts according to the provided parameters
     * @brief Set of cuts to apply on the outgoing leptons
     */
    int mcut;
    /** @brief Minimal transverse momentum of the outgoing leptons */
    double minpt;
    /** @brief Maximal transverse momentum of the outgoing leptons */
    double maxpt;
    /** @brief Minimal energy of the outgoing leptons */
    double minenergy;
    /** @brief Maximal energy of the outgoing leptons */
    double maxenergy;
    /** @brief Minimal polar angle \f$\theta\f$ of the outgoing leptons */
    double mintheta;
    /** @brief Maximal polar angle \f$\theta\f$ of the outgoing leptons */
    double maxtheta;
    double minmx;
    double maxmx;
    /**
     * @brief Maximal number of iterations to perform by VEGAS
     */
    int itmx;
    /**
     * @brief Are we generating events ? (true) or are we only computing the
     * cross-section ? (false)
     */
    bool generation;
    /**
     * @brief Are the events generated in this run to be stored in the output
     * file ?
     */
    bool store;
    /**
     * Enables or disables the production of control plots for several kinematic
     * quantities in this process
     * @brief Do we need control plots all along the process?
     */
    bool debug;
    /**
     * @brief Number of events already generated in this run
     */
    int ngen;
    /**
     * @brief The file in which to store the events generation's output
     */
    std::ofstream* file;
    std::ofstream* file_debug; //FIXME dropme!
    /**
     * List of Gnuplot objects which can be used to produce control plots
     * all along the cross-section determination and events generation process
     * @note Maximum number of these can be raised in the utils.h file, but pay
     * attention to the memory load since these Gnuplot objects are still under
     * development!
     * @brief Control plots objects
     */
    Gnuplot* plot[MAX_HISTOS];
    //GamGam* gamgam;
};

/**
 * @brief List of kinematic cuts to apply on the central and outgoing phase space.
 */
class Cuts {
  public:
    Cuts();
    ~Cuts();
    /** @brief Sets of cuts to apply on the final phase space */
    int mode;
    /** @brief Minimal transverse momentum of the single outgoing leptons */
    double ptmin;
    /** @brief Maximal transverse momentum of the single outgoing leptons */
    double ptmax;
    /** @brief Minimal energy of the central two-photons system */
    double emin;
    /** @brief Maximal energy of the central two-photons system */
    double emax;
    /** @brief Minimal polar (\f$\theta_\mathrm{min}\f$) angle of the outgoing leptons, expressed in degrees */
    double thetamin;
    /** @brief Maximal polar (\f$\theta_\mathrm{max}\f$) angle of the outgoing leptons, expressed in degrees */
    double thetamax;
    double mxmin;
    double mxmax;
};

#endif
