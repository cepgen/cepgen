#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <fstream>

//#include "particle.h"
#include "event.h"

#include "pythia6hadroniser.h"
#include "jetset7hadroniser.h"

/**
 * @brief List of parameters used to start and run the simulation job.
 * @note The default parameters are derived from GMUINI in LPAIR
 */
class Parameters {
  public:
    Parameters();
    ~Parameters();
    /**
     * Defines the range to cover in pseudo-rapidity for the outgoing leptons produced in this process. This method converts this range into a range in \f$\theta\f$, the polar angle.
     * @brief Sets the pseudo-rapidity range for the produced leptons
     * @param[in] etamin_ The minimal value of \f$\eta\f$ for the outgoing leptons
     * @param[in] etamax_ The maximal value of \f$\eta\f$ for the outgoing leptons
     */
    void SetEtaRange(double etamin_, double etamax_);
    /**
     * @brief Dumps the input parameters in the console
     */
    void Dump();
    /**
     * Reads the list of parameters to be used in this cross-section computation/events generation from an external input card.
     * @brief Reads content from config file to load the variables
     * @param[in] inFile_ Name of the configuration file to load
     * @return A boolean stating whether this input configuration file is correct or not
     */
    bool ReadConfigFile(std::string inFile_);
    /**
     * @brief Stores the full run configuration to an external config file
     * @param[in] outFile_ Name of the configuration file to create
     * @return A boolean stating whether this output configuration file is correctly written or not
     */
    bool StoreConfigFile(std::string outFile_);
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
     * The particle code of produced leptons, as defined by the PDG convention :
     * - 11 - for \f$e^+e^-\f$ pairs
     * - 13 - for \f$\mu^+\mu^-\f$ pairs
     * - 15 - for \f$\tau^+\tau^-\f$ pairs
     * @brief PDG id of the outgoing leptons
     */
    int pair;
    /**
     * Set of cuts to apply on the outgoing leptons in order to restrain the available kinematic phase space :
     * - 0 - No cuts at all (for the total cross section)
     * - 1 - Vermaserens' hypothetical detector cuts : for both leptons,
     *   + \f$\frac{|p_z|}{|\mathbf p|}\leq\f$ 0.75 and \f$p_T\geq\f$ 1 GeV/c,
     *   or
     *   + 0.75 \f$<\frac{|p_z|}{|\mathbf p|}\leq\f$ 0.95 and \f$p_z>\f$ 1 GeV/c,
     * - 2 - Cuts on both the outgoing leptons, according to the provided cuts parameters
     * - 3 - Cuts on at least one outgoing lepton, according to the provided cut parameters
     * @brief Set of cuts to apply on the outgoing leptons
     */
    int mcut;
    /**
     * Minimal transverse momentum cut to apply on the outgoing lepton(s)
     * @brief Minimal \f$p_T\f$ of the outgoing leptons
     */
    double minpt;
    /**
     * Maximal transverse momentum cut to apply on the outgoing lepton(s)
     * @brief Maximal \f$p_T\f$ of the outgoing leptons
     */
    double maxpt;
    /** @brief Minimal energy of the outgoing leptons */
    double minenergy;
    /** @brief Maximal energy of the outgoing leptons */
    double maxenergy;
    /** @brief Minimal polar angle \f$\theta\f$ of the outgoing leptons */
    double mintheta;
    /** @brief Maximal polar angle \f$\theta\f$ of the outgoing leptons */
    double maxtheta;
    /**
     * @brief Minimal value of \f$Q^2\f$, the internal photons lines' virtuality
     */
    double minq2;
    /**
     * @brief Maximal value of \f$Q^2\f$, the internal photons lines' virtuality
     */
    double maxq2;
    /**
     * Minimal mass of the outgoing proton remnants, \f$M_X\f$, in GeV/c\f${}^{2}\f$.
     * @brief Minimal \f$M_X\f$ of the outgoing proton remnants
     */
    double minmx;
    /**
     * Maximal mass of the outgoing proton remnants, \f$M_X\f$, in GeV/c\f${}^{2}\f$.
     * @brief Maximal \f$M_X\f$ of the outgoing proton remnants
     */
    double maxmx;
    int ncvg; // ??
    /**
     * @brief Maximal number of iterations to perform by VEGAS
     */
    int itvg;
    /**
     * @brief Maximal number of TREAT calls
     * @note Is it correctly implemented ?
     */
    int ntreat;
    /**
     * @brief Number of points to "shoot" in each integration bin by the algorithm
     */
    int npoints;
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
     * Enables or disables the production of control plots for several kinematic quantities in this process
     * @brief Do we need control plots all along the process?
     */
    bool debug;
    /**
     * @brief Maximal number of events to generate in this run
     */
    int maxgen;
    /**
     * @brief Number of events already generated in this run
     */
    int ngen;
    /**
     * @brief PDFLIB group to use
     */
    int gpdf;
    /**
     * @brief PDFLIB set to use
     */
    int spdf;
    /**
     * @brief Number of quarks
     */
    int qpdf;
    /**
     * @brief The file in which to store the events generation's output
     */
    std::ofstream* file;
    /**
     * List of Gnuplot objects which can be used to produce control plots all along the cross-section determination and events generation process
     * @note Maximum number of these can be raised in the utils.h file, but pay attention to the memory load since these Gnuplot objects are still under development!
     * @brief Control plots objects
     */
    bool symmetrise;
    /**
     * @brief The pointer to the last event produced in this run
     */
    Event* last_event;
    /**
     * @brief Hadronisation algorithm to use for the proton(s) remnants fragmentation
     */
    Hadroniser* hadroniser;
};

#endif
