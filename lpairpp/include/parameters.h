#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>

#include "../processes/processes.h"

#include "pythia6hadroniser.h"
#include "pythia8hadroniser.h"
#include "jetset7hadroniser.h"
#include "herwig6hadroniser.h"

/**
 * @brief List of parameters used to start and run the simulation job.
 * @note The default parameters are derived from GMUINI in LPAIR
 */
class Parameters {
  public:
    Parameters();
    ~Parameters();
    /**
     * Defines the range to cover in polar angle for the outgoing leptons produced in this process. This method converts this range into a range in rapidity.
     * @brief Sets the polar angle range for the produced leptons
     * @param[in] thetamin_ The minimal value of \f$\theta\f$ for the outgoing leptons
     * @param[in] thetamax_ The maximal value of \f$\theta\f$ for the outgoing leptons
     */
    void SetThetaRange(double thetamin_, double thetamax_);
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
    /** @brief First incoming particle's momentum (in \f$\text{GeV}/c\f$) */
    double in1p;
    /** @brief Second incoming particle's momentum (in \f$\text{GeV}/c\f$) */
    double in2p;
    /**
     * @brief First beam/primary particle's PDG identifier
     */
    ParticleId in1pdg;
    /**
     * @brief Second beam/primary particle's PDG identifier
     */
    ParticleId in2pdg;
    /**
     * The first incoming particle type and kind of interaction :
     * - 1 - electron,
     * - 2 - proton elastic,
     * - 3 - proton inelastic without parton treatment,
     * - 4 - proton inelastic in parton model
     * @brief First particle's mode
     * @note Was named PMOD/EMOD in ILPAIR
     */
    int remnant_mode;
    /**
     * The particle code of produced leptons, as defined by the PDG convention :
     * - 11 - for \f$e^+e^-\f$ pairs
     * - 13 - for \f$\mu^+\mu^-\f$ pairs
     * - 15 - for \f$\tau^+\tau^-\f$ pairs
     * @brief PDG id of the outgoing leptons
     */
    ParticleId pair;
    int process_mode;
    /**
     * Set of cuts to apply on the outgoing leptons in order to restrain the available kinematic phase space :
     * - 0 - No cuts at all (for the total cross section)
     * - 1 - Vermaserens' hypothetical detector cuts : for both leptons,
     *   + \f$\frac{|p_z|}{|\mathbf p|}\leq\f$ 0.75 and \f$p_T\geq 1~\text{GeV}/c\f$,
     *   or
     *   + 0.75 \f$<\frac{|p_z|}{|\mathbf p|}\leq\f$ 0.95 and \f$p_z> 1~\text{GeV}/c\f$,
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
    /** @brief Minimal pseudorapidity \f$\eta\f$ of the outgoing leptons */
    double mineta;
    /** @brief Maximal pseudorapidity \f$\eta\f$ of the outgoing leptons */
    double maxeta;
    /**
     * @brief Minimal value of \f$Q^2\f$, the internal photons lines' virtuality
     */
    double minq2;
    /**
     * @brief Maximal value of \f$Q^2\f$, the internal photons lines' virtuality
     */
    double maxq2;
    /**
     * Minimal mass of the outgoing proton remnants, \f$M_X\f$, in \f$\text{GeV}/c^{2}\f$.
     * @brief Minimal \f$M_X\f$ of the outgoing proton remnants
     */
    double minmx;
    /**
     * Maximal mass of the outgoing proton remnants, \f$M_X\f$, in \f$\text{GeV}/c^{2}\f$.
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
     * @brief Maximal number of trials for the hadronisation of the proton(s) remnants
     */
    int hadroniser_max_trials;
    /**
     * @brief The file in which to store the events generation's output
     */
    std::ofstream* file;
    /**
     * @brief Type of format the event will be stored into
     */
    std::string output_format;
    /**
     * @brief Do we want the events to be symmetrised with respect to the \f$z\f$-axis ?
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
    /**
     * @brief The process for which the cross-section will be computed and the events will be generated
     */
    Process* process;
};

#endif
