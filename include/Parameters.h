#ifndef Parameters_h
#define Parameters_h

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>

#include "../processes/processes.h"
#include "../hadronisers/hadronisers.h"

/// List of parameters used to start and run the simulation job
/// \note The default parameters are derived from GMUINI in LPAIR
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
    /// Read content from config file to load the variables
    /// \param[in] inFile_ Name of the configuration file to load
    /// \return A boolean stating whether this input configuration file is correct or not
    bool ReadConfigFile(const char* inFile_);
    /// Store the full run configuration to an external config file
    /// \param[in] outFile_ Name of the configuration file to create
    /// \return A boolean stating whether this output configuration file is correctly written or not
    bool StoreConfigFile(std::string outFile_);
    /// First incoming particle's momentum (in \f$\text{GeV}/c\f$)
    double in1p;
    /// Second incoming particle's momentum (in \f$\text{GeV}/c\f$)
    double in2p;
    /// First beam/primary particle's PDG identifier
    Particle::ParticleCode in1pdg;
    /// Second beam/primary particle's PDG identifier
    Particle::ParticleCode in2pdg;
    /// Type of remnant fragmentation algorithm to use
    /// \note Was named PMOD/EMOD in ILPAIR
    GenericProcess::StructureFunctions remnant_mode;
    /**
     * The particle code of produced leptons, as defined by the PDG convention :
     * - 11 - for \f$e^+e^-\f$ pairs
     * - 13 - for \f$\mu^+\mu^-\f$ pairs
     * - 15 - for \f$\tau^+\tau^-\f$ pairs
     * @brief PDG id of the outgoing leptons
     */
    Particle::ParticleCode pair;
    /// Type of outgoing state to consider for the incoming primary particles
    GenericProcess::ProcessMode process_mode;
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
    /// Minimal \f$p_T\f$ of the outgoing central particles
    double minpt;
    /// Maximal \f$p_T\f$ of the outgoing central particles
    double maxpt;
    /// Minimal \f$\Delta p_T\f$ of the outgoing central particles
    double minptdiff;
    /// Maximal \f$\Delta p_T\f$ of the outgoing central particles
    double maxptdiff;
    /// Minimal energy of the outgoing central particles
    double minenergy;
    /// Maximal energy of the outgoing central particles
    double maxenergy;
    /// Minimal pseudorapidity \f$\eta\f$ of the outgoing central particles
    double mineta;
    /// Maximal pseudorapidity \f$\eta\f$ of the outgoing central particles
    double maxeta;
    double minqt, maxqt;
    /// Minimal value of \f$Q^2\f$, the internal photons lines' virtuality
    double minq2;
    /// Maximal value of \f$Q^2\f$, the internal photons lines' virtuality
    double maxq2;
    /// Minimal \f$M_X\f$ of the outgoing proton remnants
    double minmx;
    /// Maximal \f$M_X\f$ of the outgoing proton remnants
    double maxmx;
    int ncvg; // ??
    /// Maximal number of iterations to perform by VEGAS
    int itvg;
    /// Number of points to "shoot" in each integration bin by the algorithm
    int npoints;
    /// Are we generating events ? (true) or are we only computing the cross-section ? (false)
    bool generation;
    /// Are the events generated in this run to be stored in the output file ?
    bool store;
    /// Maximal number of events to generate in this run
    int maxgen;
    /// Number of events already generated in this run
    int ngen;
    /// PDFLIB group to use
    int gpdf;
    /// PDFLIB set to use
    int spdf;
    /// Number of quarks to consider in the hadronisation part
    int qpdf;
    /// Maximal number of trials for the hadronisation of the proton(s) remnants
    int hadroniser_max_trials;
    /// File in which to store the events generation's output
    std::ofstream* file;
    /// Type of format the event will be stored into
    std::string output_format;
    /// Do we want the events to be symmetrised with respect to the \f$z\f$-axis ?
    bool symmetrise;
    /// Is it the first time the integrator is run?
    bool first_run;
    /// Pointer to the last event produced in this run
    Event* last_event;
    /// Process for which the cross-section will be computed and the events will be generated
    GenericProcess* process;
    /// Hadronisation algorithm to use for the proton(s) fragmentation
    GenericHadroniser* hadroniser;
};

#endif
