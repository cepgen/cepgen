#ifndef Parameters_h
#define Parameters_h

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>

#include "processes/GamGamLL.h"
#include "processes/PPtoLL.h"
#include "processes/PPtoWW.h"

#include "hadronisers/Pythia6Hadroniser.h"
#include "hadronisers/Jetset7Hadroniser.h"
#include "hadronisers/Herwig6Hadroniser.h"

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
    void SetThetaRange( float thetamin_, float thetamax_ );
    /**
     * @brief Dumps the input parameters in the console
     */
    void Dump();
    /// Read content from config file to load the variables
    /// \param[in] inFile_ Name of the configuration file to load
    /// \return A boolean stating whether this input configuration file is correct or not
    bool ReadConfigFile( const char* inFile_ );
    /// Store the full run configuration to an external config file
    /// \param[in] outFile_ Name of the configuration file to create
    /// \return A boolean stating whether this output configuration file is correctly written or not
    bool StoreConfigFile( const char* outFile_ );

    //----- process to compute
    /// Process for which the cross-section will be computed and the events will be generated
    GenericProcess* process;
    /// Type of outgoing state to consider for the incoming primary particles
    Kinematics::ProcessMode process_mode;

    /// Type of remnant fragmentation algorithm to use
    /// \note Was named PMOD/EMOD in ILPAIR
    GenericProcess::StructureFunctions remnant_mode;

    //----- events kinematics
    /// First incoming particle's momentum (in \f$\text{GeV}/c\f$)
    float in1p;
    /// Second incoming particle's momentum (in \f$\text{GeV}/c\f$)
    float in2p;
    /// First beam/primary particle's PDG identifier
    Particle::ParticleCode in1pdg;
    /// Second beam/primary particle's PDG identifier
    Particle::ParticleCode in2pdg;
    /// PDG id of the outgoing central particles
    Particle::ParticleCode pair;
    /// Set of cuts to apply on the outgoing central system
    Kinematics::Cuts mcut;
    /// Minimal \f$p_T\f$ of the outgoing central particles
    float minpt;
    /// Maximal \f$p_T\f$ of the outgoing central particles
    float maxpt;
    /// Minimal \f$\Delta p_T\f$ of the outgoing central particles
    float minptdiff;
    /// Maximal \f$\Delta p_T\f$ of the outgoing central particles
    float maxptdiff;
    /// Minimal energy of the outgoing central particles
    float minenergy;
    /// Maximal energy of the outgoing central particles
    float maxenergy;
    /// Minimal pseudorapidity \f$\eta\f$ of the outgoing central particles
    float mineta;
    /// Maximal pseudorapidity \f$\eta\f$ of the outgoing central particles
    float maxeta;
    float minqt, maxqt;
    /// Minimal value of \f$Q^2\f$, the internal photons lines' virtuality
    float minq2;
    /// Maximal value of \f$Q^2\f$, the internal photons lines' virtuality
    float maxq2;
    /// Minimal \f$M_X\f$ of the outgoing proton remnants
    float minmx;
    /// Maximal \f$M_X\f$ of the outgoing proton remnants
    float maxmx;

    //----- VEGAS
    int ncvg; // ??
    /// Maximal number of iterations to perform by VEGAS
    int itvg;
    /// Number of points to "shoot" in each integration bin by the algorithm
    int npoints;
    /// Is it the first time the integrator is run?
    bool first_run;

    //----- events generation
    /// Are we generating events ? (true) or are we only computing the cross-section ? (false)
    bool generation;
    /// Are the events generated in this run to be stored in the output file ?
    bool store;
    /// Maximal number of events to generate in this run
    int maxgen;
    /// Pointer to the last event produced in this run
    Event* last_event;
    /// File in which to store the events generation's output
    std::ofstream* file;
    /// Do we want the events to be symmetrised with respect to the \f$z\f$-axis ?
    bool symmetrise;

    //----- PDFLIB information
    /// Number of events already generated in this run
    unsigned int ngen;
    /// PDFLIB group to use
    unsigned int gpdf;
    /// PDFLIB set to use
    unsigned int spdf;
    /// Number of quarks to consider in the hadronisation part
    unsigned int qpdf;

    //----- hadronisation
    /// Hadronisation algorithm to use for the proton(s) fragmentation
    GenericHadroniser* hadroniser;
    /// Maximal number of trials for the hadronisation of the proton(s) remnants
    int hadroniser_max_trials;
};

#endif
