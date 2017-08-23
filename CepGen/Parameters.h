#ifndef CepGen_Parameters_h
#define CepGen_Parameters_h

#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Core/Functional.h"

#include <memory>

namespace CepGen
{
  /// List of parameters used to start and run the simulation job
  /// \note The default parameters are derived from GMUINI in LPAIR
  class Parameters {
    public:
      Parameters();
      /// Copy constructor (transfers ownership to the process!)
      Parameters( Parameters& );
      /// Const copy constructor (all but the process and the hadroniser)
      Parameters( const Parameters& );
      ~Parameters();
      /// Set the polar angle range for the produced leptons
      /// \param[in] thetamin The minimal value of \f$\theta\f$ for the outgoing leptons
      /// \param[in] thetamax The maximal value of \f$\theta\f$ for the outgoing leptons
      void setThetaRange( float thetamin, float thetamax );
      /// Dump the input parameters in the console
      void dump( std::ostream& os=Logger::get().outputStream, bool pretty=true ) const;

      //----- process to compute

      /// Process for which the cross-section will be computed and the events will be generated
      Process::GenericProcess* process() { return process_.get(); }
      std::string processName() const { return process_->name(); }
      /// Set the process to study
      void setProcess( Process::GenericProcess* proc ) { process_.reset( proc ); }

      /// Type of remnant fragmentation algorithm to use
      /// \note Was named PMOD/EMOD in ILPAIR
      StructureFunctions remnant_mode;

      //----- events kinematics

      Kinematics kinematics;

      //----- VEGAS

      /// Collection of Vegas integrator parameters
      struct Vegas
      {
        Vegas() : ncvg( 100000 ), itvg( 10 ), npoints( 100 ), first_run( true ) {}
        unsigned int ncvg; // ??
        /// Maximal number of iterations to perform by VEGAS
        unsigned int itvg;
        /// Number of points to "shoot" in each integration bin by the algorithm
        unsigned int npoints;
        /// Is it the first time the integrator is run?
        bool first_run;
      };
      Vegas vegas;

      //----- events generation

      struct Generation
      {
        Generation() : enabled( false ), maxgen( 0 ), symmetrise( false ), ngen( 0 ), gen_print_every( 1 ) {}
        /// Are we generating events ? (true) or are we only computing the cross-section ? (false)
        bool enabled;
        /// Maximal number of events to generate in this run
        unsigned int maxgen;
        /// Pointer to the last event produced in this run
        std::shared_ptr<Event> last_event;
        /// Do we want the events to be symmetrised with respect to the \f$z\f$-axis ?
        bool symmetrise;
        /// Number of events already generated in this run
        unsigned int ngen;
        /// Frequency at which the events are displayed to the end-user
        unsigned int gen_print_every;
      };
      Generation generation;

      void setStorage( bool store ) { store_ = store; }
      /// Are the events generated in this run to be stored in the output file ?
      bool storage() const { return store_; }

      //----- PDFLIB information

      /// Collection of PDFLIB steering parameters
      struct PDFLIB
      {
        PDFLIB() : gpdf( 5 ), spdf( 4 ), qpdf( 12 ) {}
        /// PDFLIB group to use
        unsigned int gpdf;
        /// PDFLIB set to use
        unsigned int spdf;
        /// Number of quarks to consider in the hadronisation part
        unsigned int qpdf;
      };
      PDFLIB pdflib;

      //----- hadronisation

      /// Hadronisation algorithm to use for the proton(s) fragmentation
      Hadroniser::GenericHadroniser* hadroniser() { return hadroniser_.get(); }
      /// Set the hadronisation algorithm
      void setHadroniser( Hadroniser::GenericHadroniser* hadr ) { hadroniser_.reset( hadr ); }
      /// Maximal number of trials for the hadronisation of the proton(s) remnants
      unsigned int hadroniser_max_trials;

      //----- rescattering corrections
      /// Functionals to be used to account for rescattering corrections (implemented within the process)
      std::map<std::string,Functional<1> > taming_functions;

    private:
      std::unique_ptr<Process::GenericProcess> process_;
      std::unique_ptr<Hadroniser::GenericHadroniser> hadroniser_;
      bool store_;
  };
}

#endif
