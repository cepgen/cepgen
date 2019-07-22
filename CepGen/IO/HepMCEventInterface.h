#ifndef CepGen_IO_HepMCEventInterface_h
#define CepGen_IO_HepMCEventInterface_h

namespace cepgen {
  class Parameters;
  class Event;
  class Particle;
}

#ifdef HEPMC3
#  include "HepMC3/Version.h"
#  include "HepMC3/GenEvent.h"
#  define HepMC HepMC3
#else
#  include "HepMC/Version.h"
#  if !defined( HEPMC_VERSION_CODE ) // HepMC v2
#    include "HepMC/GenEvent.h"
#  else
#    include "HepMC/GenEvent.h"
#    define HEPMC3
#  endif
#endif

namespace HepMC
{
  /// Interfacing between CepGen and HepMC event definitions
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2019
  class CepGenEvent : public GenEvent
  {
    public:
      explicit CepGenEvent();
      /// Initialise this conversion object with CepGen parameters
      void initialise( const cepgen::Parameters& );
      /// Feed a new CepGen event to this conversion object
      /// \param[in] ev CepGen event to be fed
      void feedEvent( const cepgen::Event& ev );
      /// Set the cross section for a given process
      /// \param[in] id Process identifier
      /// \param[in] xsec Process cross section, in pb
      /// \param[in] xsec_err Uncertainty on process cross section, in pb
      void setCrossSection( int id, double xsec, double xsec_err );
      /// Specify new process attributes
      /// \param[in] id Process identifier
      /// \param[in] xsec Process cross section, in pb
      /// \param[in] q2_scale Hard event scale \f$Q^2\f$, in GeV\f$^2\f$
      /// \param[in] alpha_qed \f$\alpha_{\rm em}\f$ for this process
      /// \param[in] alpha_qcd \f$\alpha_{\rm s}\f$ for this process
      void setProcess( int id, double xsec, double q2_scale, double alpha_qed, double alpha_qcd );

    private:
      const cepgen::Parameters* params_; // borrowed
  };
}
#endif

