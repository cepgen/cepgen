#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Parameters.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"

#ifdef LIBHEPMC
#  include "HepMC/Version.h"
#  ifndef HEPMC_VERSION_CODE // HepMC v2
#    include "HepMC/IO_GenEvent.h"
#    include "HepMC/SimpleVector.h"
#  else // HepMC v3+
#    define HEPMC_VERSION3
#    include "HepMC/WriterAscii.h"
#    include "HepMC/FourVector.h"
#  endif
#  include "HepMC/GenEvent.h"
#  include "HepMC/GenVertex.h"
#  include "HepMC/GenParticle.h"
#endif

#include <memory>

namespace cepgen
{
  namespace output
  {
    /**
     * \brief Handler for the HepMC file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class HepMCHandler : public GenericExportHandler
    {
      public:
        /// Class constructor
        /// \param[in] filename Output file path
        /// \param[in] type Output type
        HepMCHandler( const char* filename, const GenericExportHandler::OutputType& type = GenericExportHandler::HepMC );
        /// Class constructor
        HepMCHandler( const ParametersList& );
        void initialise( const Parameters& /*params*/ ) override {}
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override;

      private:
        /// Clear the associated HepMC event content
        void clearEvent();
        /// Populate the associated HepMC event with a Event object
        void fillEvent( const Event& );
#ifdef LIBHEPMC
#  ifdef HEPMC_VERSION3
        /// Writer object (from HepMC v3+)
        std::unique_ptr<HepMC::WriterAscii> output_;
        /// Generator cross section and error
        HepMC::GenCrossSectionPtr xs_;
#  else
        /// Writer object (from HepMC v<3)
        std::unique_ptr<HepMC::IO_GenEvent> output_;
        /// Generator cross section and error
        HepMC::GenCrossSection xs_;
#  endif
        /// Associated HepMC event
        std::shared_ptr<HepMC::GenEvent> event_;
#endif
    };

    HepMCHandler::HepMCHandler( const char* filename, const GenericExportHandler::OutputType& type ) :
      GenericExportHandler( type )
#ifdef LIBHEPMC
#  ifdef HEPMC_VERSION3
      , output_( new HepMC::WriterAscii( filename ) ),
#  else
      , output_( new HepMC::IO_GenEvent( filename ) ),
#  endif
      event_( new HepMC::GenEvent() )
#endif
    {}

    HepMCHandler::HepMCHandler( const ParametersList& params ) :
      GenericExportHandler( (GenericExportHandler::OutputType)params.get<int>( "type" ) )
#ifdef LIBHEPMC
#  ifdef HEPMC_VERSION3
      , output_( new HepMC::WriterAscii( params.get<std::string>( "filename" ) ) ),
#  else
      , output_( new HepMC::IO_GenEvent( params.get<std::string>( "filename" ) ) ),
#  endif
      event_( new HepMC::GenEvent() )
#endif
    {}

    void
    HepMCHandler::operator<<( const Event& evt )
    {
      fillEvent( evt );
#ifdef LIBHEPMC
#  ifdef HEPMC_VERSION3
      output_->write_event( *event_ );
#  else
      output_->write_event( event_.get() );
#  endif
#endif
    }

    void
    HepMCHandler::setCrossSection( double xsect, double xsect_err )
    {
#ifdef LIBHEPMC
#  ifdef HEPMC_VERSION3
      xs_->set_cross_section( xsect, xsect_err );
      event_->add_attribute( "AlphaQCD", HepMC::make_shared<HepMC::DoubleAttribute>( constants::ALPHA_QCD ) );
      event_->add_attribute( "AlphaEM", HepMC::make_shared<HepMC::DoubleAttribute>( constants::ALPHA_EM ) );
#  else
      xs_.set_cross_section( xsect, xsect_err );
      event_->set_alphaQCD( constants::ALPHA_QCD );
      event_->set_alphaQED( constants::ALPHA_EM );
#  endif
#endif
    }

    void
    HepMCHandler::fillEvent( const Event& evt )
    {
#ifdef LIBHEPMC
      event_->clear();

      // general information
      event_->set_cross_section( xs_ );

      event_->set_event_number( event_num_ );
      event_->weights().push_back( 1. ); // unweighted events

      // filling the particles content
      const HepMC::FourVector origin( 0., 0., 0., 0. );
      Particles part_vec = evt.particles();

      int cm_id = 0, idx = 1;

#  ifdef HEPMC_VERSION3
      HepMC::GenVertexPtr v1 = HepMC::make_shared<HepMC::GenVertex>( origin ),
                          v2 = HepMC::make_shared<HepMC::GenVertex>( origin ),
                          vcm = HepMC::make_shared<HepMC::GenVertex>( origin );
#  else
      HepMC::GenVertex* v1 = new HepMC::GenVertex( origin ),
                       *v2 = new HepMC::GenVertex( origin ),
                       *vcm = new HepMC::GenVertex( origin );
#  endif
      for ( unsigned int i = 0; i < part_vec.size(); ++i ) {
        const Particle part_orig = part_vec.at( i );
        HepMC::FourVector pmom( part_orig.momentum().px(),
                                part_orig.momentum().py(),
                                part_orig.momentum().pz(),
                                part_orig.energy() );
#  ifdef HEPMC_VERSION3
        HepMC::GenParticlePtr part = HepMC::make_shared<HepMC::GenParticle>( pmom, part_orig.integerPdgId(), (int)part_orig.status() );
#  else
        HepMC::GenParticle* part = new HepMC::GenParticle( pmom, part_orig.integerPdgId(), (int)part_orig.status() );
        part->suggest_barcode( idx++ );
#  endif
        const ParticlesIds moth = part_orig.mothers();

        switch ( part_orig.role() ) {
          case Particle::IncomingBeam1: { v1->add_particle_in( part ); } break;
          case Particle::IncomingBeam2: { v2->add_particle_in( part ); } break;
          case Particle::OutgoingBeam1: { v1->add_particle_out( part ); } break;
          case Particle::OutgoingBeam2: { v2->add_particle_out( part ); } break;
          case Particle::Parton1:       { v1->add_particle_out( part ); vcm->add_particle_in( part ); } break;
          case Particle::Parton2:       { v2->add_particle_out( part ); vcm->add_particle_in( part ); } break;
          case Particle::Intermediate:  { cm_id = i; continue; } break;
          case Particle::CentralSystem:
          default: {
            if ( moth.empty() ) continue;
            if ( *moth.begin() == cm_id ) vcm->add_particle_out( part );
            else
              throw CG_FATAL( "HepMCHandler:fillEvent" )
                << "Other particle requested! Not yet implemented!";
          } break;
        }
        idx++;
      }
      event_->add_vertex( v1 );
      event_->add_vertex( v2 );
      event_->add_vertex( vcm );

#  ifndef HEPMC_VERSION3
      event_->set_beam_particles( *v1->particles_in_const_begin(), *v2->particles_in_const_begin() );
      event_->set_signal_process_vertex( *v1->vertices_begin() );
      event_->set_beam_particles( *v1->particles_in_const_begin(), *v2->particles_in_const_end() );
#  endif
#endif
      event_num_++;
    }
  }
}

REGISTER_IO_MODULE( hepmc, HepMCHandler )
