#ifdef PYTHIA8
#include "Pythia8Hadroniser.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Version.h"

namespace CepGen
{
  namespace Hadroniser
  {
    Pythia8Hadroniser::Pythia8Hadroniser() :
      GenericHadroniser( "pythia8" )
    {
#ifdef PYTHIA8
      pythia_.reset( new Pythia8::Pythia );

      //--- start by disabling some unnecessary output
      pythia_->readString( "Next:numberCount = 0" );
#endif
    }

    Pythia8Hadroniser::~Pythia8Hadroniser()
    {}

    bool
    Pythia8Hadroniser::init()
    {
      bool res = pythia_->init();
      if ( !res ) {
        FatalError( "Failed to initialise the Pythia8 core!\n\t"
                    "See the message above for more details." );
      }
      return res;
    }

    void
    Pythia8Hadroniser::readString( const char* param )
    {
      if ( !pythia_->readString( param ) ) {
        FatalError( Form( "The Pythia8 core failed to parse the following setting:\n\t%s", param ) );
      }
    }

    void
    Pythia8Hadroniser::setSeed( long long seed )
    {
#ifdef PYTHIA8
      if ( seed == -1ll ) {
        pythia_->settings.flag( "Random:setSeed", false );
        return;
      }
      pythia_->settings.flag( "Random:setSeed", true );
      pythia_->settings.parm( "Random:seed", seed );
#endif
    }

    bool
    Pythia8Hadroniser::hadronise( Event& ev, double& weight )
    {
      weight = 1.;
#ifdef PYTHIA8
      //--- start by cleaning up the previous runs leftovers

      const double sqrt_s = ev.cmEnergy();

      pythia_->event.reset();
      pythia_->event[0].e( sqrt_s );
      pythia_->event[0].m( sqrt_s );

      const unsigned short num_before = ev.numParticles();
      std::map<short,short> py_cg_corresp, cg_py_corresp;

      //--- loop to add the particles

      for ( unsigned short i = 0; i < num_before; ++i ) {
        const Particle& part = ev.getConstById( i );

        const Particle::Momentum mom = part.momentum();
        Pythia8::Particle py8part( part.integerPdgId(), 0, 0, 0, 0, 0, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), part.mass() );
        switch ( part.role() ) {
          case Particle::IncomingBeam1:
          case Particle::IncomingBeam2:
          case Particle::Parton1:
          case Particle::Parton2:
          case Particle::Parton3:
          case Particle::Intermediate:
          case Particle::UnknownRole:
            continue;
          case Particle::CentralSystem: {
            if ( !pythia_->particleData.canDecay( py8part.id() ) ) continue;
            if ( !pythia_->particleData.mayDecay( py8part.id() ) ) continue;
            py8part.status( 93 );
          } break;
          case Particle::OutgoingBeam1:
          case Particle::OutgoingBeam2: {
            //if ( part.status() != Particle::Unfragmented ) continue;
            if ( part.status() == Particle::Unfragmented ) {
              //py8part.id( 9902210 ); // diffractive proton
              py8part.id( 2212 );
              //py8part.status( 15 );
              //py8part.status( 12 );
              //py8part.status( 2 );
              py8part.status( 93 );
              //py8part.status( 63 );
            }
            else {
              py8part.status( 23 );
            }
          } break;
        }
        const unsigned short py_id = pythia_->event.append( py8part );
        cg_py_corresp[part.id()] = py_id;
        py_cg_corresp[py_id] = part.id();
      }

      //pythia_->event.scale( ev.getOneByRole( Particle::Intermediate ).mass() );
      //pythia_->event.list(true,true);
      //pythia_->event.listJunctions();
      const unsigned short num_py_parts = pythia_->event.size();

      if ( !pythia_->next() ) {
        InWarning( "Pythia8 failed to process the event." );
        return false;
      }

//pythia_->event.list(true,true);
//exit(0);

      // check if something happened in the event processing by Pythia
      // if not, return the event as it is...
      if ( pythia_->event.size() == num_py_parts ) return true;

      for ( unsigned short i = 1; i < pythia_->event.size(); ++i ) { // skip the central system
        const Pythia8::Particle p = pythia_->event[i];
        std::map<short,short>::const_iterator it = py_cg_corresp.find( i );
        if ( it != py_cg_corresp.end() ) {
          if ( p.daughterList().size() > 0 ) {
            weight *= p.particleDataEntry().pickChannel().bRatio();
            ev.getById( it->second ).setStatus( Particle::Resonance );
          }
          continue;
        }

        const std::vector<int> mothers = p.motherList();
        if ( mothers.size() == 0 ) continue;
        Particle& op = ev.addParticle( Particle::CentralSystem );
        cg_py_corresp[op.id()] = i;
        py_cg_corresp[i] = op.id();

        op.setPdgId( static_cast<ParticleCode>( abs( p.id() ) ), p.charge() );
        if ( p.isFinal() ) op.setStatus( Particle::FinalState );
        else op.setStatus( Particle::Propagator );

        op.setMomentum( Particle::Momentum( p.px(), p.py(), p.pz(), p.e() ) );
        for ( std::vector<int>::const_iterator moth = mothers.begin(); moth != mothers.end(); ++moth ) {
          if ( *moth != 0 && py_cg_corresp.count( *moth ) == 0 ) {
            FatalError( Form( "Particle with id=%d was not found in the event content!", *moth ) );
          }
          op.addMother( ev.getById( py_cg_corresp[*moth] ) );
        }
      }
      //ev.dump();
#else
      FatalError( "Pythia8 is not linked to this instance!" );
#endif
      return true;
    }
  }
}

#endif
