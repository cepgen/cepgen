#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Hadronisers/HadronisersHandler.h"

#include "CepGen/IO/PythiaEventInterface.h"

#include "CepGen/Parameters.h"
#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

#include <Pythia8/Pythia.h>

#include <memory>
#include <unordered_map>
#include <vector>

namespace cepgen
{
  namespace hadr
  {
    /**
     * Full interface to the Pythia8 hadronisation algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     * \brief Pythia8 hadronisation algorithm
     */
    class Pythia8Hadroniser : public GenericHadroniser
    {
      public:
        explicit Pythia8Hadroniser( const ParametersList& );
        ~Pythia8Hadroniser();

        void setParameters( const Parameters& ) override;
        void readString( const char* param ) override;
        void init() override;
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override;

      private:
        static constexpr unsigned short PYTHIA_STATUS_IN_BEAM = 12;
        static constexpr unsigned short PYTHIA_STATUS_IN_PARTON_KT = 61;

        std::vector<unsigned int> min_ids_;
        std::unordered_map<short,short> py_cg_corresp_;
        unsigned short findRole( const Event& ev, const Pythia8::Particle& p ) const;
        void updateEvent( Event& ev, double& weight ) const;
        Particle& addParticle( Event& ev, const Pythia8::Particle&, const Pythia8::Vec4& mom, unsigned short ) const;
        /// A Pythia8 core to be wrapped
        std::unique_ptr<Pythia8::Pythia> pythia_;
        std::unique_ptr<Pythia8::CepGenEvent> cg_evt_;
        const bool correct_central_;
        const bool debug_lhef_;
        bool res_decay_;
        bool enable_hadr_;
        unsigned short offset_;
        bool first_evt_;
    };

    Pythia8Hadroniser::Pythia8Hadroniser( const ParametersList& plist ) :
      GenericHadroniser( plist, "pythia8" ),
      pythia_( new Pythia8::Pythia ), cg_evt_( new Pythia8::CepGenEvent ),
      correct_central_( plist.get<bool>( "correctCentralSystem", false ) ),
      debug_lhef_( plist.get<bool>( "debugLHEF", false ) ),
      res_decay_( true ), enable_hadr_( false ),
      offset_( 0 ), first_evt_( true )
    {}

    void
    Pythia8Hadroniser::setParameters( const Parameters& params )
    {
      params_ = &params;
      cg_evt_->initialise( params );
      pythia_->setLHAupPtr( (Pythia8::LHAup*)cg_evt_.get() );
      pythia_->settings.parm( "Beams:idA", (short)params_->kinematics.incoming_beams.first.pdg );
      pythia_->settings.parm( "Beams:idB", (short)params_->kinematics.incoming_beams.second.pdg );
      // specify we will be using a LHA input
      pythia_->settings.mode( "Beams:frameType", 5 );
      pythia_->settings.parm( "Beams:eCM", params_->kinematics.sqrtS() );
      min_ids_ = params_->kinematics.minimum_final_state;
      if ( debug_lhef_ )
        cg_evt_->openLHEF( "debug.lhe" );
    }

    Pythia8Hadroniser::~Pythia8Hadroniser()
    {
      pythia_->settings.writeFile( "last_pythia_config.cmd", false );
      if ( debug_lhef_ )
        cg_evt_->closeLHEF( true );
    }

    void
    Pythia8Hadroniser::readString( const char* param )
    {
      if ( !pythia_->readString( param ) )
        throw CG_FATAL( "Pythia8Hadroniser" ) << "The Pythia8 core failed to parse the following setting:\n\t" << param;
    }

    void
    Pythia8Hadroniser::init()
    {
      pythia_->settings.flag( "ProcessLevel:resonanceDecays", res_decay_ );
      if ( pythia_->settings.flag( "ProcessLevel:all" ) != enable_hadr_ )
        pythia_->settings.flag( "ProcessLevel:all", enable_hadr_ );

      if ( seed_ == -1ll )
        pythia_->settings.flag( "Random:setSeed", false );
      else {
        pythia_->settings.flag( "Random:setSeed", true );
        pythia_->settings.mode( "Random:seed", seed_ );
      }

#if defined( PYTHIA_VERSION_INTEGER ) && PYTHIA_VERSION_INTEGER >= 8226
      switch ( params_->process()->mode() ) {
        case KinematicsMode::ElasticElastic: {
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 3 );
          pythia_->settings.flag( "PartonLevel:all", false );
        } break;
        case KinematicsMode::InelasticElastic: {
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 2 );
        } break;
        case KinematicsMode::ElasticInelastic: {
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 1 );
        } break;
        case KinematicsMode::InelasticInelastic: default: {
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 0 );
        } break;
      }
#else
      CG_WARNING( "Pythia8Hadroniser" )
        << "Beam remnants framework for this version of Pythia "
        << "(" << Form( "%.3f", PYTHIA_VERSION ) << ")\n\t"
        << "does not support mixing of unresolved hadron states.\n\t"
        << "The proton remnants output might hence be wrong.\n\t"
        << "Please update the Pythia version or disable this part.";
#endif
      if ( correct_central_ && res_decay_ )
        CG_WARNING( "Pythia8Hadroniser" )
          << "Central system's kinematics correction enabled while resonances are\n\t"
          << "expected to be decayed. Please check that this is fully intended.";

      if ( !pythia_->init() )
        throw CG_FATAL( "Pythia8Hadroniser" )
          << "Failed to initialise the Pythia8 core!\n\t"
          << "See the message above for more details.";

      if ( debug_lhef_ )
        cg_evt_->initLHEF();
    }

    void
    Pythia8Hadroniser::setCrossSection( double xsec, double xsec_err )
    {
      cg_evt_->setCrossSection( 0, xsec, xsec_err );
    }

    bool
    Pythia8Hadroniser::run( Event& ev, double& weight, bool full )
    {
      //--- initialise the event weight before running any decay algorithm
      weight = 1.;

      //--- only launch Pythia if:
      // 1) the full event kinematics (i.e. with remnants) is to be specified,
      // 2) the remnants are to be fragmented, or
      // 3) the resonances are to be decayed.
      if ( full && !remn_fragm_ && !res_decay_ )
        return true;
      if ( !full && !res_decay_ )
        return true;

      //--- switch full <-> partial event
      if ( full != enable_hadr_ ) {
        enable_hadr_ = full;
        init();
      }

      //===========================================================================================
      // convert our event into a custom LHA format
      //===========================================================================================

      cg_evt_->feedEvent( ev, full ? Pythia8::CepGenEvent::Type::centralAndBeamRemnants
                                   : Pythia8::CepGenEvent::Type::centralAndPartons );
      //if ( full ) cg_evt_->listEvent();
      if ( debug_lhef_ && full )
        cg_evt_->eventLHEF();

      //===========================================================================================
      // launch the hadronisation / resonances decays, and update the event accordingly
      //===========================================================================================

      ev.num_hadronisation_trials = 0;
      while ( true ) {
        if ( ev.num_hadronisation_trials++ > max_trials_ )
          return false;
        //--- run the hadronisation/fragmentation algorithm
        if ( pythia_->next() ) {
          //--- hadronisation successful
          if ( first_evt_ && full ) {
            offset_ = 0;
            for ( unsigned short i = 1; i < pythia_->event.size(); ++i )
              if ( pythia_->event[i].status() == -PYTHIA_STATUS_IN_BEAM )
                //--- no incoming particles in further stages
                offset_++;
            first_evt_ = false;
          }
          break;
        }
      }
      CG_DEBUG( "Pythia8Hadroniser" )
        << "Pythia8 hadronisation performed successfully.\n\t"
        << "Number of trials: " << ev.num_hadronisation_trials << "/" << max_trials_ << ".\n\t"
        << "Particles multiplicity: " << ev.particles().size() << " → " << pythia_->event.size() << ".\n\t"
        << "  indices offset: " << offset_ << ".";

      //===========================================================================================
      // update the event content with Pythia's output
      //===========================================================================================

      updateEvent( ev, weight );
      return true;
    }

    Particle&
    Pythia8Hadroniser::addParticle( Event& ev, const Pythia8::Particle& py_part, const Pythia8::Vec4& mom, unsigned short role ) const
    {
      Particle& op = ev.addParticle( (Particle::Role)role );
      ParticleProperties prop;
      const pdgid_t pdg_id = abs( py_part.id() );
      try { prop = PDG::get()( pdg_id ); } catch ( const Exception& ) {
        prop = ParticleProperties{ pdg_id,
          py_part.name(), py_part.name(),
          (short)py_part.col(), // colour factor
          py_part.m0(), py_part.mWidth(),
          (short)py_part.charge(), // charge
          py_part.isLepton()
        };
        PDG::get().define( prop );
      }
      op.setPdgId( (short)py_part.id() );
      op.setStatus( py_part.isFinal()
        ? Particle::Status::FinalState
        : (Particle::Role)role == Particle::Role::CentralSystem
          ? Particle::Status::Propagator
          : Particle::Status::Fragmented );
      op.setMomentum( Particle::Momentum( mom.px(), mom.py(), mom.pz(), mom.e() ) );
      op.setMass( mom.mCalc() );
      cg_evt_->addCorresp( py_part.index()-offset_, op.id() );
      return op;
    }

    void
    Pythia8Hadroniser::updateEvent( Event& ev, double& weight ) const
    {
      std::vector<unsigned short> central_parts;

      for ( unsigned short i = 1+offset_; i < pythia_->event.size(); ++i ) {
        const Pythia8::Particle& p = pythia_->event[i];
        const unsigned short cg_id = cg_evt_->cepgenId( i-offset_ );
        if ( cg_id != Pythia8::CepGenEvent::INVALID_ID ) {
          //----- particle already in the event
          Particle& cg_part = ev[cg_id];
          //--- fragmentation result
          if ( cg_part.role() == Particle::OutgoingBeam1
            || cg_part.role() == Particle::OutgoingBeam2 ) {
            cg_part.setStatus( Particle::Status::Fragmented );
            continue;
          }
          //--- resonance decayed; apply branching ratio for this decay
          if ( cg_part.role() == Particle::CentralSystem && p.status() < 0 ) {
            if ( res_decay_ )
              weight *= p.particleDataEntry().pickChannel().bRatio();
            cg_part.setStatus( Particle::Status::Resonance );
            central_parts.emplace_back( i );
          }
          //--- particle is not what we expect
          if ( p.idAbs() != abs( cg_part.integerPdgId() ) ) {
            CG_INFO( "Pythia8Hadroniser:update" ) << "LHAEVT event content:";
            cg_evt_->listEvent();
            CG_INFO( "Pythia8Hadroniser:update" ) << "Pythia event content:";
            pythia_->event.list();
            CG_INFO( "Pythia8Hadroniser:update" ) << "CepGen event content:";
            ev.dump();
            CG_INFO( "Pythia8Hadroniser:update" ) << "Correspondence:";
            cg_evt_->dumpCorresp();

            throw CG_FATAL( "Pythia8Hadroniser:update" )
              << "Event list corruption detected for (Pythia/CepGen) particle " << i << "/" << cg_id << ":\n\t"
              << "should be " << abs( p.id() ) << ", "
              << "got " << cg_part.integerPdgId() << "!";
          }
        }
        //--- check for messed up particles parentage and discard incoming beam particles
        /*else if ( p.mother1() > i || p.mother1() <= offset_ )
          continue;
        else if ( p.mother2() > i || p.mother2() <= offset_ )
          continue;*/
        else {
          //----- new particle to be added
          const unsigned short role = findRole( ev, p );
          switch ( (Particle::Role)role ) {
            case Particle::OutgoingBeam1:
              ev[Particle::OutgoingBeam1][0].setStatus( Particle::Status::Fragmented );
              break;
            case Particle::OutgoingBeam2:
              ev[Particle::OutgoingBeam2][0].setStatus( Particle::Status::Fragmented );
              break;
            default: break;
          }
          // found the role ; now we can add the particle
          Particle& cg_part = addParticle( ev, p, p.p(), role );
          if ( correct_central_ && (Particle::Role)role == Particle::CentralSystem ) {
            const auto& ip = std::find( central_parts.begin(), central_parts.end(), p.mother1() );
            if ( ip != central_parts.end() )
              cg_part.setMomentum( ev[cg_evt_->cepgenId( *ip-offset_ )].momentum() );
          }
          for ( const auto& moth_id : p.motherList() ) {
            if ( moth_id <= offset_ )
              continue;
            const unsigned short moth_cg_id = cg_evt_->cepgenId( moth_id-offset_ );
            if ( moth_cg_id != Pythia8::CepGenEvent::INVALID_ID )
              cg_part.addMother( ev[moth_cg_id] );
            else
              cg_part.addMother( addParticle( ev, pythia_->event[moth_id], p.p(), role ) );
            if ( !p.isFinal() ) {
              if ( p.isResonance() || !p.daughterList().empty() )
                cg_part.setStatus( Particle::Status::Resonance );
              else
                cg_part.setStatus( Particle::Status::Undefined );
            }
          }
        }
      }
    }

    unsigned short
    Pythia8Hadroniser::findRole( const Event& ev, const Pythia8::Particle& p ) const
    {
      for ( const auto& par_id : p.motherList() ) {
        if ( par_id == 1 && offset_ > 0 )
          return (unsigned short)Particle::OutgoingBeam1;
        if ( par_id == 2 && offset_ > 0 )
          return (unsigned short)Particle::OutgoingBeam2;
        const unsigned short par_cg_id = cg_evt_->cepgenId( par_id-offset_ );
        if ( par_cg_id != Pythia8::CepGenEvent::INVALID_ID )
          return (unsigned short)ev[par_cg_id].role();
        return findRole( ev, pythia_->event[par_id] );
      }
      return (unsigned short)Particle::UnknownRole;
    }
  }
}
// register hadroniser and define alias
REGISTER_HADRONISER( pythia8, Pythia8Hadroniser )
