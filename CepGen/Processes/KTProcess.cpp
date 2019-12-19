#include "CepGen/Processes/KTProcess.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen
{
  namespace proc
  {
    KTProcess::KTProcess( const ParametersList& params,
                          const std::string& name,
                          const std::string& description,
                          const std::array<pdgid_t,2>& partons,
                          const std::vector<pdgid_t>& central ) :
      Process( params, name, description+" (kT-factor.)" ),
      qt1_( 0. ), phi_qt1_( 0. ), qt2_( 0. ), phi_qt2_( 0. ),
      kIntermediateParts( partons ), kProducedParts( central )
    {}

    void
    KTProcess::addEventContent()
    {
      Process::setEventContent(
        { // incoming state
          { Particle::IncomingBeam1, PDG::proton },
          { Particle::IncomingBeam2, PDG::proton },
          { Particle::Parton1, kIntermediateParts[0] },
          { Particle::Parton2, kIntermediateParts[1] }
        },
        { // outgoing state
          { Particle::OutgoingBeam1, { PDG::proton } },
          { Particle::OutgoingBeam2, { PDG::proton } },
          { Particle::CentralSystem, kProducedParts }
        }
      );
      setExtraContent();
    }

    void
    KTProcess::prepareKinematics()
    {
      const KTFlux flux1 = (KTFlux)kin_.incoming_beams.first.kt_flux,
                   flux2 = (KTFlux)kin_.incoming_beams.second.kt_flux;

      if ( kin_.mode == KinematicsMode::invalid ) {
        //--- try to extrapolate kinematics mode from unintegrated fluxes
        bool el1 = ( flux1 == KTFlux::P_Photon_Elastic
                  || flux1 == KTFlux::HI_Photon_Elastic
                  || flux1 == KTFlux::P_Gluon_KMR );
        bool el2 = ( flux2 == KTFlux::P_Photon_Elastic
                  || flux2 == KTFlux::HI_Photon_Elastic
                  || flux2 == KTFlux::P_Gluon_KMR );
        if ( el1 && el2 )
          kin_.mode = KinematicsMode::ElasticElastic;
        else if ( el1 )
          kin_.mode = KinematicsMode::ElasticInelastic;
        else if ( el2 )
          kin_.mode = KinematicsMode::InelasticElastic;
        else
          kin_.mode = KinematicsMode::InelasticInelastic;
      }
      else {
        //--- try to extrapolate unintegrated fluxes from kinematics mode
        const HeavyIon hi1( kin_.incoming_beams.first.pdg ), hi2( kin_.incoming_beams.second.pdg );
        //==========================================================================================
        // ensure the first incoming flux is compatible with the kinematics mode
        //==========================================================================================
        if ( ( kin_.mode == KinematicsMode::ElasticElastic
            || kin_.mode == KinematicsMode::ElasticInelastic )
          && ( flux1 != KTFlux::P_Photon_Elastic ) ) {
          kin_.incoming_beams.first.kt_flux = hi1
            ? KTFlux::HI_Photon_Elastic
            : KTFlux::P_Photon_Elastic;
          CG_DEBUG( "KTProcess:kinematics" )
            << "KT flux for first incoming parton set to \""
            << kin_.incoming_beams.first.kt_flux << "\".";
        }
        else if ( flux1 != KTFlux::P_Photon_Inelastic
               && flux1 != KTFlux::P_Photon_Inelastic_Budnev ) {
          if ( hi1 )
            throw CG_FATAL( "KTProcess:kinematics" )
              << "Inelastic photon emission from HI not yet supported!";
          kin_.incoming_beams.first.kt_flux = KTFlux::P_Photon_Inelastic_Budnev;
          CG_INFO( "KTProcess:kinematics" )
            << "KT flux for first incoming parton set to \""
            << kin_.incoming_beams.first.kt_flux << "\".";
        }
        //==========================================================================================
        // ensure the second incoming flux is compatible with the kinematics mode
        //==========================================================================================
        if ( ( kin_.mode == KinematicsMode::ElasticElastic
            || kin_.mode == KinematicsMode::InelasticElastic )
          && ( flux2 != KTFlux::P_Photon_Elastic ) ) {
          kin_.incoming_beams.second.kt_flux = hi2
            ? KTFlux::HI_Photon_Elastic
            : KTFlux::P_Photon_Elastic;
          CG_DEBUG( "KTProcess:kinematics" )
            << "KT flux for second incoming parton set to \""
            << kin_.incoming_beams.second.kt_flux << "\".";
        }
        else if ( flux2 != KTFlux::P_Photon_Inelastic
               && flux2 != KTFlux::P_Photon_Inelastic_Budnev ) {
          if ( hi2 )
            throw CG_FATAL( "KTProcess:kinematics" )
              << "Inelastic photon emission from HI not yet supported!";
          kin_.incoming_beams.second.kt_flux = KTFlux::P_Photon_Inelastic_Budnev;
          CG_INFO( "KTProcess:kinematics" )
            << "KT flux for second incoming parton set to \""
            << kin_.incoming_beams.second.kt_flux << "\".";
        }
      }

      //============================================================================================
      // register the incoming partons' variables
      //============================================================================================

      defineVariable( qt1_, Mapping::exponential, kin_.cuts.initial.qt, { 1.e-10, 500. }, "First incoming parton virtuality" );
      defineVariable( qt2_, Mapping::exponential, kin_.cuts.initial.qt, { 1.e-10, 500. }, "Second incoming parton virtuality" );
      defineVariable( phi_qt1_, Mapping::linear, kin_.cuts.initial.phi_qt, { 0., 2.*M_PI }, "First incoming parton azimuthal angle" );
      defineVariable( phi_qt2_, Mapping::linear, kin_.cuts.initial.phi_qt, { 0., 2.*M_PI }, "Second incoming parton azimuthal angle" );

      //============================================================================================
      // register the incoming partons
      //============================================================================================

      switch ( kin_.incoming_beams.first.kt_flux ) {
        case KTFlux::P_Gluon_KMR:
          event_->oneWithRole( Particle::Parton1 ).setPdgId( PDG::gluon ); break;
        case KTFlux::P_Photon_Elastic:
        case KTFlux::P_Photon_Inelastic:
        case KTFlux::P_Photon_Inelastic_Budnev:
        case KTFlux::HI_Photon_Elastic:
          event_->oneWithRole( Particle::Parton1 ).setPdgId( PDG::photon ); break;
        case KTFlux::invalid: default:
          throw CG_FATAL( "KTProcess:kinematics" )
            << "Invalid flux for 2nd incoming parton: " << kin_.incoming_beams.first.kt_flux << "!";
      }
      switch ( kin_.incoming_beams.second.kt_flux ) {
        case KTFlux::P_Gluon_KMR:
          event_->oneWithRole( Particle::Parton2 ).setPdgId( PDG::gluon ); break;
        case KTFlux::P_Photon_Elastic:
        case KTFlux::P_Photon_Inelastic:
        case KTFlux::P_Photon_Inelastic_Budnev:
        case KTFlux::HI_Photon_Elastic:
          event_->oneWithRole( Particle::Parton2 ).setPdgId( PDG::photon ); break;
        case KTFlux::invalid: default:
          throw CG_FATAL( "KTProcess:kinematics" )
            << "Invalid flux for 2nd incoming parton: " << kin_.incoming_beams.second.kt_flux << "!";
      }

      //============================================================================================
      // register all process-dependent variables
      //============================================================================================

      preparePhaseSpace();

      //============================================================================================
      // register the outgoing remnants' variables
      //============================================================================================

      mX2_ = event_->oneWithRole( Particle::IncomingBeam1 ).mass2();
      mY2_ = event_->oneWithRole( Particle::IncomingBeam2 ).mass2();
      if ( kin_.mode == KinematicsMode::InelasticElastic || kin_.mode == KinematicsMode::InelasticInelastic )
        defineVariable( mX2_, Mapping::square, kin_.cuts.remnants.mass_single, { 1.07, 1000. }, "Positive z proton remnant squared mass" );
      if ( kin_.mode == KinematicsMode::ElasticInelastic || kin_.mode == KinematicsMode::InelasticInelastic )
        defineVariable( mY2_, Mapping::square, kin_.cuts.remnants.mass_single, { 1.07, 1000. }, "Negative z proton remnant squared mass" );
    }

    double
    KTProcess::computeWeight()
    {
      return std::max( 0., computeKTFactorisedMatrixElement() );
    }

    void
    KTProcess::fillKinematics( bool )
    {
      fillCentralParticlesKinematics(); // process-dependent!
      fillPrimaryParticlesKinematics();
    }

    void
    KTProcess::fillPrimaryParticlesKinematics()
    {
      //============================================================================================
      //     outgoing protons
      //============================================================================================

      Particle& op1 = event_->oneWithRole( Particle::OutgoingBeam1 ),
               &op2 = event_->oneWithRole( Particle::OutgoingBeam2 );

      op1.setMomentum( pX_ );
      op2.setMomentum( pY_ );

      switch ( kin_.mode ) {
        case KinematicsMode::ElasticElastic:
          op1.setStatus( Particle::Status::FinalState );
          op2.setStatus( Particle::Status::FinalState );
          break;
        case KinematicsMode::ElasticInelastic:
          op1.setStatus( Particle::Status::FinalState );
          op2.setStatus( Particle::Status::Unfragmented ).setMass( sqrt( mY2_ ) );
          break;
        case KinematicsMode::InelasticElastic:
          op1.setStatus( Particle::Status::Unfragmented ).setMass( sqrt( mX2_ ) );
          op2.setStatus( Particle::Status::FinalState );
          break;
        case KinematicsMode::InelasticInelastic:
          op1.setStatus( Particle::Status::Unfragmented ).setMass( sqrt( mX2_ ) );
          op2.setStatus( Particle::Status::Unfragmented ).setMass( sqrt( mY2_ ) );
          break;
        default: {
          throw CG_FATAL( "KTProcess" )
            << "This kT factorisation process is intended for p-on-p collisions! Aborting.";
        } break;
      }

      //============================================================================================
      //     incoming partons (photons, pomerons, ...)
      //============================================================================================

      Particle& g1 = event_->oneWithRole( Particle::Parton1 );
      g1.setMomentum( event_->oneWithRole( Particle::IncomingBeam1 ).momentum()-pX_, true );

      Particle& g2 = event_->oneWithRole( Particle::Parton2 );
      g2.setMomentum( event_->oneWithRole( Particle::IncomingBeam2 ).momentum()-pY_, true );

      //============================================================================================
      //     two-parton system
      //============================================================================================

      event_->oneWithRole( Particle::Intermediate ).setMomentum( g1.momentum()+g2.momentum() );
    }
  }
}
