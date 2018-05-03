#include "CepGen/Processes/GenericKTProcess.h"

#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PDG.h"

#ifdef KMR_FLUX
extern "C" {
  extern void f_inter_kmr_fg_( double& rx, double& rkt2, double& rmu2, int& iread, double& parton);
}
#endif

namespace CepGen
{
  namespace Process
  {
    GenericKTProcess::GenericKTProcess( const std::string& name,
                                        const std::string& description,
                                        const std::array<PDG,2>& partons,
                                        const std::vector<PDG>& central ) :
      GenericProcess( name, description+" (kT-factorisation approach)" ),
      num_dimensions_( 0 ), kt_jacobian_( 0. ),
      qt1_( 0. ), phi_qt1_( 0. ),
      qt2_( 0. ), phi_qt2_( 0. ),
      flux1_( 0. ), flux2_( 0. ),
      kIntermediateParts( partons ), kProducedParts( central )
    {}

    void
    GenericKTProcess::addEventContent()
    {
      GenericProcess::setEventContent(
        { // incoming state
          { Particle::IncomingBeam1, PDG::Proton },
          { Particle::IncomingBeam2, PDG::Proton },
          { Particle::Parton1, kIntermediateParts[0] },
          { Particle::Parton2, kIntermediateParts[1] }
        },
        { // outgoing state
          { Particle::OutgoingBeam1, { PDG::Proton } },
          { Particle::OutgoingBeam2, { PDG::Proton } },
          { Particle::CentralSystem, kProducedParts }
        }
      );
      setExtraContent();
    }

    unsigned int
    GenericKTProcess::numDimensions( const Kinematics::Mode& ) const
    {
      return num_dimensions_;
    }

    void
    GenericKTProcess::setKinematics( const Kinematics& kin )
    {
      cuts_ = kin;

      //============================================================================================
      // initialise the "constant" (wrt x) part of the Jacobian
      //============================================================================================

      kt_jacobian_ = 1.;
      num_dimensions_ = 0;
      mapped_variables_.clear();

      //============================================================================================
      // register the incoming partons' variables
      //============================================================================================

      registerVariable( qt1_, Mapping::logarithmic, cuts_.cuts.initial.qt, { 1.e-10, 500. }, "First incoming parton virtuality" );
      registerVariable( qt2_, Mapping::logarithmic, cuts_.cuts.initial.qt, { 1.e-10, 500. }, "Second incoming parton virtuality" );
      registerVariable( phi_qt1_, Mapping::linear, cuts_.cuts.initial.phi_qt, { 0., 2.*M_PI }, "First incoming parton azimuthal angle" );
      registerVariable( phi_qt2_, Mapping::linear, cuts_.cuts.initial.phi_qt, { 0., 2.*M_PI }, "Second incoming parton azimuthal angle" );

      //============================================================================================
      // register all process-dependent variables
      //============================================================================================

      preparePhaseSpace();

      //============================================================================================
      // register the outgoing remnants' variables
      //============================================================================================

      MX_ = MY_ = event_->getOneByRole( Particle::IncomingBeam1 ).mass();
      if ( cuts_.mode == Kinematics::Mode::InelasticElastic || cuts_.mode == Kinematics::Mode::InelasticInelastic )
        registerVariable( MX_, Mapping::square, cuts_.cuts.remnants.mass_single, { 1.07, 1000. }, "Positive z proton remnant mass" );
      if ( cuts_.mode == Kinematics::Mode::ElasticInelastic || cuts_.mode == Kinematics::Mode::InelasticInelastic )
        registerVariable( MY_, Mapping::square, cuts_.cuts.remnants.mass_single, { 1.07, 1000. }, "Negative z proton remnant mass" );

      prepareKinematics();
    }

    double
    GenericKTProcess::computeWeight()
    {
      if ( mapped_variables_.size() == 0 )
        throw CG_FATAL( "GenericKTProcess:weight" )
          << "No variables are mapped with this process!";
      if ( kt_jacobian_ == 0. )
        throw CG_FATAL( "GenericKTProcess:weight" )
          << "Point-independant component of the Jacobian for this "
          << "kt-factorised process is null.\n\t"
          << "Please check the validity of the phase space!";

      //============================================================================================
      // generate and initialise all variables, and auxiliary (x-dependent) part of the Jacobian
      // for this phase space point.
      //============================================================================================

      const double aux_jacobian = generateVariables();
      if ( aux_jacobian <= 0. )
        return 0.;

      //============================================================================================
      // compute the integrand and combine together into a single weight for the phase space point.
      //============================================================================================

      const double integrand = computeKTFactorisedMatrixElement();
      if ( integrand <= 0. )
        return 0.;

      const double weight = ( kt_jacobian_*aux_jacobian ) * integrand;

      CG_DEBUG_LOOP( "GenericKTProcess:weight" )
        << "Jacobian = " << kt_jacobian_ << " * " << aux_jacobian
        << "\n\tIntegrand = " << integrand
        << "\n\tdW = " << weight << ".";

      return weight;
    }

    const double GenericKTProcess::kMinFlux = 1.e-20;

    void
    GenericKTProcess::computeIncomingFluxes( double x1, double q1t2, double x2, double q2t2 )
    {
      flux1_ = flux2_ = 0.;
      switch ( cuts_.mode ) {
        case Kinematics::Mode::ElasticElastic:
          flux1_ = flux( Flux::ElasticBudnev, q1t2, x1 );
          flux2_ = flux( Flux::ElasticBudnev, q2t2, x2 );
          break;
        case Kinematics::Mode::ElasticInelastic:
          flux1_ = flux( Flux::ElasticBudnev, q1t2, x1 );
          flux2_ = flux( Flux::InelasticBudnev, q2t2, x2, cuts_.structure_functions, MY_ );
          break;
        case Kinematics::Mode::InelasticElastic:
          flux1_ = flux( Flux::InelasticBudnev, q1t2, x1, cuts_.structure_functions, MX_ );
          flux2_ = flux( Flux::ElasticBudnev, q2t2, x2 );
          break;
        case Kinematics::Mode::InelasticInelastic:
          flux1_ = flux( Flux::InelasticBudnev, q1t2, x1, cuts_.structure_functions, MX_ );
          flux2_ = flux( Flux::InelasticBudnev, q2t2, x2, cuts_.structure_functions, MY_ );
          break;
        default:
          throw CG_FATAL( "GenericKTProcess" ) << "Invalid kinematics mode selected!";
      }
      flux1_ = std::max( flux1_, kMinFlux );
      flux2_ = std::max( flux2_, kMinFlux );
      CG_DEBUG_LOOP( "GenericKTProcess" ) << "Form factors: " << flux1_ << " / " << flux2_ << ".";
    }

    void
    GenericKTProcess::registerVariable( double& out, const Mapping& type,
                                        const Limits& in, Limits default_limits,
                                        const char* description )
    {
      Limits lim = in;
      out = 0.; // reset the variable
      if ( !in.valid() ) {
        CG_DEBUG( "GenericKTProcess:registerVariable" )
          << description << " could not be retrieved from the user configuration!\n\t"
          << "Setting it to the default value: " << default_limits << ".";
        lim = default_limits;
      }
      if ( type == Mapping::logarithmic )
        lim = {
          std::max( log( lim.min() ), -10. ),
          std::min( log( lim.max() ), +10. )
        };
      mapped_variables_.emplace_back( MappingVariable{ description, lim, out, type, num_dimensions_++ } );
      switch ( type ) {
        case Mapping::square:
          kt_jacobian_ *= 2.*lim.range();
          break;
        default:
          kt_jacobian_ *= lim.range();
          break;
      }
      CG_DEBUG( "GenericKTProcess:registerVariable" )
        << description << " has been mapped to variable " << num_dimensions_ << ".\n\t"
        << "Allowed range for integration: " << lim << ".\n\t"
        << "Variable integration mode: " << type << ".";
    }

    void
    GenericKTProcess::dumpVariables() const
    {
      std::ostringstream os;
      for ( const auto& var : mapped_variables_ )
        os << "\n\t(" << var.index << ") " << var.type << " mapping (" << var.description << ") in range " << var.limits;
      CG_INFO( "GenericKTProcess:dumpVariables" )
        << "List of variables handled by this kt-factorised process:"
        << os.str();
    }

    double
    GenericKTProcess::generateVariables() const
    {
      double jacobian = 1.;
      for ( const auto& cut : mapped_variables_ ) {
        if ( !cut.limits.valid() )
          continue;
        const double xv = x( cut.index ); // between 0 and 1
        switch ( cut.type ) {
          case Mapping::linear: {
            cut.variable = cut.limits.x( xv );
          } break;
          case Mapping::logarithmic: {
            cut.variable = exp( cut.limits.x( xv ) );
            jacobian *= cut.variable;
          } break;
          case Mapping::square: {
            cut.variable = cut.limits.x( xv );
            jacobian *= cut.variable;
          } break;
        }
      }
      if ( CG_EXCEPT_MATCH( "KtProcess:vars", debugInsideLoop ) ) {
        std::ostringstream oss;
        for ( const auto& cut : mapped_variables_ ) {
          oss << "variable " << cut.index
              << " in range " << std::left << std::setw( 20 ) << cut.limits << std::right
              << " has value " << cut.variable << "\n\t";
        }
        CG_DEBUG_LOOP( "KtProcess:vars" ) << oss.str();
      }
      return jacobian;
    }

    void
    GenericKTProcess::fillKinematics( bool )
    {
      fillCentralParticlesKinematics(); // process-dependent!
      fillPrimaryParticlesKinematics();
    }

    void
    GenericKTProcess::fillPrimaryParticlesKinematics()
    {
      //============================================================================================
      //     outgoing protons
      //============================================================================================

      Particle& op1 = event_->getOneByRole( Particle::OutgoingBeam1 ),
               &op2 = event_->getOneByRole( Particle::OutgoingBeam2 );

      op1.setMomentum( PX_ );
      op2.setMomentum( PY_ );

      switch ( cuts_.mode ) {
        case Kinematics::Mode::ElasticElastic:
          op1.setStatus( Particle::FinalState );
          op2.setStatus( Particle::FinalState );
          break;
        case Kinematics::Mode::ElasticInelastic:
          op1.setStatus( Particle::FinalState );
          op2.setStatus( Particle::Unfragmented ); op2.setMass( MY_ );
          break;
        case Kinematics::Mode::InelasticElastic:
          op1.setStatus( Particle::Unfragmented ); op1.setMass( MX_ );
          op2.setStatus( Particle::FinalState );
          break;
        case Kinematics::Mode::InelasticInelastic:
          op1.setStatus( Particle::Unfragmented ); op1.setMass( MX_ );
          op2.setStatus( Particle::Unfragmented ); op2.setMass( MY_ );
          break;
        default: {
          throw CG_FATAL( "GenericKTProcess" )
            << "This kT factorisation process is intended for p-on-p collisions! Aborting.";
        } break;
      }

      //============================================================================================
      //     incoming partons (photons, pomerons, ...)
      //============================================================================================

      Particle& g1 = event_->getOneByRole( Particle::Parton1 );
      g1.setMomentum( event_->getOneByRole( Particle::IncomingBeam1 ).momentum()-PX_, true );
      g1.setStatus( Particle::Incoming );

      Particle& g2 = event_->getOneByRole( Particle::Parton2 );
      g2.setMomentum( event_->getOneByRole( Particle::IncomingBeam2 ).momentum()-PY_, true );
      g2.setStatus( Particle::Incoming );

      //============================================================================================
      //     two-parton system
      //============================================================================================

      event_->getOneByRole( Particle::Intermediate ).setMomentum( g1.momentum()+g2.momentum() );
    }

    double
    GenericKTProcess::flux( const Flux& type, double kt2, double x, const StructureFunctions::Type& sf, double mx )
    {
      switch ( type ) {
        case Flux::ElasticBudnev: {
          const double Q2_min = x*x*mp2_/( 1.-x ), Q2_ela = ( kt2+x*x*mp2_ )/( 1.-x );
          const FormFactors ff = FormFactors::ProtonElastic( Q2_ela );
          const double ela1 = ( 1.-x )*( 1.-Q2_min/Q2_ela );
          //const double ela3 = 1.-( Q2_ela-kt2 )/Q2_ela;
          const double f_ela = Constants::alphaEM*M_1_PI/kt2*( ela1*ff.FE + 0.5*x*x*ff.FM );
          // last factor below the Jacobian from dQ^2/Q^2 --> dkT^2/kT^2*(kT^2/Q^2)
          return f_ela*( 1.-x )*kt2 / Q2_ela;
        }
        case Flux::InelasticBudnev: {
          const double mx2 = mx*mx;
          // F2 structure function
          const double Q2min = 1. / ( 1.-x )*( x*( mx2-mp2_ ) + x*x*mp2_ ),
                       Q2 = Q2min + kt2/( 1.-x );
          float xbj = Q2 / ( Q2 + mx2 - mp2_ );
          const StructureFunctions str_fun = StructureFunctionsBuilder::get( sf, Q2, xbj );
          const double term1 = ( 1.-x )*( 1.-Q2min/Q2 );
          const double f_D = str_fun.F2/( mx2 + Q2 - mp2_ ) * term1;
          const double f_C = str_fun.F1( Q2, xbj ) * 2./Q2;
          return Constants::alphaEM*M_1_PI*( 1.-x )/Q2*( f_D+0.5*x*x*f_C );
        }
        case Flux::GluonKMR: {
#ifdef KMR_FLUX
          double logx = log10( x ), logq2 = log10( kt2 ), logmu2 = 2.*log10( mx ), fg = 0.;
          int iread = 1;
          f_inter_kmr_fg_( logx, logq2, logmu2, iread, fg );
          return fg;
#else
          throw CG_FATAL("GenericKTProcess:flux") << "KMR gluon fluxes are not linked to this instance!";
#endif
        }
        default:
          throw CG_FATAL("GenericKTProcess:flux") << "Invalid flux type: " << type;
      }
      return 0.;
    }

    double
    GenericKTProcess::flux( const Flux& type, double kt2, double x, const Kinematics::HeavyIon& hi )
    {
      switch ( type ) {
        case Flux::HIElastic: {
          const double r_a = 1.1*std::pow( hi.A, 1./3 ), a0 = 0.7, m_a = hi.A*mp_;
          const double q2_ela = ( kt2+x*x*m_a*m_a )/( 1.-x ), cons = sqrt( q2_ela )/0.1973;
          const double tau = cons*r_a, tau1 = cons*a0;
          // "Realistic nuclear form-factor" as used in STARLIGHT
          const double ff1 = 3.*( sin( tau )-tau*cos( tau ) )/pow( tau+1.e-10, 3 ), ff2 = 1./( 1.+tau1*tau1 );
          const double ela1 = pow( kt2/( kt2+x*x*m_a*m_a ), 2 ), ela2 = pow( ff1*ff2, 2 )/*, ela3 = 1.-( q2_ela-kt2 )/q2_ela*/;
          return hi.Z*hi.Z*Constants::alphaEM*M_1_PI*ela1*ela2/q2_ela;
        }
        default:
          throw CG_FATAL("GenericKTProcess:flux") << "Invalid flux type: " << type;
      }
      return 0.;
    }

    std::ostream&
    operator<<( std::ostream& os, const GenericKTProcess::Flux& type )
    {
      switch ( type ) {
        case GenericKTProcess::Flux::Elastic: return os << "[p] Elastic";
        case GenericKTProcess::Flux::Inelastic: return os << "[p] Inelastic";
        case GenericKTProcess::Flux::ElasticBudnev: return os << "[p] Budnev elastic";
        case GenericKTProcess::Flux::InelasticBudnev: return os << "[p] Budnev inelastic";
        case GenericKTProcess::Flux::GluonKMR: return os << "[p] Gluon (KMR)";
        case GenericKTProcess::Flux::HIElastic: return os << "[HI] Elastic";
      }
      return os;
    }

    std::ostream&
    operator<<( std::ostream& os, const GenericKTProcess::Mapping& type )
    {
      switch ( type ) {
        case GenericKTProcess::Mapping::linear: return os << "linear";
        case GenericKTProcess::Mapping::logarithmic: return os << "logarithmic";
        case GenericKTProcess::Mapping::square: return os << "squared";
      }
      return os;
    }
  }
}
