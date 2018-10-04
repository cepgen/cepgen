#include "CepGen/Processes/FortranKTProcess.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

extern "C"
{
  /// General physics constants
  struct Constants {
    double m_p; ///< Proton mass
    double units; ///< Conversion factor GeV\f${}^2\to\f$ barn
    double pi; ///< \f$\pi\f$
    double alpha_em; ///< Electromagnetic coupling constant
  };
  /// Generic run parameters
  struct Parameters {
    int icontri; ///< Kinematics mode
    int iflux1; ///< Type of kT-factorised flux for first incoming parton
    int iflux2; ///< Type of kT-factorised flux for second incoming parton
    int imethod; ///< Computation method for matrix element
    int sfmod; ///< Structure functions modelling
    int pdg_l; ///< Central system PDG id
    int a_nuc1; ///< First beam mass number
    int z_nuc1; ///< First beam atomic number
    int a_nuc2; ///< Second beam mass number
    int z_nuc2; ///< Second beam atomic number
    double inp1; ///< First beam momentum, in GeV/c
    double inp2; ///< Second beam momentum, in GeV/c
  };
  /// Kinematics properties of the kT-factorised process
  struct KtKinematics {
   double q1t; ///< Transverse momentum of the first incoming parton
   double q2t; ///< Transverse momentum of the second incoming parton
   double phiq1t; ///< Azimutal angle of the first incoming parton
   double phiq2t; ///< Azimutal angle of the second incoming parton
   double y1; ///< First incoming parton rapidity
   double y2; ///< Second incoming parton rapidity
   double ptdiff; ///< Central system pT balance
   double phiptdiff; ///< Central system azimutal angle difference
   double m_x; ///< Invariant mass for the first diffractive state
   double m_y; ///< Invariant mass for the second diffractive state
  };
  /// Phase space cuts for event kinematics
  struct KinematicsCuts {
    int ipt; ///< Switch for cut on single particle transverse momentum
    int iene; ///< Switch for cut on single particle energy
    int ieta; ///< Switch for cut on single particle pseudo-rapidity
    int iinvm; ///< Switch for cut on central system invariant mass
    int iptsum; ///< Switch for cut on central system transverse momentum
    int idely; ///< Switch for cut on rapididty difference
    double pt_min; ///< Minimal single particle transverse momentum
    double pt_max; ///< Maximal single particle transverse momentum
    double ene_min; ///< Minimal single particle energy
    double ene_max; ///< Maximal single particle energy
    double eta_min; ///< Minimal single particle pseudo-rapidity
    double eta_max; ///< Maximal single particle pseudo-rapidity
    double invm_min; ///< Minimal central system invariant mass
    double invm_max; ///< Maximal central system invariant mass
    double ptsum_min; ///< Minimal central system transverse momentum
    double ptsum_max; ///< Maximal central system transverse momentum
    double dely_min; ///< Minimal rapidity difference for central system
    double dely_max; ///< Maximal rapidity difference for central system
  };
  /// Single event kinematics
  struct EventKinematics {
    double px[4]; ///< 4-momentum of first outgoing proton state
    double py[4]; ///< 4-momentum of second outgoing proton state
    int nout; ///< Number of particles in central system
    int idum; ///< Placeholder for blocks alignment
    int pdg[4]; ///< PDG ids of all particles in central system
    double pc[4][4]; ///< 4-momenta of all particles in central system
  };

  extern Constants constants_;
  extern Parameters params_;
  extern KtKinematics ktkin_;
  extern KinematicsCuts kincuts_;
  extern EventKinematics evtkin_;
}

namespace CepGen
{
  namespace Process
  {
    FortranKTProcess::FortranKTProcess( const ParametersList& params, const char* name, const char* descr, std::function<void( double& )> func ) :
      GenericKTProcess( params, name, descr, { { PDG::photon, PDG::photon } }, { PDG::muon, PDG::muon } ),
      pair_( params.get<int>( "pair", 13 ) ),
      method_( params.get<int>( "method", 1 ) ),
      func_( func )
    {
      constants_.m_p = GenericProcess::mp_;
      constants_.units = Constants::GeV2toBarn;
      constants_.pi = M_PI;
      constants_.alpha_em = Constants::alphaEM;
    }

    void
    FortranKTProcess::preparePhaseSpace()
    {
      mom_ip1_ = event_->getOneByRole( Particle::IncomingBeam1 ).momentum();
      mom_ip2_ = event_->getOneByRole( Particle::IncomingBeam2 ).momentum();

      registerVariable( y1_, Mapping::linear, cuts_.cuts.central.rapidity_single, { -6., 6. }, "First central particle rapidity" );
      registerVariable( y2_, Mapping::linear, cuts_.cuts.central.rapidity_single, { -6., 6. }, "Second central particle rapidity" );
      registerVariable( pt_diff_, Mapping::linear, cuts_.cuts.central.pt_diff, { 0., 50. }, "Transverse momentum difference between central particles" );
      registerVariable( phi_pt_diff_, Mapping::linear, cuts_.cuts.central.phi_pt_diff, { 0., 2.*M_PI }, "Central particles azimuthal angle difference" );

      //===========================================================================================
      // feed phase space cuts to the common block
      //===========================================================================================

      cuts_.cuts.central.pt_single.save( (bool&)kincuts_.ipt, kincuts_.pt_min, kincuts_.pt_max );
      cuts_.cuts.central.energy_single.save( (bool&)kincuts_.iene, kincuts_.ene_min, kincuts_.ene_max );
      cuts_.cuts.central.eta_single.save( (bool&)kincuts_.ieta, kincuts_.eta_min, kincuts_.eta_max );
      cuts_.cuts.central.mass_sum.save( (bool&)kincuts_.iinvm, kincuts_.invm_min, kincuts_.invm_max );
      cuts_.cuts.central.pt_sum.save( (bool&)kincuts_.iptsum, kincuts_.ptsum_min, kincuts_.ptsum_max );
      cuts_.cuts.central.rapidity_diff.save( (bool&)kincuts_.idely, kincuts_.dely_min, kincuts_.dely_max );

      //===========================================================================================
      // feed run parameters to the common block
      //===========================================================================================

      params_.icontri = (int)cuts_.mode;
      params_.imethod = method_;
      params_.sfmod = (int)cuts_.structure_functions->type;
      params_.pdg_l = pair_;

      //-------------------------------------------------------------------------------------------
      // incoming beams information
      //-------------------------------------------------------------------------------------------

      params_.inp1 = cuts_.incoming_beams.first.pz;
      params_.inp2 = cuts_.incoming_beams.second.pz;
      const HeavyIon in1 = (HeavyIon)cuts_.incoming_beams.first.pdg;
      if ( in1 ) {
        params_.a_nuc1 = in1.A;
        params_.z_nuc1 = (unsigned short)in1.Z;
        if ( params_.z_nuc1 > 1 ) {
          event_->getOneByRole( Particle::IncomingBeam1 ).setPdgId( (PDG)in1 );
          event_->getOneByRole( Particle::OutgoingBeam1 ).setPdgId( (PDG)in1 );
        }
      }
      else
        params_.a_nuc1 = params_.z_nuc1 = 1;

      const HeavyIon in2 = (HeavyIon)cuts_.incoming_beams.second.pdg;
      if ( in2 ) {
        params_.a_nuc2 = in2.A;
        params_.z_nuc2 = (unsigned short)in2.Z;
        if ( params_.z_nuc2 > 1 ) {
          event_->getOneByRole( Particle::IncomingBeam2 ).setPdgId( (PDG)in2 );
          event_->getOneByRole( Particle::OutgoingBeam2 ).setPdgId( (PDG)in2 );
        }
      }
      else
        params_.a_nuc2 = params_.z_nuc2 = 1;

      //-------------------------------------------------------------------------------------------
      // intermediate partons information
      //-------------------------------------------------------------------------------------------

      params_.iflux1 = (int)cuts_.incoming_beams.first.kt_flux;
      params_.iflux2 = (int)cuts_.incoming_beams.second.kt_flux;
      if ( (KTFlux)params_.iflux1 == KTFlux::P_Gluon_KMR )
        event_->getOneByRole( Particle::Parton1 ).setPdgId( PDG::gluon );
      if ( (KTFlux)params_.iflux2 == KTFlux::P_Gluon_KMR )
        event_->getOneByRole( Particle::Parton2 ).setPdgId( PDG::gluon );
    }

    double
    FortranKTProcess::computeKTFactorisedMatrixElement()
    {
      ktkin_.q1t = qt1_;
      ktkin_.q2t = qt2_;
      ktkin_.phiq1t = phi_qt1_;
      ktkin_.phiq2t = phi_qt2_;
      ktkin_.y1 = y1_;
      ktkin_.y2 = y2_;
      ktkin_.ptdiff = pt_diff_;
      ktkin_.phiptdiff = phi_pt_diff_;
      ktkin_.m_x = MX_;
      ktkin_.m_y = MY_;

      //--- compute the event weight
      double weight = 0.;
      func_( weight );
      return weight;
    }

    void
    FortranKTProcess::fillCentralParticlesKinematics()
    {
      //===========================================================================================
      // outgoing beam remnants
      //===========================================================================================

      PX_ = Particle::Momentum( evtkin_.px );
      PY_ = Particle::Momentum( evtkin_.py );
      // express these momenta per nucleon
      PX_ *= 1./params_.a_nuc1;
      PY_ *= 1./params_.a_nuc2;

      //===========================================================================================
      // intermediate partons
      //===========================================================================================

      const Particle::Momentum mom_par1 = mom_ip1_-PX_, mom_par2 = mom_ip2_-PY_;
      event_->getOneByRole( Particle::Parton1 ).setMomentum( mom_par1 );
      event_->getOneByRole( Particle::Parton2 ).setMomentum( mom_par2 );
      event_->getOneByRole( Particle::Intermediate ).setMomentum( mom_par1+mom_par2 );

      //===========================================================================================
      // central system
      //===========================================================================================

      Particles& oc = event_->getByRole( Particle::CentralSystem );
      for ( int i = 0; i < evtkin_.nout; ++i ) {
        Particle& p = oc[i];
        p.setPdgId( evtkin_.pdg[i] );
        p.setStatus( Particle::Status::FinalState );
        p.setMomentum( Particle::Momentum( evtkin_.pc[i] ) );
      }
    }
  }
}

