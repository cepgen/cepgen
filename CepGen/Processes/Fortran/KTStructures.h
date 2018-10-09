#ifndef CepGen_Processes_Fortran_KTStructures_h
#define CepGen_Processes_Fortran_KTStructures_h

namespace cepgen
{
  /// Collection of common blocks for Fortran kT-processes
  namespace ktblock
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
    struct KTKinematics {
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
    struct Cuts {
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
    struct Event {
      int nout; ///< Number of particles in central system
      int pdg[10]; ///< PDG ids of all particles in central system
      double pc[10][4]; ///< 4-momenta of all particles in central system
      double px[4]; ///< 4-momentum of first outgoing proton state
      double py[4]; ///< 4-momentum of second outgoing proton state
    };
  }
}

#endif
