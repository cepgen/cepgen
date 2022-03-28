/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Process_Fortran_KTStructures_h
#define CepGen_Process_Fortran_KTStructures_h

namespace cepgen {
  /// Collection of common blocks for Fortran \f$k_{\rm T}\f$-processes
  namespace ktblock {
    /// General physics constants
    struct Constants {
      double m_p{0.};       ///< Proton mass
      double units{0.};     ///< Conversion factor GeV\f$^2\to\f$ barn
      double pi{0.};        ///< \f$\pi\f$
      double alpha_em{0.};  ///< Electromagnetic coupling constant
    };
    /// Generic run parameters
    struct Parameters {
      int icontri{0};   ///< Kinematics mode
      int idum{0};      ///< Dummy padding variable
      int iflux1{0};    ///< Type of \f$k_{\rm T}\f$-factorised flux for first incoming parton
      int iflux2{0};    ///< Type of \f$k_{\rm T}\f$-factorised flux for second incoming parton
      int a_nuc1{0};    ///< First beam mass number
      int z_nuc1{0};    ///< First beam atomic number
      int a_nuc2{0};    ///< Second beam mass number
      int z_nuc2{0};    ///< Second beam atomic number
      double inp1{0.};  ///< First beam momentum, in GeV/c
      double inp2{0.};  ///< Second beam momentum, in GeV/c
    };
    /// Kinematics properties of the \f$k_{\rm T}\f$-factorised process
    struct KTKinematics {
      double q1t{0.};        ///< Transverse momentum of the first incoming parton
      double q2t{0.};        ///< Transverse momentum of the second incoming parton
      double phiq1t{0.};     ///< Azimuthal angle of the first incoming parton
      double phiq2t{0.};     ///< Azimuthal angle of the second incoming parton
      double y1{0.};         ///< First incoming parton rapidity
      double y2{0.};         ///< Second incoming parton rapidity
      double ptdiff{0.};     ///< Central system pT balance
      double phiptdiff{0.};  ///< Central system azimuthal angle difference
      double m_x{0.};        ///< Invariant mass for the first diffractive state
      double m_y{0.};        ///< Invariant mass for the second diffractive state
    };
    /// Phase space cuts for event kinematics
    struct Cuts {
      bool ipt{false};       ///< Switch for cut on single particle transverse momentum
      bool iene{false};      ///< Switch for cut on single particle energy
      bool ieta{false};      ///< Switch for cut on single particle pseudo-rapidity
      bool iinvm{false};     ///< Switch for cut on central system invariant mass
      bool iptsum{false};    ///< Switch for cut on central system transverse momentum
      bool idely{false};     ///< Switch for cut on rapidity difference
      double pt_min{0.};     ///< Minimal single particle transverse momentum
      double pt_max{0.};     ///< Maximal single particle transverse momentum
      double ene_min{0.};    ///< Minimal single particle energy
      double ene_max{0.};    ///< Maximal single particle energy
      double eta_min{0.};    ///< Minimal single particle pseudo-rapidity
      double eta_max{0.};    ///< Maximal single particle pseudo-rapidity
      double invm_min{0.};   ///< Minimal central system invariant mass
      double invm_max{0.};   ///< Maximal central system invariant mass
      double ptsum_min{0.};  ///< Minimal central system transverse momentum
      double ptsum_max{0.};  ///< Maximal central system transverse momentum
      double dely_min{0.};   ///< Minimal rapidity difference for central system
      double dely_max{0.};   ///< Maximal rapidity difference for central system
    };
    /// Single event kinematics
    struct Event {
      static constexpr size_t MAX_PART = 10;
      int nout{0};             ///< Number of particles in central system
      int pdg[MAX_PART];       ///< PDG ids of all particles in central system
      int idum{0};             ///< Padding
      double pc[MAX_PART][4];  ///< 4-momenta of all particles in central system
      double px[4];            ///< 4-momentum of first outgoing proton state
      double py[4];            ///< 4-momentum of second outgoing proton state
    };
  }  // namespace ktblock
}  // namespace cepgen

#endif
