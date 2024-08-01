/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

/// Collection of common blocks for Fortran \f$k_{\rm T}\f$-processes
namespace cepgen::ktblock {
  /// General physics constants
  struct Constants {
    double m_p;    ///< Proton mass
    double units;  ///< Conversion factor GeV\f$^2\to\f$ barn
    double pi;     ///< \f$\pi\f$
  };
  /// Generic run parameters
  struct GenParameters {
    int icontri;  ///< Kinematics mode
    int idum;     ///< Dummy padding variable
    int iflux1;   ///< Type of \f$k_{\rm T}\f$-factorised flux for first incoming parton
    int iflux2;   ///< Type of \f$k_{\rm T}\f$-factorised flux for second incoming parton
    int a_nuc1;   ///< First beam mass number
    int z_nuc1;   ///< First beam atomic number
    int a_nuc2;   ///< Second beam mass number
    int z_nuc2;   ///< Second beam atomic number
    double inp1;  ///< First beam momentum, in GeV/c
    double inp2;  ///< Second beam momentum, in GeV/c
  };
  /// Kinematics properties of the \f$k_{\rm T}\f$-factorised process
  struct KTKinematics {
    double q1t;        ///< Transverse momentum of the first incoming parton
    double q2t;        ///< Transverse momentum of the second incoming parton
    double phiq1t;     ///< Azimuthal angle of the first incoming parton
    double phiq2t;     ///< Azimuthal angle of the second incoming parton
    double y1;         ///< First incoming parton rapidity
    double y2;         ///< Second incoming parton rapidity
    double ptdiff;     ///< Central system pT balance
    double phiptdiff;  ///< Central system azimuthal angle difference
    double m_x;        ///< Invariant mass for the first diffractive state
    double m_y;        ///< Invariant mass for the second diffractive state
  };
  /// Phase space cuts for event kinematics
  struct KinCuts {
    int ipt;           ///< Switch for cut on single particle transverse momentum
    int iene;          ///< Switch for cut on single particle energy
    int ieta;          ///< Switch for cut on single particle pseudo-rapidity
    int iinvm;         ///< Switch for cut on central system invariant mass
    int iptsum;        ///< Switch for cut on central system transverse momentum
    int idely;         ///< Switch for cut on rapidity difference
    double pt_min;     ///< Minimal single particle transverse momentum
    double pt_max;     ///< Maximal single particle transverse momentum
    double ene_min;    ///< Minimal single particle energy
    double ene_max;    ///< Maximal single particle energy
    double eta_min;    ///< Minimal single particle pseudo-rapidity
    double eta_max;    ///< Maximal single particle pseudo-rapidity
    double invm_min;   ///< Minimal central system invariant mass
    double invm_max;   ///< Maximal central system invariant mass
    double ptsum_min;  ///< Minimal central system transverse momentum
    double ptsum_max;  ///< Maximal central system transverse momentum
    double dely_min;   ///< Minimal rapidity difference for central system
    double dely_max;   ///< Maximal rapidity difference for central system
  };
  /// Single event kinematics
  struct EventKinematics {
    static constexpr size_t MAX_PART = 10;
    int nout;                ///< Number of particles in central system
    int pdg[MAX_PART];       ///< PDG ids of all particles in central system
    int idum;                ///< Padding
    double pc[MAX_PART][4];  ///< 4-momenta of all particles in central system
    double px[4];            ///< 4-momentum of first outgoing proton state
    double py[4];            ///< 4-momentum of second outgoing proton state
  };
}  // namespace cepgen::ktblock

#endif
