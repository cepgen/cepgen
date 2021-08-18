/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

#ifdef __cplusplus
extern "C" {
#endif
/// Expose structure functions calculators to Fortran
void cepgen_structure_functions_(int& sfmode, double& xbj, double& q2, double& f2, double& fl) {
  using namespace cepgen;
  static auto sf = strfun::StructureFunctionsFactory::get().build(sfmode);
  const auto& val = (*sf)(xbj, q2);
  f2 = val.F2;
  fl = val.FL;
}

/// Compute a \f$k_{\rm T}\f$-dependent flux for single nucleons
/// \param[in] fmode Flux mode (see cepgen::KTFlux)
/// \param[in] x Fractional momentum loss
/// \param[in] kt2 The \f$k_{\rm T}\f$ transverse momentum norm
/// \param[in] sfmode Structure functions set for dissociative emission
/// \param[in] min Incoming particle mass
/// \param[in] mout Diffractive state mass for dissociative emission
double cepgen_kt_flux_(int& fmode, double& x, double& kt2, int& sfmode, double& min, double& mout) {
  using namespace cepgen;
  static auto ff = formfac::FormFactorsFactory::get().build(
      formfac::gFFStandardDipoleHandler);  // use another argument for the modelling?
  static auto sf = strfun::StructureFunctionsFactory::get().build(sfmode);
  ff->setStructureFunctions(sf.get());
  return ktFlux((KTFlux)fmode, x, kt2, *ff, min * min, mout * mout);
}

/// Compute a \f$k_{\rm T}\f$-dependent flux for heavy ions
/// \param[in] fmode Flux mode (see cepgen::KTFlux)
/// \param[in] x Fractional momentum loss
/// \param[in] kt2 The \f$k_{\rm T}\f$ transverse momentum norm
/// \param[in] a Mass number for the heavy ion
/// \param[in] z Atomic number for the heavy ion
double cepgen_kt_flux_hi_(int& fmode, double& x, double& kt2, int& a, int& z) {
  using namespace cepgen;
  return ktFlux((KTFlux)fmode, x, kt2, HeavyIon{(unsigned short)a, (Element)z});
}

/// Mass of a particle, in GeV/c^2
double cepgen_particle_mass_(int& pdg_id) {
  try {
    return cepgen::PDG::get().mass((cepgen::pdgid_t)pdg_id);
  } catch (const cepgen::Exception& e) {
    e.dump();
    exit(0);
  }
}

/// Charge of a particle, in e
double cepgen_particle_charge_(int& pdg_id) {
  try {
    return cepgen::PDG::get().charge((cepgen::pdgid_t)pdg_id);
  } catch (const cepgen::Exception& e) {
    e.dump();
    exit(0);
  }
}

/// Colour factor of a particle
double cepgen_particle_colour_(int& pdg_id) {
  try {
    return cepgen::PDG::get().colours((cepgen::pdgid_t)pdg_id);
  } catch (const cepgen::Exception& e) {
    e.dump();
    exit(0);
  }
}

void cepgen_init_() { cepgen::Generator gen; }

void cepgen_debug_(char* str, int size) { CG_DEBUG("fortran_process") << std::string(str, size); }

void cepgen_warning_(char* str, int size) { CG_WARNING("fortran_process") << std::string(str, size); }

void cepgen_error_(char* str, int size) { CG_ERROR("fortran_process") << std::string(str, size); }

void cepgen_fatal_(char* str, int size) { throw CG_FATAL("fortran_process") << std::string(str, size); }

#ifdef __cplusplus
}
#endif
