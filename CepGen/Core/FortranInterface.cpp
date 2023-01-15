/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

#ifdef __cplusplus
extern "C" {
#endif
/// Expose structure functions calculators to Fortran
void cepgen_structure_functions_(int& sfmode, double& xbj, double& q2, double& f2, double& fl) {
  using namespace cepgen;
  static auto sf = strfun::StructureFunctionsFactory::get().build(sfmode);
  f2 = sf->F2(xbj, q2);
  fl = sf->FL(xbj, q2);
}

/// Compute a \f$k_{\rm T}\f$-dependent flux for single nucleons
/// \param[in] fmode Flux mode
/// \param[in] x Fractional momentum loss
/// \param[in] kt2 The \f$k_{\rm T}\f$ transverse momentum norm
/// \param[in] sfmode Structure functions set for dissociative emission
/// \param[in] min Incoming particle mass
/// \param[in] mout Diffractive state mass for dissociative emission
double cepgen_kt_flux_(int& fmode, double& x, double& kt2, int& sfmode, double& min, double& mout) {
  using namespace cepgen;
  const auto params =
      ParametersList()
          .set<double>("mass", min)
          .set<ParametersList>("structureFunctions",
                               strfun::StructureFunctionsFactory::get().describeParameters(sfmode).parameters())
          .set<ParametersList>(
              "formFactors",
              FormFactorsFactory::get()
                  .describeParameters(formfac::gFFStandardDipoleHandler)  // use another argument for the modelling?
                  .parameters());
  auto flux_name = [](int mode) -> std::string {
    switch (mode) {
      case 0:
        return "ElasticKT";
      case 10:
        return "BudnevElasticKT";
      case 1:
        return "InelasticKT";
      case 11:
        return "BudnevInelasticKT";
      case 100:
        return "ElasticHeavyIonKT";
      case 20:
        return "KMRElasticGluonKT";
      default:
        throw CG_FATAL("cepgen_kt_flux") << "Invalid flux modelling: " << mode << ".";
    }
  };
  static auto flux = PartonFluxFactory::get().build(flux_name(fmode), params);
  return (*flux)(x, kt2, mout);
}

/// Compute a \f$k_{\rm T}\f$-dependent flux for heavy ions
/// \param[in] fmode Flux mode
/// \param[in] x Fractional momentum loss
/// \param[in] kt2 The \f$k_{\rm T}\f$ transverse momentum norm
/// \param[in] a Mass number for the heavy ion
/// \param[in] z Atomic number for the heavy ion
double cepgen_kt_flux_hi_(int& fmode, double& x, double& kt2, int& a, int& z) {
  using namespace cepgen;
  (void)fmode;
  static auto flux = PartonFluxFactory::get().build(
      "ElasticHeavyIonKT",
      ParametersList().setAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon{(unsigned short)a, (Element)z}));
  return (*flux)(x, kt2, 0.);
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

double cepgen_alphas_(double& q) {
  static std::unique_ptr<cepgen::Coupling> kAlphaSPtr;
  if (!kAlphaSPtr) {
    CG_INFO("fortran_process") << "Initialisation of the alpha(S) evolution algorithm.";
    kAlphaSPtr = cepgen::AlphaSFactory::get().build("pegasus");
  }
  return (*kAlphaSPtr)(q);
}

#ifdef __cplusplus
}
#endif
