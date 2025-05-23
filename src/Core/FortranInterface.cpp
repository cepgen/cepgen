/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2025  Laurent Forthomme
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
#include "CepGen/Generator.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/PartonFluxes/KTFlux.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

using namespace std::string_literals;

static std::unordered_map<int, std::unique_ptr<cepgen::strfun::Parameterisation> > kBuiltStrFunctionsParameterisations;
static std::unordered_map<int, std::unique_ptr<cepgen::KTFlux> > kBuiltKtFluxParameterisations;
static std::unordered_map<std::string, std::unique_ptr<cepgen::Coupling> > kBuiltAlphaSParameterisations,
    kBuiltAlphaEMParameterisations;

#ifdef __cplusplus
extern "C" {
#endif
/// Expose structure functions calculators to Fortran
void cepgen_structure_functions_(int& sfmode, double& xbj, double& q2, double& f2, double& fl) {
  using namespace cepgen;
  if (kBuiltStrFunctionsParameterisations.count(sfmode) == 0)
    kBuiltStrFunctionsParameterisations[sfmode] = StructureFunctionsFactory::get().build(sfmode);
  const auto& sf = kBuiltStrFunctionsParameterisations.at(sfmode);
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
  if (kBuiltKtFluxParameterisations.count(fmode) == 0)
    kBuiltKtFluxParameterisations[fmode] = KTFluxFactory::get().build(
        fmode,
        ParametersList()
            .set("mass", min)
            .set("structureFunctions", StructureFunctionsFactory::get().describeParameters(sfmode).parameters())
            .set("formFactors",
                 FormFactorsFactory::get()
                     .describeParameters(formfac::gFFStandardDipoleHandler)  // use another argument for the modelling?
                     .parameters()));
  return kBuiltKtFluxParameterisations.at(fmode)->fluxMX2(x, kt2, mout * mout);
}

/// Compute a \f$k_{\rm T}\f$-dependent flux for heavy ions
/// \param[in] fmode Flux mode
/// \param[in] x Fractional momentum loss
/// \param[in] kt2 The \f$k_{\rm T}\f$ transverse momentum norm
/// \param[in] a Mass number for the heavy ion
/// \param[in] z Atomic number for the heavy ion
double cepgen_kt_flux_hi_(int& fmode, double& x, double& kt2, int& a, int& z) {
  using namespace cepgen;
  if (kBuiltKtFluxParameterisations.count(fmode) == 0)
    kBuiltKtFluxParameterisations[fmode] =
        KTFluxFactory::get().build(fmode,
                                   ParametersList().setAs<pdgid_t, HeavyIon>(
                                       "heavyIon", HeavyIon{static_cast<unsigned short>(a), static_cast<Element>(z)}));
  return kBuiltKtFluxParameterisations.at(fmode)->fluxMX2(x, kt2, 0.);
}

/// Mass of a particle, in GeV/c^2
double cepgen_particle_mass_(int& pdg_id) {
  try {
    return cepgen::PDG::get().mass(static_cast<cepgen::pdgid_t>(pdg_id));
  } catch (const cepgen::Exception& exception) {
    exception.dump();
    std::exit(EXIT_FAILURE);
  }
}

/// Charge of a particle, in e
double cepgen_particle_charge_(int& pdg_id) {
  try {
    return cepgen::PDG::get().charge(static_cast<cepgen::pdgid_t>(pdg_id));
  } catch (const cepgen::Exception& exception) {
    exception.dump();
    std::exit(EXIT_FAILURE);
  }
}

/// Colour factor of a particle
double cepgen_particle_colour_(int& pdg_id) {
  try {
    return cepgen::PDG::get().colours(static_cast<cepgen::pdgid_t>(pdg_id));
  } catch (const cepgen::Exception& exception) {
    exception.dump();
    std::exit(EXIT_FAILURE);
  }
}

void cepgen_init_() { cepgen::Generator gen; }

void cepgen_debug_(char* str, int size) { CG_DEBUG("fortran_process") << std::string(str, size); }

void cepgen_warning_(char* str, int size) { CG_WARNING("fortran_process") << std::string(str, size); }

void cepgen_error_(char* str, int size) { CG_ERROR("fortran_process") << std::string(str, size); }

void cepgen_fatal_(char* str, int size) { throw CG_FATAL("fortran_process") << std::string(str, size); }

double cepgen_alphas_(char* str, double& q, int size) {
  const auto name = std::string(str, size);
  if (name.empty())
    throw CG_FATAL("cepgen_alphas"s) << "Invalid name retrieved for alphaS(q) parameterisation: '" << name << "'.";
  if (kBuiltAlphaSParameterisations.count(name) == 0)
    kBuiltAlphaSParameterisations[name] = cepgen::AlphaSFactory::get().build(name);
  return kBuiltAlphaSParameterisations.at(name)->operator()(q);
}

double cepgen_alphaem_(char* str, double& q, int size) {
  const auto name = std::string(str, size);
  if (name.empty())
    throw CG_FATAL("cepgen_alphaem"s) << "Invalid name retrieved for alphaEM(q) parameterisation: '" << name << "'.";
  if (kBuiltAlphaEMParameterisations.count(name) == 0)
    kBuiltAlphaEMParameterisations[name] = cepgen::AlphaEMFactory::get().build(name);
  return kBuiltAlphaEMParameterisations.at(name)->operator()(q);
}

#ifdef __cplusplus
}
#endif
