/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#ifndef CepGen_Validation_InitialStates_h
#define CepGen_Validation_InitialStates_h

#include <cstdlib>
#include <string>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/HeavyIon.h"

namespace cepgen::validation {
  /// Description of a named set of incoming beams, to be used by the validation (comparison) executables
  /// \note Adding a new initial state is done by appending to InitialState::all().
  struct InitialState {
    std::string name;                 ///< short name, used on the command line and in output file names
    std::string label;                ///< plot label
    std::vector<int> pdg_ids;         ///< beams PDG ids (positive-z beam first)
    std::vector<double> pz;           ///< beams longitudinal momenta (GeV; per nucleon for heavy ions)
    std::vector<std::string> ff;      ///< elastic form factors modelling, per beam
    double ptmin{5.};                 ///< single-lepton transverse momentum threshold (GeV)
    double mass_max{50.};             ///< upper bound for dilepton invariant mass plots (GeV)
    double pt_pair_max{5.};           ///< upper bound for pair transverse momentum plots (GeV)
    double ww_mass_max{500.};         ///< upper bound for W pair invariant mass plots (GeV)

    /// Kinematics parameters (beams) common to all processes
    ParametersList kinematics() const {
      std::vector<ParametersList> form_factors;
      for (const auto& ff_name : ff)
        form_factors.emplace_back(ParametersList().setName(ff_name));
      return ParametersList()
          .set<std::vector<int> >("pdgIds", pdg_ids)
          .set<std::vector<double> >("pz", pz)
          .set<std::vector<ParametersList> >("formFactors", form_factors)
          .set<int>("mode", 1 /* elastic-elastic: coherent/ultra-peripheral photon emission */);
    }

    /// Name of the kt-factorised parton flux modelling for one beam, for a given family of modellings
    /// \param[in] family modelling family ("default" lets CepGen choose; "elastic"; "budnev"; "klein")
    /// \param[in] beam_pdg_id PDG id of the beam
    static std::string ktFlux(const std::string& family, int beam_pdg_id) {
      const bool lepton = std::abs(beam_pdg_id) >= 11 && std::abs(beam_pdg_id) <= 16;
      const bool ion = HeavyIon::isHI(std::abs(beam_pdg_id)) && std::abs(beam_pdg_id) != 2212;
      if (family == "default")
        return "";
      if (lepton)
        return "BudnevElasticLepton";
      if (family == "elastic")
        return ion ? "ElasticHeavyIon" : "Elastic";
      if (family == "budnev")
        return ion ? "BudnevElasticHeavyIon" : "BudnevElastic";
      if (family == "klein")
        return "KleinElasticHI";
      throw CG_FATAL("InitialState:ktFlux") << "Unknown flux modelling family: '" << family << "'. "
                                            << "Allowed values: default, elastic, budnev, klein.";
    }

    /// Kinematics parameters, including a given family of kt-factorised flux modelling
    ParametersList kinematics(const std::string& flux_family) const {
      auto params = kinematics();
      if (flux_family == "default")
        return params;
      std::vector<ParametersList> fluxes;
      for (const auto& pdg_id : pdg_ids)
        fluxes.emplace_back(ParametersList().setName(ktFlux(flux_family, pdg_id)));
      return params.set<std::vector<ParametersList> >("partonFluxes", fluxes);
    }

    /// All known initial states (proton, lead, electron beams)
    static const std::vector<InitialState>& all() {
      // Pb momenta are given per nucleon, i.e. Z/A times the proton momentum for equal-rigidity beams
      static const auto pb = static_cast<int>(HeavyIon::Pb());
      static const std::vector<InitialState> states = {
          {"pp", "pp, 13.6 TeV", {2212, 2212}, {6800., 6800.}, {"StandardDipole", "StandardDipole"}, 10., 500., 50., 1000.},
          {"pPb", "pPb, 8.16 TeV", {2212, pb}, {6500., 2562.5}, {"StandardDipole", "HeavyIonDipole"}, 5., 50., 5., 500.},
          {"PbPb", "PbPb, 5.36 TeV", {pb, pb}, {2680.8, 2680.8}, {"HeavyIonDipole", "HeavyIonDipole"}, 5., 50., 2., 300.},
          {"ep", "ep, 50 GeV $\\times$ 7 TeV", {2212, 11}, {7000., 50.}, {"StandardDipole", "PointLikeFermion"}, 2.5, 50., 5., 300.},
      };
      return states;
    }
    /// Retrieve a state by its name
    static const InitialState& get(const std::string& name) {
      for (const auto& state : all())
        if (state.name == name)
          return state;
      std::string known;
      for (const auto& state : all())
        known += (known.empty() ? "" : ", ") + state.name;
      throw CG_FATAL("InitialState:get") << "Unknown initial state: '" << name << "'. Known states: " << known << ".";
    }
  };
}  // namespace cepgen::validation

#endif
