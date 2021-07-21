#include "CepGen/Physics/Kinematics.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  const double Kinematics::MX_MIN = 1.07;  // mp+mpi+-

  Kinematics::Kinematics(const ParametersList& params) {
    //----- per-incoming beam kinematics
    incoming_beams = IncomingBeams(params);

    //----- phase space definition
    setParameters(params);
  }

  void Kinematics::setParameters(const ParametersList& params) {
    //--- initial partons
    cuts.initial.setParameters(params);

    //--- central system
    cuts.central.setParameters(params);
    if (params.has<Limits>("phiptdiff")) {
      CG_WARNING("Kinematics") << "\"phiptdiff\" parameter is deprecated! "
                               << "Please use \"phidiff\" instead.";
      params.fill<Limits>("phiptdiff", cuts.central.phi_diff());  //legacy
    }
    if (params.has<std::vector<int> >("minFinalState"))
      for (const auto& pdg : params.get<std::vector<int> >("minFinalState"))
        minimum_final_state.emplace_back((pdgid_t)pdg);
    if (params.has<ParametersList>("cuts")) {  // per-particle cuts
      const auto& per_parts = params.get<ParametersList>("cuts");
      for (const auto& part : per_parts.keys())
        cuts.central_particles[(pdgid_t)stoi(part)].setParameters(per_parts.get<ParametersList>(part));
    }

    //--- outgoing remnants
    cuts.remnants.setParameters(params);
    // sanity check
    if (cuts.remnants.mx().min() < MX_MIN) {
      CG_WARNING("Kinematics:setParameters") << "Minimum diffractive mass set to " << MX_MIN << " GeV.";
      cuts.remnants.mx().min() = MX_MIN;
    }

    //--- specify where to look for the grid path for gluon emission
    if (params.has<std::string>("kmrGridPath"))
      kmr::GluonGrid::get(params.get<std::string>("kmrGridPath"));
  }

  ParametersList Kinematics::parameters() const {
    ParametersList params;
    params += incoming_beams.parameters();
    for (const auto& lim : cuts.initial.list())
      params.set<Limits>(lim.name, lim.limits);
    for (auto& lim : cuts.central.list())
      params.set<Limits>(lim.name, lim.limits);
    if (!minimum_final_state.empty()) {
      std::vector<int> min_pdgs;
      for (const auto& pdg : minimum_final_state)
        min_pdgs.emplace_back((int)pdg);
      params.set<std::vector<int> >("minFinalState", min_pdgs);
    }
    if (!cuts.central_particles.empty()) {
      ParametersList per_part;
      for (const auto& cuts_vs_part : cuts.central_particles) {
        ParametersList cuts_vs_id;
        for (const auto& lim : cuts_vs_part.second.list())
          params.set<Limits>(lim.name, lim.limits);
        per_part.set<ParametersList>(std::to_string(cuts_vs_part.first), cuts_vs_id);
      }
      params.set<ParametersList>("cuts", per_part);
    }
    for (const auto& lim : cuts.remnants.list())
      params.set<Limits>(lim.name, lim.limits);
    return params;
  }

  std::ostream& operator<<(std::ostream& os, const Kinematics::CutsList& kin) {
    std::string sep;
    os << "initial: {";
    for (const auto& cut : kin.initial.list())
      os << sep << cut, sep = ", ";
    os << "}, central: {";
    sep.clear();
    for (const auto& cut : kin.central.list())
      os << sep << cut, sep = ", ";
    os << "}, remnants: {";
    sep.clear();
    for (const auto& cut : kin.remnants.list())
      os << sep << cut, sep = ", ";
    return os << "}";
  }

  //--------------------------------------------------------------------
  // List of kinematics limits
  //--------------------------------------------------------------------

  Kinematics::CutsList::CutsList()
      : initial(ParametersList().set<Limits>("q2", {0., 1.e5})),
        central(ParametersList().set<double>("ptmin", 0.)),
        remnants(ParametersList().set<Limits>("mx", {MX_MIN, 1000.})) {}
}  // namespace cepgen
