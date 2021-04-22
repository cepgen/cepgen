#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/GluonGrid.h"

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  const double Kinematics::MX_MIN = 1.07;  // mp+mpi+-

  Kinematics::Kinematics(const ParametersList& params) {
    //----- per-incoming beam kinematics

    incoming_beams.first.pdg = params.get<int>("beam1id", (int)PDG::proton);
    params.fill<double>("beam1pz", incoming_beams.first.pz);
    const int hi_A1 = params.get<int>("beam1A", 1);
    const int hi_Z1 = params.get<int>("beam1Z", 0);
    if (hi_Z1 != 0)
      incoming_beams.first.pdg = HeavyIon(hi_A1, (Element)hi_Z1);
    if (params.has<std::vector<int> >("heavyIonA")) {
      const auto& hi_beam = params.get<std::vector<int> >("heavyIonA");
      if (hi_beam.size() != 2)
        throw CG_FATAL("Kinematics") << "Invalid format for first incoming beam's HI specification!\n\t"
                                     << "A pair of (A,Z) is required, got " << hi_beam << ".";
      incoming_beams.first.pdg = HeavyIon{(unsigned short)hi_beam.at(0), (Element)hi_beam.at(1)};
    }

    incoming_beams.second.pdg = params.get<int>("beam2id", (int)PDG::proton);
    params.fill<double>("beam2pz", incoming_beams.second.pz);
    const int hi_A2 = params.get<int>("beam2A", 1);
    const int hi_Z2 = params.get<int>("beam2Z", 0);
    if (hi_Z2 != 0)
      incoming_beams.second.pdg = HeavyIon(hi_A2, (Element)hi_Z2);
    if (params.has<std::vector<int> >("heavyIonB")) {
      const auto& hi_beam = params.get<std::vector<int> >("heavyIonB");
      if (hi_beam.size() != 2)
        throw CG_FATAL("Kinematics") << "Invalid format for second incoming beam's HI specification!\n\t"
                                     << "A pair of (A,Z) is required, got " << hi_beam << ".";
      incoming_beams.second.pdg = HeavyIon{(unsigned short)hi_beam.at(0), (Element)hi_beam.at(1)};
    }

    //----- combined two-beam system

    //--- beams PDG ids
    if (params.has<std::vector<ParametersList> >("pdgIds")) {
      const auto& beams_pdg = params.get<std::vector<ParametersList> >("pdgIds");
      if (beams_pdg.size() != 2)
        throw CG_FATAL("Kinematics") << "Invalid list of PDG ids retrieved for incoming beams:\n\t"
                                     << "2 PDG ids are expected, " << beams_pdg << " provided.";
      incoming_beams.first.pdg = beams_pdg.at(0).get<int>("pdgid");
      incoming_beams.second.pdg = beams_pdg.at(1).get<int>("pdgid");
    }
    //--- beams longitudinal momentum
    if (params.has<std::vector<double> >("pz")) {
      const auto& beams_pz = params.get<std::vector<double> >("pz");
      if (beams_pz.size() != 2)
        throw CG_FATAL("Kinematics") << "Invalid format for beams pz specification!\n\t"
                                     << "A vector of two pz's is required.";
      incoming_beams.first.pz = beams_pz.at(0);
      incoming_beams.second.pz = beams_pz.at(1);
    }
    //--- centre-of-mass energy
    const double sqrt_s = params.get<double>("sqrtS", params.get<double>("cmEnergy", -1.));
    if (sqrt_s > 0.)
      setSqrtS(sqrt_s);
    //--- form factors
    if (params.has<std::string>("formFactors") || !form_factors_) {
      std::string ff_mod = params.get<std::string>("formFactors");
      if (ff_mod.empty())
        ff_mod = "StandardDipole";
      form_factors_ = formfac::FormFactorsFactory::get().build(ff_mod);
    }
    if (params.get<int>("mode", (int)mode::Kinematics::invalid) != (int)mode::Kinematics::invalid)
      setMode((mode::Kinematics)params.get<int>("mode"));
    //--- structure functions
    auto strfun = params.get<ParametersList>("structureFunctions");
    if (!strfun.empty() || !str_fun_) {
      if (strfun.name<int>(-999) == -999)
        strfun.setName<int>(11);  // default is Suri-Yennie
      str_fun_ = strfun::StructureFunctionsFactory::get().build(strfun);
      if (form_factors_)
        form_factors_->setStructureFunctions(str_fun_.get());
    }
    //--- parton fluxes for kt-factorisation
    if (params.has<std::vector<int> >("ktFluxes")) {
      auto kt_fluxes = params.get<std::vector<int> >("ktFluxes");
      if (!kt_fluxes.empty()) {
        incoming_beams.first.kt_flux = (KTFlux)kt_fluxes.at(0);
        incoming_beams.second.kt_flux = (kt_fluxes.size() > 1) ? (KTFlux)kt_fluxes.at(1) : (KTFlux)kt_fluxes.at(0);
      }
    } else if (params.has<int>("ktFluxes"))
      incoming_beams.first.kt_flux = incoming_beams.second.kt_flux = (KTFlux)params.get<int>("ktFluxes");

    //----- phase space definition

    //--- initial partons
    for (auto& cut : cuts.initial.rawList()) {
      Limits buf;
      params.fill<Limits>(cut.name, cut.limits)
          .fill<double>(cut.name + "min", buf.min())
          .fill<double>(cut.name + "max", buf.max());
      if (buf.valid())
        cut.limits = buf.validate();
    }

    //--- central system
    for (auto& cut : cuts.central.rawList()) {
      Limits buf;
      params.fill<Limits>(cut.name, buf)
          .fill<double>(cut.name + "min", buf.min())
          .fill<double>(cut.name + "max", buf.max());
      if (buf.valid())
        cut.limits = buf.validate();
    }
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
      for (const auto& part : per_parts.keys()) {
        const auto& part_cuts = per_parts.get<ParametersList>(part);
        for (auto& cut : cuts.central_particles[(pdgid_t)stoi(part)].rawList()) {
          Limits buf;
          part_cuts.fill<Limits>(cut.name, buf)
              .fill<double>(cut.name + "min", buf.min())
              .fill<double>(cut.name + "max", buf.max());
          if (buf.valid())
            cut.limits = buf.validate();
        }
      }
    }

    //--- outgoing remnants
    for (auto& cut : cuts.remnants.rawList()) {
      Limits buf;
      params.fill<Limits>(cut.name, buf)
          .fill<double>(cut.name + "min", buf.min())
          .fill<double>(cut.name + "max", buf.max());
      if (buf.valid())
        cut.limits = buf.validate();
    }
    cuts.remnants.mx().min() = std::max(cuts.remnants.mx().min(), MX_MIN);

    //--- specify where to look for the grid path for gluon emission
    if (params.has<std::string>("kmrGridPath"))
      kmr::GluonGrid::get(params.get<std::string>("kmrGridPath"));
  }

  ParametersList Kinematics::parameters() const {
    ParametersList params;
    params.set<ParametersList>("structureFunctions", str_fun_->parameters())
        .set<int>("mode", (int)mode())
        .set<int>("beam1id", incoming_beams.first.pdg)
        .set<double>("beam1pz", incoming_beams.first.pz)
        .set<int>("beam2id", incoming_beams.second.pdg)
        .set<double>("beam2pz", incoming_beams.second.pz)
        .set<std::vector<int> >("ktFluxes", {(int)incoming_beams.first.kt_flux, (int)incoming_beams.second.kt_flux})
        .set<double>("sqrtS", sqrtS());
    const HeavyIon hi1(incoming_beams.first.pdg), hi2(incoming_beams.second.pdg);
    if (hi1)
      params.set<int>("beam1A", hi1.A).set<int>("beam1Z", (int)hi1.Z);
    if (hi2)
      params.set<int>("beam2A", hi2.A).set<int>("beam2Z", (int)hi2.Z);
    for (const auto& lim : cuts.initial.list()) {
      params.set<Limits>(lim.name, lim.limits);
      if (lim.limits.hasMin())
        params.set<double>(lim.name + "min", lim.limits.min());
      if (lim.limits.hasMax())
        params.set<double>(lim.name + "max", lim.limits.max());
    }
    for (auto& lim : cuts.central.list()) {
      params.set<Limits>(lim.name, lim.limits);
      if (lim.limits.hasMin())
        params.set<double>(lim.name + "min", lim.limits.min());
      if (lim.limits.hasMax())
        params.set<double>(lim.name + "max", lim.limits.max());
    }
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
        for (const auto& lim : cuts_vs_part.second.list()) {
          params.set<Limits>(lim.name, lim.limits);
          if (lim.limits.hasMin())
            params.set<double>(lim.name + "min", lim.limits.min());
          if (lim.limits.hasMax())
            params.set<double>(lim.name + "max", lim.limits.max());
        }
        per_part.set<ParametersList>(std::to_string(cuts_vs_part.first), cuts_vs_id);
      }
      params.set<ParametersList>("cuts", per_part);
    }
    for (const auto& lim : cuts.remnants.list()) {
      params.set<Limits>(lim.name, lim.limits);
      if (lim.limits.hasMin())
        params.set<double>(lim.name + "min", lim.limits.min());
      if (lim.limits.hasMax())
        params.set<double>(lim.name + "max", lim.limits.max());
    }
    return params;
  }

  Kinematics& Kinematics::setSqrtS(double sqrts) {
    if (incoming_beams.first.pdg != incoming_beams.second.pdg)
      throw CG_FATAL("Kinematics") << "Trying to set √s with asymmetric beams"
                                   << " (" << incoming_beams.first.pdg << "/" << incoming_beams.second.pdg << ").\n"
                                   << "Please fill incoming beams objects manually!";
    incoming_beams.first.pz = incoming_beams.second.pz = 0.5 * sqrts;
    return *this;
  }

  double Kinematics::sqrtS() const {
    const HeavyIon hi1(incoming_beams.first.pdg), hi2(incoming_beams.second.pdg);
    const double m1 = hi1 ? HeavyIon::mass(hi1) : PDG::get().mass(incoming_beams.first.pdg);
    const double m2 = hi2 ? HeavyIon::mass(hi2) : PDG::get().mass(incoming_beams.second.pdg);
    const auto p1 = Momentum::fromPxPyPzM(0., 0., +incoming_beams.first.pz, m1);
    const auto p2 = Momentum::fromPxPyPzM(0., 0., -incoming_beams.second.pz, m2);
    return (p1 + p2).mass();
  }

  Kinematics& Kinematics::setMode(const mode::Kinematics& mode) {
    switch (mode) {
      case mode::Kinematics::ElasticElastic:
        incoming_beams.first.mode = mode::Beam::ProtonElastic;
        incoming_beams.second.mode = mode::Beam::ProtonElastic;
        break;
      case mode::Kinematics::ElasticInelastic:
        incoming_beams.first.mode = mode::Beam::ProtonElastic;
        incoming_beams.second.mode = mode::Beam::ProtonInelastic;
        break;
      case mode::Kinematics::InelasticElastic:
        incoming_beams.first.mode = mode::Beam::ProtonInelastic;
        incoming_beams.second.mode = mode::Beam::ProtonElastic;
        break;
      case mode::Kinematics::InelasticInelastic:
        incoming_beams.first.mode = mode::Beam::ProtonInelastic;
        incoming_beams.second.mode = mode::Beam::ProtonInelastic;
        break;
      default:
        throw CG_FATAL("Kinematics:mode") << "Unsupported kinematics mode: " << mode << "!";
    }
    return *this;
  }

  mode::Kinematics Kinematics::mode() const {
    switch (incoming_beams.first.mode) {
      case mode::Beam::ProtonElastic: {
        if (incoming_beams.second.mode == mode::Beam::ProtonElastic)
          return mode::Kinematics::ElasticElastic;
        else
          return mode::Kinematics::ElasticInelastic;
      }
      case mode::Beam::ProtonInelastic: {
        if (incoming_beams.second.mode == mode::Beam::ProtonElastic)
          return mode::Kinematics::InelasticElastic;
        else
          return mode::Kinematics::InelasticInelastic;
      }
      default:
        throw CG_FATAL("Kinematics:mode") << "Unsupported kinematics mode for beams with modes:\n\t"
                                          << incoming_beams.first.mode << " / " << incoming_beams.second.mode << "!";
    }
  }

  Kinematics& Kinematics::setStructureFunctions(int sf_model, int sr_model) {
    const unsigned long kLHAPDFCodeDec = 10000000, kLHAPDFPartDec = 1000000;
    sf_model = (sf_model == 0 ? (int)strfun::Type::SuriYennie : sf_model);
    sr_model = (sr_model == 0 ? (int)sigrat::Type::E143 : sr_model);
    auto sf_params = ParametersList().setName<int>(sf_model).set<ParametersList>(
        "sigmaRatio", ParametersList().setName<int>(sr_model));
    if (sf_model / kLHAPDFCodeDec == 1) {  // SF from parton
      const unsigned long icode = sf_model % kLHAPDFCodeDec;
      sf_params.setName<int>((int)strfun::Type::Partonic)
          .set<int>("pdfId", icode % kLHAPDFPartDec)
          .set<int>("mode", icode / kLHAPDFPartDec);  // 0, 1, 2
    }
    return setStructureFunctions(strfun::StructureFunctionsFactory::get().build(sf_params));
  }

  Kinematics& Kinematics::setStructureFunctions(std::unique_ptr<strfun::Parameterisation> param) {
    str_fun_ = std::move(param);
    form_factors_->setStructureFunctions(str_fun_.get());
    return *this;
  }

  //--------------------------------------------------------------------
  // User-friendly display of various members
  //--------------------------------------------------------------------

  std::ostream& operator<<(std::ostream& os, const Kinematics::Beam& beam) {
    if ((HeavyIon)beam.pdg)
      os << (HeavyIon)beam.pdg;
    else
      os << PDG::get().name(beam.pdg);
    os << " (" << beam.pz << " GeV/c), " << beam.mode;
    if (beam.kt_flux != KTFlux::invalid)
      os << " [unint.flux: " << beam.kt_flux << "]";
    return os;
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
  // Default beam parameters
  //--------------------------------------------------------------------

  Kinematics::Beam::Beam() : pz(0.), pdg(PDG::proton), mode(mode::Beam::invalid), kt_flux(KTFlux::invalid) {}

  //--------------------------------------------------------------------
  // List of kinematics limits
  //--------------------------------------------------------------------

  Kinematics::CutsList::CutsList() {
    initial.q2() = {0., 1.e5};
    central.pt_single().min() = 0.;
    remnants.mx() = {MX_MIN, 320.};
  }
}  // namespace cepgen
