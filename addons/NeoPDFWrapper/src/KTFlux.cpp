/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025-2026  Laurent Forthomme
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

#include <CepGen/Core/Exception.h>
#include <CepGen/Modules/PartonFluxFactory.h>
#include <CepGen/PartonFluxes/KTFlux.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/Limits.h>

#include <NeoPDF.hpp>
#include <cmath>

using namespace std::string_literals;

namespace cepgen::neopdf {
  /// Extract limits in kT for any NeoPDF set
  static Limits ktLimits(const ::neopdf::NeoPDF& pdf) {
    if (const auto range = pdf.param_range(NEOPDF_SUBGRID_PARAMS_KT); range.size() == 2 && range.at(1) > range.at(0))
      return Limits{range.at(0), range.at(1)};
    throw CG_FATAL("neopdf:KTFlux") << "NeoPDF set has no kt dimension; use a collinear flux instead.";
  }
  /// kt-dependent parton flux evaluated from a NeoPDF TMD grid
  class KTFlux final : public cepgen::KTFlux {
  public:
    explicit KTFlux(const ParametersList& params)
        : cepgen::KTFlux(params), parton_pdg_id_(steer<int>("partonPdgId")), fixed_scale_(steer<double>("fixedScale")) {
      const auto set = steer<std::string>("set");
      if (set.empty())
        throw CG_FATAL("neopdf:KTFlux") << "A NeoPDF set name must be provided through the 'name' parameter.";
      const auto member = steer<int>("member");
      if (member < 0)
        throw CG_FATAL("neopdf:KTFlux") << "Invalid member index: " << member << ".";
      neopdf_ = std::make_unique<::neopdf::NeoPDF>(set, static_cast<size_t>(member));
      x_limits_ = Limits{neopdf_->x_min(), neopdf_->x_max()};
      q2_limits_ = Limits{neopdf_->q2_min(), neopdf_->q2_max()};
      kt_limits_ = ktLimits(*neopdf_);
      if (fixed_scale_ > 0. && !q2_limits_.contains(fixed_scale_))
        throw CG_FATAL("neopdf:KTFlux") << "Fixed factorisation scale mu^2=" << fixed_scale_
                                        << " GeV^2 is outside the grid validity range " << q2_limits_ << " GeV^2.";
      CG_INFO("neopdf:KTFlux") << "NeoPDF kt-dependent flux for parton " << static_cast<PDG::Id>(parton_pdg_id_)
                               << " initialised.\n\t"
                               << "Set: " << set << ", member: " << member << ".\n\t"
                               << "x range: " << x_limits_ << ", kt range: " << kt_limits_ << " GeV, "
                               << "mu^2 range: " << q2_limits_ << " GeV^2.\n\t"
                               << "Scale: "
                               << (fixed_scale_ > 0. ? std::to_string(fixed_scale_) + " GeV^2 (fixed)"
                                                     : "parton virtuality"s)
                               << ".";
    }

    static ParametersDescription description() {
      auto desc = cepgen::KTFlux::description();
      desc.setDescription("NeoPDF kt-dependent flux");
      desc.add("set", ""s).setDescription("NeoPDF set name (TMD grid with a kt dimension)");
      desc.add("member", 0).setDescription("set member");
      desc.addAs<int>("partonPdgId", PDG::gluon).setDescription("parton PDG identifier (must be present in the set)");
      desc.add("fixedScale", -1.)
          .setDescription("fixed factorisation scale mu^2 (in GeV^2); (<=0 -> use parton virtuality");
      return desc;
    }

    bool fragmenting() const override { return true; }
    double mass2() const override { return mp2_; }
    spdgid_t partonPdgId() const override { return parton_pdg_id_; }

    double fluxQ2(double x, double kt2, double q2) const override {
      const auto kt = std::sqrt(kt2), mu2 = fixed_scale_ > 0. ? fixed_scale_ : q2;
      if (!x_limits_.contains(x) || !kt_limits_.contains(kt) || !q2_limits_.contains(mu2))
        return 0.;  // outside grid
      return neopdf_->xfxQ2_ND(parton_pdg_id_, {kt, x, mu2 /*TMDlib convention*/});
    }

  private:
    const spdgid_t parton_pdg_id_;
    const double fixed_scale_;
    std::unique_ptr<::neopdf::NeoPDF> neopdf_;
    Limits x_limits_, q2_limits_, kt_limits_;
  };
}  // namespace cepgen::neopdf
using NeoPDFKTFlux = cepgen::neopdf::KTFlux;
REGISTER_KT_FLUX("NeoPDF", 200, NeoPDFKTFlux);
