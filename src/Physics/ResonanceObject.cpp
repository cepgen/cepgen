/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ResonanceObject.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Utils/Message.h"

using namespace cepgen;

/// General definition for a resonance
ResonanceObject::ResonanceObject(const ParametersList& params)
    : SteeredObject(params),
      br_(steer<ParametersList>("branchingRatios")),
      ang_mom_(steer<int>("angularMomentum")),
      x0_(steer<double>("x0")),
      mass_(steer<double>("mass")),
      width_(steer<double>("width")),
      mp_(PDG::get().mass(PDG::proton)),
      mp2_(mp_ * mp_),
      mpi2_(std::pow(PDG::get().mass(PDG::piZero), 2)),
      meta2_(std::pow(PDG::get().mass(PDG::eta), 2)),
      x02_(x0_ * x0_) {}

ParametersDescription ResonanceObject::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Set of physical properties for one resonance");
  desc.add<ParametersDescription>("branchingRatios", BranchingRatios::description());
  desc.add<int>("angularMomentum", 0).setDescription("meson angular momentum");
  desc.add<double>("x0", 0.).setDescription("damping parameter");
  desc.add<double>("mass", 0.).setDescription("mass, in GeV/c^2");
  desc.add<double>("width", 0.).setDescription("full width, in GeV");
  return desc;
}

double ResonanceObject::ecmr(double m2) const { return mass_ == 0 ? 0. : utils::energyFromW(mass_, mp2_, m2); }

double ResonanceObject::partialWidth(const KinematicsBlock& kin) const {
  double par_width = 0.;
  if (br_.single_pion > 0.) {
    //----- 1-pion decay mode
    const double pcmrpi = pcmr(mpi2_);
    par_width += br_.single_pion * (std::pow(kin.ppicm / pcmrpi, 2. * ang_mom_ + 1.) *
                                    std::pow((pcmrpi * pcmrpi + x02_) / (kin.ppicm * kin.ppicm + x02_), ang_mom_));
  }
  if (br_.double_pion > 0.) {
    //----- 2-pion decay mode
    const double pcmrpi2 = pcmr(4. * mpi2_);
    par_width +=
        br_.double_pion *
        (std::pow(kin.ppi2cm / pcmrpi2, 2. * (ang_mom_ + 2.)) *
         std::pow((pcmrpi2 * pcmrpi2 + x02_) / (kin.ppi2cm * kin.ppi2cm + x02_), ang_mom_ + 2) * kin.w / mass_);
  }
  if (br_.eta > 0.) {
    //----- eta decay mode
    const double pcmreta = pcmr(meta2_);
    par_width += br_.eta * (std::pow(kin.petacm / pcmreta, 2. * ang_mom_ + 1.) *
                            std::pow((pcmreta * pcmreta + x02_) / (kin.petacm * kin.petacm + x02_), ang_mom_));
  }
  return width_ * par_width;
}

double ResonanceObject::photonWidth(const KinematicsBlock& kin) const {
  const double kcm2 = kin.kcm * kin.kcm;
  const double kcmr2 = std::pow(kcmr(), 2);
  return width_ * kcm2 / kcmr2 * (kcmr2 + x02_) / (kcm2 + x02_);
}

/// kinematics needed for threshold relativistic B-W
ResonanceObject::KinematicsBlock::KinematicsBlock(double w2, double q2, double mp2, double mpi2, double meta2)
    : w2(w2),
      w(std::sqrt(w2)),
      q2(q2),
      k(0.5 * (w2 - mp2) / std::sqrt(mp2)),
      kcm(utils::energyFromW(w, mp2, 0.)),
      ppicm(mom(utils::energyFromW(w, mp2, mpi2), mpi2)),
      ppi2cm(mom(utils::energyFromW(w, mp2, 4. * mpi2), 4. * mpi2)),
      petacm(mom(utils::energyFromW(w, mp2, meta2), meta2)) {}

ResonanceObject::BranchingRatios::BranchingRatios(const ParametersList& params)
    : SteeredObject(params),
      single_pion(steer<double>("singlePi")),
      double_pion(steer<double>("doublePi")),
      eta(steer<double>("eta")) {
  if (!valid())
    CG_WARNING("ResonanceObject:BranchingRatios")
        << "Invalid branching fractions. Sum = " << (single_pion + double_pion + eta) << " != 1.";
}

ParametersDescription ResonanceObject::BranchingRatios::description() {
  auto desc = ParametersDescription();
  desc.add<double>("singlePi", 0.).setDescription("branching fraction for a resonance decay into a single pion");
  desc.add<double>("doublePi", 0.).setDescription("branching fraction for a resonance decay into a pion pair");
  desc.add<double>("eta", 0.).setDescription("branching fraction for a resonance decay into an eta");
  return desc;
}
