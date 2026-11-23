/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2026  Laurent Forthomme
 *                2017-2019  Wolfgang Schaefer
 *                2019       Marta Luszczak
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

#ifndef CepGen_HelicityAmplitudes_HelicityAmplitude_h
#define CepGen_HelicityAmplitudes_HelicityAmplitude_h

#include <complex>

#include "CepGen/HelicityAmplitudes/Kinematics.h"
#include "CepGen/Modules/ModuleFactory.h"
#include "CepGen/Modules/NamedModule.h"

/// Add an amplitudes computation algorithm definition to the factory
#define REGISTER_HELICITY_AMPLITUDES(name, obj)                                           \
  namespace cepgen::amplitudes {                                                          \
    struct BUILDER_NAME(obj) {                                                            \
      BUILDER_NAME(obj)() { HelicityAmplitudesFactory::get().registerModule<obj>(name); } \
    };                                                                                    \
    static const BUILDER_NAME(obj) gNachAmp##obj;                                         \
  }                                                                                       \
  static_assert(true, "")

namespace cepgen::amplitudes {
  /// Amplitudes computational tool, as developed by Nachtmann et al.
  /// \cite Nachtmann:2005en
  class HelicityAmplitudes : public NamedModule<HelicityAmplitudes> {
  public:
    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.setDescription("generic helicity model");
      desc.add("s1", 0.);
      return desc;
    }

    /// Compute the amplitude for a given kinematics and a given set of helicity components
    virtual std::complex<double> amplitude(const Kinematics&, const HelicityState&) const = 0;
    virtual bool hasAAWW() const { return false; }
    virtual bool hasAAZZ() const { return false; }
    virtual bool hasAAFF() const { return false; }
    virtual bool hasAASS() const { return false; }

  protected:
    explicit HelicityAmplitudes(const ParametersList& params) : NamedModule(params), s1_(steer<double>("s1")) {}

    double s1() const { return s1_; }
    double c1() const { return std::sqrt(1. - s1_ * s1_); }

  private:
    const double s1_;
  };
  /// A helicity-dependent amplitudes calculation algorithms factory
  DEFINE_FACTORY(HelicityAmplitudesFactory, HelicityAmplitudes, "Helicity-dependent amplitudes factory");
}  // namespace cepgen::amplitudes

#endif
