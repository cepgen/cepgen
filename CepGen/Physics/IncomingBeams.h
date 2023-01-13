/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#ifndef CepGen_Physics_IncomingBeams_h
#define CepGen_Physics_IncomingBeams_h

#include <iosfwd>
#include <memory>
#include <vector>

#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Modes.h"

namespace cepgen {
  /// Beam/primary particle's kinematics
  class IncomingBeams : public SteeredObject<IncomingBeams> {
  public:
    explicit IncomingBeams(const ParametersList&);

    static ParametersDescription description();
    void setParameters(const ParametersList&) override;
    /// List containing all parameters handled
    const ParametersList& parameters() const override;

    const Beam& positive() const { return pos_beam_; }
    Beam& positive() { return pos_beam_; }
    const Beam& negative() const { return neg_beam_; }
    Beam& negative() { return neg_beam_; }

    /// Type of kinematics to consider for the phase space
    mode::Kinematics mode() const;

    /// Form factors evaluator parameters
    ParametersList formFactors() const { return steer<ParametersList>("formFactors"); }
    /// Structure functions evaluator parameters
    ParametersList structureFunctions() const { return steer<ParametersList>("structureFunctions"); }
    /// Set the integer-type of structure functions evaluator to build
    void setStructureFunctions(int, int);

    /// Set the incoming beams centre of mass energy (in GeV)
    void setSqrtS(double);
    /// Incoming beams squared centre of mass energy (in GeV^2)
    double s() const;
    /// Incoming beams centre of mass energy (in GeV)
    double sqrtS() const;

  private:
    Beam pos_beam_;
    Beam neg_beam_;
  };
}  // namespace cepgen

#endif
