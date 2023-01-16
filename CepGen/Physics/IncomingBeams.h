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

    /// Initialise the beam parameterisation objects
    void initialise();

    /// Find the type of kinematics from the positive/negative beams
    static mode::Kinematics modeFromBeams(const Beam::Mode&, const Beam::Mode&);
    /// Type of kinematics to consider for the phase space
    mode::Kinematics mode() const;

    /// Constant reference to the positive-z beam information
    const Beam& positive() const { return pos_beam_; }
    /// Reference to the positive-z beam information
    Beam& positive() { return pos_beam_; }
    /// Constant reference to the negative-z beam information
    const Beam& negative() const { return neg_beam_; }
    /// Reference to the negative-z beam information
    Beam& negative() { return neg_beam_; }

    /// Form factors evaluator parameters
    const ParametersList& formFactors() const { return formfac_; }
    /// Structure functions evaluator parameters
    const ParametersList& structureFunctions() const { return strfun_; }
    /// Set the integer-type of structure functions evaluator to build
    void setStructureFunctions(int, int);

    /// Set the incoming beams centre of mass energy (in GeV)
    void setSqrtS(double);
    /// Incoming beams squared centre of mass energy (in GeV^2)
    double s() const;
    /// Incoming beams centre of mass energy (in GeV)
    double sqrtS() const;

  private:
    ParametersList formfac_;
    ParametersList strfun_;
    Beam pos_beam_{ParametersList()};
    Beam neg_beam_{ParametersList()};
  };
}  // namespace cepgen

#endif
