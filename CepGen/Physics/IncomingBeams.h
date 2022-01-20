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

    /// Set the incoming beams centre of mass energy (in GeV)
    void setSqrtS(double);
    /// Incoming beams squared centre of mass energy (in GeV^2)
    double s() const;
    /// Incoming beams centre of mass energy (in GeV)
    double sqrtS() const;

    /// Form factors evaluator
    formfac::Parameterisation* formFactors() const { return form_factors_.get(); }
    /// Set a form factors evaluator object
    void setFormFactors(std::unique_ptr<formfac::Parameterisation>);

    /// Structure functions evaluator
    const std::vector<std::shared_ptr<strfun::Parameterisation> >& structureFunctions() const { return str_funs_; }
    /// Set a structure functions evaluator object
    void addStructureFunctions(std::unique_ptr<strfun::Parameterisation>);
    /// Set the integer-type of structure functions evaluator to build
    void addStructureFunctions(int, int);

  private:
    Beam pos_beam_;
    Beam neg_beam_;
    /// Type of form factors to consider
    std::shared_ptr<formfac::Parameterisation> form_factors_;
    /// Type of structure functions to consider
    std::vector<std::shared_ptr<strfun::Parameterisation> > str_funs_;
  };
}  // namespace cepgen

#endif
