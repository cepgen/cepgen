/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Modes.h"

namespace cepgen {
  /// Beam/primary particle's kinematics
  class IncomingBeams : public SteeredObject<IncomingBeams> {
  public:
    explicit IncomingBeams(const ParametersList&);

    static ParametersDescription description();
    void setParameters(const ParametersList&) override;
    const ParametersList& parameters() const override;  ///< List containing all parameters handled

    static mode::Kinematics modeFromBeams(const Beam&, const Beam&);  ///< Extract kinematics type from both beams
    mode::Kinematics mode() const;  ///< Type of kinematics to consider for the phase space

    inline const Beam& positive() const { return pos_beam_; }  ///< Reference to the positive-z beam information
    inline Beam& positive() { return pos_beam_; }              ///< Reference to the positive-z beam information
    inline const Beam& negative() const { return neg_beam_; }  ///< Reference to the negative-z beam information
    inline Beam& negative() { return neg_beam_; }              ///< Reference to the negative-z beam information

    inline const ParametersList& formFactors() const { return formfac_; }        ///< Form factors parameters
    inline const ParametersList& structureFunctions() const { return strfun_; }  ///< Structure functions parameters
    void setStructureFunctions(int, int);  ///< Set the integer-type of structure functions evaluator to build

    void setSqrtS(double);  ///< Set the incoming beams centre of mass energy (in GeV)
    double s() const;       ///< Incoming beams squared centre of mass energy (in GeV^2)
    double sqrtS() const;   ///< Incoming beams centre of mass energy (in GeV)

  private:
    ParametersList formfac_;
    ParametersList strfun_;
    Beam pos_beam_{ParametersList()};
    Beam neg_beam_{ParametersList()};
  };
}  // namespace cepgen

#endif
